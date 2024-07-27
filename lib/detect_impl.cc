/* -*- c++ -*- */
/*
 * Copyright 2016 Bastille Networks.
 * Copyright 2024 Zhang Maiyun <maz005@ucsd.edu>.
 *   As a project developed during Maiyun's participation in the UCSD SRIP
 *   program, The University of California may have claims to this work.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "detect_impl.h"

#define DEBUG_OFF     0
#define DEBUG_INFO    1
#define DEBUG_VERBOSE 2
#define DEBUG         DEBUG_OFF

#define DUMP_IQ       0

#define OVERLAP_DEFAULT 1
#define OVERLAP_FACTOR  16

namespace gr {
  namespace lora {

    detect::sptr
    detect::make( uint8_t   spreading_factor,
                  float     beta,
                  uint16_t  fft_factor,
                  uint8_t   peak_search_algorithm,
                  uint16_t  peak_search_phase_k,
                  float     fs_bw_ratio,
                  uint16_t  on_length,
                  const std::string &filename_template)
    {
      return gnuradio::get_initial_sptr
        (new detect_impl(spreading_factor, beta, fft_factor, peak_search_algorithm, peak_search_phase_k, fs_bw_ratio, on_length, filename_template));
    }

    /*
     * The private constructor
     */
    detect_impl::detect_impl( uint8_t   spreading_factor,
                              float     beta,
                              uint16_t  fft_factor,
                              uint8_t   peak_search_algorithm,
                              uint16_t  peak_search_phase_k,
                              float     fs_bw_ratio,
                              uint16_t  on_length,
                              const std::string &filename_template)
      : gr::block("detect",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(0, 0, 0)),
        d_sf(spreading_factor),
        d_beta(beta),
        d_fft_size_factor(fft_factor),
        d_peak_search_algorithm(peak_search_algorithm),
        d_peak_search_phase_k(peak_search_phase_k),
        d_on_length(on_length),
        d_filename_template(filename_template)
    {
      assert((d_sf > 5) && (d_sf < 13));
      if (d_sf == 6) assert(!header);
      assert(d_fft_size_factor > 0);
      assert(((int)fs_bw_ratio) == fs_bw_ratio);
      d_p = (int) fs_bw_ratio;

      d_state = S_RESET;

      d_num_symbols = (1 << d_sf);
      d_num_samples = d_p*d_num_symbols;
      d_bin_size = d_fft_size_factor*d_num_symbols;
      d_fft_size = d_fft_size_factor*d_num_samples;
      d_fft = new fft::fft_complex_fwd(d_fft_size, 1);
      d_overlaps = OVERLAP_DEFAULT;
      d_offset = 0;
      d_preamble_drift_max = d_fft_size_factor * (d_ldr ? 2 : 1);

      d_window = fft::window::build(fft::window::WIN_KAISER, d_num_samples, d_beta);

      // Create local chirp tables.  Each table is 2 chirps long to allow memcpying from arbitrary offsets.
      for (int i = 0; i < d_num_samples; i++) {
        double phase = M_PI/d_p*(i-i*i/(float)d_num_samples);
        d_downchirp.push_back(gr_complex(std::polar(1.0, phase)));
        d_upchirp.push_back(gr_complex(std::polar(1.0, -phase)));
      }

      set_history(DEMOD_HISTORY_DEPTH*d_num_samples);  // Sync is 2.25 chirp periods long
    }

    /*
     * Our virtual destructor.
     */
    detect_impl::~detect_impl()
    {
      delete d_fft;
    }

    uint32_t
    detect_impl::argmax_32f(float *fft_result, float *max_val_p)
    {
      float mag   = abs(fft_result[0]);
      float max_val = mag;
      uint32_t   max_idx = 0;

      for (uint32_t i = 0; i < d_bin_size; i++)
      {
        mag = abs(fft_result[i]);
        if (mag > max_val)
        {
          max_idx = i;
          max_val = mag;
        }
      }

      *max_val_p = max_val;
      return max_idx;
    }

    uint32_t
    detect_impl::search_fft_peak(const lv_32fc_t *fft_result,
                                float *buffer1, float *buffer2,
                                gr_complex *buffer_c, float *max_val_p)
    {
      // size of buffer1:   d_fft_size (float)
      // size of buffer2:   d_bin_size  (float)
      // size of buffer_c:  d_bin_size  (complex)
      uint32_t max_idx = 0;
      *max_val_p = 0;
      if (d_peak_search_algorithm == FFT_PEAK_SEARCH_ABS)
      {
        // fft result magnitude summation
        volk_32fc_magnitude_32f(buffer1, fft_result, d_fft_size);
        volk_32f_x2_add_32f(buffer2, buffer1, &buffer1[d_fft_size-d_bin_size], d_bin_size);

        // Take argmax of returned FFT (similar to MFSK demod)
        max_idx = argmax_32f(buffer2, max_val_p);
      }
      else if (d_peak_search_algorithm == FFT_PEAK_SEARCH_PHASE)
      {
        uint32_t tmp_max_idx;
        float tmp_max_val;
        for (int i = 0; i < d_peak_search_phase_k; i++)
        {
          float phase_offset = 2*M_PI/d_peak_search_phase_k*i;
          tmp_max_idx = fft_add(fft_result, buffer2, buffer_c, &tmp_max_val, phase_offset);
          if (tmp_max_val > *max_val_p)
          {
            *max_val_p = tmp_max_val;
            max_idx = tmp_max_idx;
          }
        }
      }
      else
      {
        max_idx = fft_add(fft_result, buffer2, buffer_c, max_val_p, 0);
      }

      return max_idx;
    }

    uint32_t
    detect_impl::fft_add(const lv_32fc_t *fft_result, float *buffer, gr_complex *buffer_c,
                        float *max_val_p, float phase_offset)
    {
      lv_32fc_t s = lv_cmake((float)std::cos(phase_offset), (float)std::sin(phase_offset));
      volk_32fc_s32fc_multiply_32fc(buffer_c, fft_result, s, d_bin_size);
      volk_32fc_x2_add_32fc(buffer_c, buffer_c, &fft_result[d_fft_size-d_bin_size], d_bin_size);
      volk_32fc_magnitude_32f(buffer, buffer_c, d_bin_size);
      return argmax_32f(buffer, max_val_p);
    }

    uint16_t
    detect_impl::argmax(gr_complex *fft_result)
    {
      float magsq   = pow(real(fft_result[0]), 2) + pow(imag(fft_result[0]), 2);
      float max_val = magsq;
      uint16_t   max_idx = 0;


      for (uint16_t i = 0; i < d_fft_size; i++)
      {
        magsq = pow(real(fft_result[i]), 2) + pow(imag(fft_result[i]), 2);
        if (magsq > max_val)
        {
          max_idx = i;
          max_val = magsq;
        }
      }

      return max_idx;
    }

    int
    detect_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      if (ninput_items[0] < DEMOD_HISTORY_DEPTH*d_num_samples) return 0;
      const gr_complex *in0 = (const gr_complex *) input_items[0];
      const gr_complex *in  = &in0[(DEMOD_HISTORY_DEPTH-1)*d_num_samples];

      uint32_t num_consumed   = d_num_samples;
      uint32_t max_idx        = 0;
      uint32_t max_idx_sfd    = 0;
      bool         preamble_found = false;
      bool         sfd_found      = false;
      float        max_val        = 0;
      float        max_val_sfd    = 0;

      // Nomenclature:
      //  up_block   == de-chirping buffer to contain upchirp features: the preamble, sync word, and data chirps
      //  down_block == de-chirping buffer to contain downchirp features: the SFD
      gr_complex *up_block   = (gr_complex *)volk_malloc(d_fft_size*sizeof(gr_complex), volk_get_alignment());
      gr_complex *down_block = (gr_complex *)volk_malloc(d_fft_size*sizeof(gr_complex), volk_get_alignment());
      float *fft_res_mag = (float*)volk_malloc(d_fft_size*sizeof(float), volk_get_alignment());
      float *fft_res_add = (float*)volk_malloc(d_bin_size*sizeof(float), volk_get_alignment());
      gr_complex *fft_res_add_c = (gr_complex*)volk_malloc(d_bin_size*sizeof(gr_complex), volk_get_alignment());

      if (up_block == NULL || down_block == NULL ||
          fft_res_mag == NULL || fft_res_add == NULL || fft_res_add_c == NULL)
      {
        std::cerr << "Unable to allocate processing buffer!" << std::endl;
      }

      // Dechirp the incoming signal
      volk_32fc_x2_multiply_32fc(up_block, in, &d_downchirp[0], d_num_samples);

      // Enable to write IQ to disk for debugging
      #if DUMP_IQ
        f_up_windowless.write((const char*)&up_block[0], d_num_samples*sizeof(gr_complex));
      #endif

      // Windowing
      // volk_32fc_32f_multiply_32fc(up_block, up_block, &d_window[0], d_num_samples);

      #if DUMP_IQ
        if (d_state != S_SFD_SYNC) f_down.write((const char*)&down_block[0], d_num_samples*sizeof(gr_complex));
        f_up.write((const char*)&up_block[0], d_num_samples*sizeof(gr_complex));
      #endif

      // Preamble and Data FFT
      // If d_fft_size_factor is greater than 1, the rest of the sample buffer will be zeroed out and blend into the window
      memset(d_fft->get_inbuf(),            0, d_fft_size*sizeof(gr_complex));
      memcpy(d_fft->get_inbuf(), &up_block[0], d_num_samples*sizeof(gr_complex));
      d_fft->execute();
      #if DUMP_IQ
        f_fft.write((const char*)d_fft->get_outbuf(), d_fft_size*sizeof(gr_complex));
      #endif

      // Take argmax of returned FFT (similar to MFSK demod)
      max_idx = search_fft_peak(d_fft->get_outbuf(), fft_res_mag, fft_res_add, fft_res_add_c, &max_val);

      d_argmax_history.insert(d_argmax_history.begin(), max_idx);

      if (d_argmax_history.size() > REQUIRED_PREAMBLE_CHIRPS)
      {
        d_argmax_history.pop_back();
      }

      switch (d_state) {
      case S_RESET:
      {
        d_overlaps = OVERLAP_DEFAULT;
        d_offset = 0;
        d_symbols.clear();
        d_argmax_history.clear();
        d_sfd_history.clear();
        d_sync_recovery_counter = 0;
        d_header_received = false;

        d_state = S_PREFILL;

        #if DEBUG >= DEBUG_INFO
          std::cout << "Next state: S_PREFILL" << std::endl;
        #endif

        break;
      }



      case S_PREFILL:
      {
        if (d_argmax_history.size() >= REQUIRED_PREAMBLE_CHIRPS)
        {
          d_state = S_DETECT_PREAMBLE;

          #if DEBUG >= DEBUG_INFO
            std::cout << "Next state: S_DETECT_PREAMBLE" << std::endl;
          #endif
        }
        break;
      }



      // Looks for the same symbol appearing consecutively, signifying the LoRa preamble
      case S_DETECT_PREAMBLE:
      {
        d_preamble_idx = d_argmax_history[0];

        #if DEBUG >= DEBUG_VERBOSE
          std::cout << "PREAMBLE " << d_argmax_history[0] << std::endl;
        #endif

        // Check for discontinuities that exceed some tolerance
        preamble_found = true;
        {
          std::time_t t = std::time(nullptr);
          std::tm tm = *std::gmtime(&t);
          std::ostringstream buffer;
          buffer << std::put_time(&tm, d_filename_template.c_str());
          d_current_filename = buffer.str();
          std::cout << "Writing to " << d_current_filename << std::endl;
        }
        for (int i = 1; i < REQUIRED_PREAMBLE_CHIRPS; i++)
        {
          uint32_t dis = gr::lora::pmod(int(d_preamble_idx) - int(d_argmax_history[i]), d_bin_size);
          if (dis > d_preamble_drift_max && dis < d_bin_size-d_preamble_drift_max)
          {
            preamble_found = false;
          }
        }

        // Advance to SFD/sync discovery if a contiguous preamble is found
        if (preamble_found)
        {
          d_state = S_OUTPUT;

          // move preamble peak to bin zero
          num_consumed = d_num_samples - d_p*d_preamble_idx/d_fft_size_factor;

          #if DEBUG >= DEBUG_INFO
            std::cout << "Next state: S_OUTPUT" << std::endl;
          #endif
        }
        break;
      }

      default:
        break;
      }

      #if DUMP_IQ
        f_raw.write((const char*)&in[0], num_consumed*sizeof(gr_complex));
      #endif

      consume_each (num_consumed);

      volk_free(down_block);
      volk_free(up_block);
      volk_free(fft_res_mag);
      volk_free(fft_res_add);
      volk_free(fft_res_add_c);

      return noutput_items;
    }

  } /* namespace lora */
} /* namespace gr */
