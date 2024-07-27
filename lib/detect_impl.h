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

#ifndef INCLUDED_LORA_DETECT_IMPL_H
#define INCLUDED_LORA_DETECT_IMPL_H

#include <cmath>
#include <cstdlib>
#include <vector>
#include <queue>
#include <complex>
#include <fstream>
#include <gnuradio/fft/fft.h>
#include <gnuradio/fft/window.h>
#include <volk/volk.h>
#include "gnuradio/lora/detect.h"
#include "utilities.h"

namespace gr {
  namespace lora {

    class detect_impl : public detect
    {
     private:
      pmt::pmt_t d_header_port;
      pmt::pmt_t d_out_port;

      detect_state_t d_state;
      uint8_t d_sf;
      uint8_t d_cr;
      uint8_t d_payload_len;
      bool d_crc;
      bool d_ldr;
      bool d_header;
      bool d_header_received;
      bool d_header_valid;

      uint16_t d_num_symbols;
      uint16_t d_fft_size_factor;
      uint32_t d_fft_size;
      uint16_t d_overlaps;
      uint16_t d_offset;
      uint16_t d_p;
      uint32_t d_num_samples;
      uint32_t d_bin_size;
      uint32_t d_preamble_drift_max;

      uint32_t d_packet_symbol_len;

      float d_cfo;

      uint32_t d_preamble_idx;
      uint16_t d_sfd_idx;
      std::vector<uint32_t> d_argmax_history;
      std::vector<uint16_t> d_sfd_history;
      uint16_t d_sync_recovery_counter;

      uint16_t d_peak_search_algorithm;
      uint16_t d_peak_search_phase_k;

      fft::fft_complex_fwd *d_fft;
      std::vector<float>   d_window;
      float                d_beta;

      std::vector<gr_complex> d_upchirp;
      std::vector<gr_complex> d_downchirp;

      std::vector<float> d_symbols;

      std::ofstream f_current_out;
      std::string d_current_filename;

      uint16_t d_on_length;
      std::string d_filename_template;

     public:
      detect_impl( uint8_t   spreading_factor,
                  float     beta,
                  uint16_t  fft_factor,
                  uint8_t   peak_search_algorithm,
                  uint16_t  peak_search_phase_k,
                  float     fs_bw_ratio,
                  uint16_t  on_length,
                  const std::string &filename_template);
      ~detect_impl();

      uint16_t argmax(gr_complex *fft_result);
      uint32_t argmax_32f(float *fft_result, float *max_val_p);
      uint32_t search_fft_peak(const lv_32fc_t *fft_result,
                                   float *buffer1, float *buffer2,
                                   gr_complex *buffer_c, float *max_val_p);
      uint32_t fft_add(const lv_32fc_t *fft_result, float *buffer, gr_complex *buffer_c,
                           float *max_val_p, float phase_offset);

      // Where all the action really happens

      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace lora
} // namespace gr

#endif /* INCLUDED_LORA_DETECT_IMPL_H */
