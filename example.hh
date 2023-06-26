/*
WebUSB and WebAssembly RTL-SDR Example

Copyright 2019 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

extern "C" {

__attribute__((visibility("default")))
void _start();

__attribute__((visibility("default")))
void toggle_rainbow();

__attribute__((visibility("default")))
void toggle_phosphor();

__attribute__((visibility("default")))
void change_irate(int);

__attribute__((visibility("default")))
void change_orate(int);

__attribute__((visibility("default")))
void change_color(int);

__attribute__((visibility("default")))
void process_input();

__attribute__((visibility("default")))
int input_length();

__attribute__((visibility("default")))
unsigned char *input_pointer();

__attribute__((visibility("default")))
void consumed_output(int);

__attribute__((visibility("default")))
int output_length();

__attribute__((visibility("default")))
float *output_pointer();

__attribute__((visibility("default")))
void change_region(int);

__attribute__((visibility("default")))
void change_type(int);

__attribute__((visibility("default")))
int scope_length();

__attribute__((visibility("default")))
int scope_width();

__attribute__((visibility("default")))
int scope_height();

__attribute__((visibility("default")))
int *scope_pointer();

__attribute__((visibility("default")))
int spectrum_length();

__attribute__((visibility("default")))
int spectrum_width();

__attribute__((visibility("default")))
int spectrum_height();

__attribute__((visibility("default")))
int *spectrum_pointer();

__attribute__((visibility("default")))
int spectrogram_length();

__attribute__((visibility("default")))
int spectrogram_width();

__attribute__((visibility("default")))
int spectrogram_height();

__attribute__((visibility("default")))
int *spectrogram_pointer();

}

