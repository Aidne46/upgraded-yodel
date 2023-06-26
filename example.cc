/*
WebUSB and WebAssembly RTL-SDR Example

Copyright 2019 Ahmet Inan <inan@aicodix.de>
*/

#include "quirks.hh"
#include "complex.hh"
#include "example.hh"
#include "decibel.hh"
#include "fft.hh"
#include "ema.hh"
#include "fmd.hh"
#include "biquad.hh"
#include "normalize.hh"
#include "calculus.hh"

typedef DSP::Complex<float> cmplx;

static const int SAMPLES = 32768;
static const int BINS = 1024;
static const int INPUT = 32768;
unsigned char buffer[2 * INPUT];
cmplx samples[SAMPLES];
cmplx input[INPUT];
cmplx output[BINS];
DSP::FastFourierTransform<BINS, cmplx, -1> fwd;

bool show_rainbow = true;
bool show_phosphor = true;
int draw_color = 0;
int cur_irate = 0;
int cur_orate = 0;
int cur_ipos = 0;
int cur_opos = 0;
int cur_knee = 15000;
int cur_cutoff = 75000;
bool mono = false;

template <typename TYPE, int WIDTH, int HEIGHT>
struct Image {
	static const int width = WIDTH;
	static const int height = HEIGHT;
	static const int length = WIDTH * HEIGHT;
	TYPE pixels[length];
	static TYPE overwrite(int, int, TYPE, TYPE c)
	{
		return c;
	}
	template <typename OP>
	void set(int x, int y, TYPE c, OP op)
	{
		if (0 <= x && x < width && 0 <= y && y < height)
			pixels[width * y + x] = op(x, y, pixels[width * y + x], c);
	}
	void set(int x, int y, TYPE c)
	{
		if (0 <= x && x < width && 0 <= y && y < height)
			pixels[width * y + x] = c;
	}
	template <typename OP>
	void fill(TYPE c, OP op)
	{
		for (int y = 0; y < height; ++y)
			for (int x = 0; x < width; ++x)
				pixels[width * y + x] = op(x, y, pixels[width * y + x], c);
	}
	void fill(TYPE c)
	{
		for (int i = 0; i < length; ++i)
			pixels[i] = c;
	}
	template <typename OP>
	void vline(int x, TYPE c, OP op)
	{
		if (0 <= x && x < width)
			for (int i = 0; i < height; ++i)
				set(x, i, c, op);
	}
	void vline(int x, TYPE c)
	{
		if (0 <= x && x < width)
			for (int i = 0; i < height; ++i)
				set(x, i, c);
	}
	template <typename OP>
	void hline(int y, TYPE c, OP op)
	{
		if (0 <= y && y < height)
			for (int i = 0; i < width; ++i)
				set(i, y, c, op);
	}
	void hline(int y, TYPE c)
	{
		if (0 <= y && y < height)
			for (int i = 0; i < width; ++i)
				set(i, y, c);
	}
	// only made for abs(x0-x1) <= 1
	template <typename OP>
	void vline(int x0, int y0, int x1, int y1, TYPE c, OP op)
	{
		int a0 = min(y0, (y0 + y1) / 2);
		int a1 = max(y0, (y0 + y1) / 2);
		for (int y = a0; y <= a1; ++y)
			set(x0, y, c, op);
		int b0 = min((y0 + y1) / 2, y1);
		int b1 = max((y0 + y1) / 2, y1);
		for (int y = b0; y <= b1; ++y)
			set(x1, y, c, op);
	}
	void vline(int x0, int y0, int x1, int y1, TYPE c)
	{
		vline(x0, y0, x1, y1, c, overwrite);
	}
	// only made for abs(y0-y1) <= 1
	template <typename OP>
	void hline(int x0, int y0, int x1, int y1, TYPE c, OP op)
	{
		int a0 = min(x0, (x0 + x1) / 2);
		int a1 = max(x0, (x0 + x1) / 2);
		for (int x = a0; x <= a1; ++x)
			set(x, y0, c, op);
		int b0 = min((x0 + x1) / 2, x1);
		int b1 = max((x0 + x1) / 2, x1);
		for (int x = b0; x <= b1; ++x)
			set(x, y1, c, op);
	}
	void hline(int x0, int y0, int x1, int y1, TYPE c)
	{
		hline(x0, y0, x1, y1, c, overwrite);
	}
	template <typename OP>
	void low_angle(int x0, int y0, int x1, int y1, TYPE c, OP op)
	{
		int dx = x1 - x0;
		int dy = y1 - y0;
		int yi = 1;
		if (dy < 0) {
			yi = -1;
			dy = -dy;
		}
		int D = 2 * dy - dx;
		int y = y0;
		for (int x = x0; x <= x1; ++x) {
			set(x, y, c, op);
			if (D > 0) {
				y += yi;
				D -= 2 * dx;
			}
			D += 2 * dy;
		}
	}
	template <typename OP>
	void high_angle(int x0, int y0, int x1, int y1, TYPE c, OP op)
	{
		int dx = x1 - x0;
		int dy = y1 - y0;
		int xi = 1;
		if (dx < 0) {
			xi = -1;
			dx = -dx;
		}
		int D = 2 * dx - dy;
		int x = x0;
		for (int y = y0; y <= y1; ++y) {
			set(x, y, c, op);
			if (D > 0) {
				x += xi;
				D -= 2 * dy;
			}
			D += 2 * dx;
		}
	}
	template <typename OP>
	void line(int x0, int y0, int x1, int y1, TYPE c, OP op)
	{
		if (abs(y0 - y1) < abs(x0 - x1)) {
			if (x0 < x1)
				low_angle(x0, y0, x1, y1, c, op);
			else
				low_angle(x1, y1, x0, y0, c, op);
		} else {
			if (y0 < y1)
				high_angle(x0, y0, x1, y1, c, op);
			else
				high_angle(x1, y1, x0, y0, c, op);
		}
	}
	void line(int x0, int y0, int x1, int y1, TYPE c)
	{
		line(x0, y0, x1, y1, c, overwrite);
	}
};

static const int SCOPE_WIDTH = 256;
static const int SCOPE_HEIGHT = 256;
Image<int, SCOPE_WIDTH, SCOPE_HEIGHT> scope;

static const int SPECTRUM_WIDTH = BINS;
static const int SPECTRUM_HEIGHT = SPECTRUM_WIDTH / 4;
Image<int, SPECTRUM_WIDTH, SPECTRUM_HEIGHT> spectrum;

static const int SPECTROGRAM_WIDTH = SPECTRUM_WIDTH;
static const int SPECTROGRAM_HEIGHT = SPECTRUM_HEIGHT;
Image<int, SPECTROGRAM_WIDTH, SPECTROGRAM_HEIGHT> spectrogram;

int rgba(float r, float g, float b, float a)
{
	r = min(max(r, 0.f), 1.f);
	g = min(max(g, 0.f), 1.f);
	b = min(max(b, 0.f), 1.f);
	a = min(max(a, 0.f), 1.f);
	r = sqrt(r);
	g = sqrt(g);
	b = sqrt(b);
	//a = sqrt(a);
	int R = nearbyint(255.f * r);
	int G = nearbyint(255.f * g);
	int B = nearbyint(255.f * b);
	int A = nearbyint(255.f * a);
	return (A << 24) | (B << 16) | (G << 8) | (R << 0);
}

int alpha(float a)
{
	a = min(max(a, 0.f), 1.f);
	int A = nearbyint(255.f * a);
	return A << 24;
}

int rainbow(float v)
{
	v = min(max(v, 0.f), 1.f);
	float r = 4.f * v - 2.f;
	float g = 1.f - 4.f * abs(v - .5f);
	float b = 2.f - 4.f * v;
	float a = 4.f * v;
	return rgba(r, g, b, a);
}

void _start()
{
	for (int i = 0; i < BINS; ++i)
		spectrogram.vline(i, rainbow((float)i/BINS));
}

cmplx convert(unsigned char *inp)
{
	float real(inp[0]), imag(inp[1]),
		bias(127.5), scale(0.0078125);
	return scale * cmplx(real - bias, imag - bias);
}

DSP::EMACascade<cmplx, float, 3> iflp;
DSP::NormalizeIQ<cmplx> normalize;
DSP::FMD5<cmplx> demod;
DSP::Biquad<float, float> notch, bandpass;
DSP::Differentiator<float> diff;
DSP::BiquadCascade<cmplx, float, 2> aalp;
DSP::EMA<cmplx, float> deemp, bias;

void produce_output()
{
	for (int i = 0; i < INPUT; ++i) {
		cmplx iq = input[i];
		iq = iflp(iq);
		iq = normalize(iq);
		float phase = demod(iq);
		float pilot = bandpass(phase);
		float level = notch(phase);
		float twice = pilot * diff(pilot);
		float l = 0.f, r = 0.f;
		if (mono)
			l = r = level;
		else if (twice >= 0.f)
			l = 2.f * level;
		else
			r = 2.f * level;
		cmplx tmp(l, r);
		tmp = aalp(tmp);
		if (cur_ipos >= cur_irate) {
			tmp = deemp(tmp);
			tmp -= bias(tmp);
			samples[cur_opos] = tmp;
			cur_opos = min(cur_opos+1, SAMPLES-1);
			cur_ipos -= cur_irate;
		}
		cur_ipos += cur_orate;
	}
}

float dB_avg[BINS], avg_pwr[BINS], avg_max;

void process_input()
{
	for (int i = 0; i < INPUT; ++i)
		input[i] = convert(buffer + 2 * i);
	produce_output();
	if (show_phosphor) {
		for (int i = 0; i < scope.length; ++i)
			scope.pixels[i] = draw_color |
				(((((unsigned)scope.pixels[i]>>24)*7)>>3)<<24);
		for (int i = 0; i < spectrum.length; ++i)
			spectrum.pixels[i] = draw_color |
				(((((unsigned)spectrum.pixels[i]>>24)*7)>>3)<<24);
	} else {
		scope.fill(draw_color);
		spectrum.fill(draw_color);
	}
	for (int i = 0; i < INPUT; ++i) {
		int real(buffer[2*i+0]), imag(buffer[2*i+1]);
		scope.pixels[scope.width*imag+real] = 0xff000000 | draw_color;
	}
	float tmp_max = -1.f;
	for (int k = 0; k < INPUT; k += BINS) {
		fwd(output, input + k);
		for (int i = 0; i < BINS; ++i) {
			float scale = 1.f / BINS;
			float pwr = norm(scale * output[i]);
			tmp_max = max(tmp_max, pwr);
			avg_pwr[i] = DSP::lerp(avg_pwr[i], pwr, 0.1f);
		}
	}
	for (int i = 0; i < (BINS+1)/2; ++i)
		dB_avg[i+BINS/2] = DSP::decibel(avg_pwr[i]);
	for (int i = (BINS+1)/2; i < BINS; ++i)
		dB_avg[i-(BINS+1)/2] = DSP::decibel(avg_pwr[i]);
	float tmp_min = dB_avg[0];
	for (int i = 1; i < BINS; ++i)
		tmp_min = min(tmp_min, dB_avg[i]);
	float dB_min = max(tmp_min, -120.f);
	avg_max = DSP::lerp(avg_max, tmp_max, avg_max < tmp_max ? 0.1f : 0.001f);
	float dB_max = DSP::decibel(avg_max);
	spectrum.vline(BINS/2, 0xff7f7f7f);
	for (int b = 0, i0, j0; b < BINS; ++b) {
		float dB = dB_avg[b];
		float scale = (spectrum.height-1) / (dB_min-dB_max);
		int v = nearbyint(scale * (dB-dB_max));
		int i1 = b;
		int j1 = min(max(v, 0), spectrum.height-1);
		if (b)
			spectrum.vline(i0, j0, i1, j1, 0xff000000 | draw_color);
		i0 = i1;
		j0 = j1;
	}
	for (int j = spectrogram.height-1; j; --j)
		for (int i = 0; i < spectrogram.width; ++i)
			spectrogram.pixels[spectrogram.width*j+i] = spectrogram.pixels[spectrogram.width*(j-1)+i];
	for (int b = 0; b < BINS; ++b) {
		float dB = dB_avg[b];
		float scale = (1.f / (dB_min-dB_max));
		float v = scale * (dB-dB_max);
		if (show_rainbow)
			spectrogram.pixels[b] = rainbow(1.f-v);
		else
			spectrogram.pixels[b] = alpha(1.f-v) | draw_color;
	}
}

void consumed_output(int len)
{
	cur_opos = max(cur_opos - len, 0);
	for (int i = 0; i < SAMPLES-len; ++i)
		samples[i] = samples[i+len];
}

void change_irate(int rate)
{
	cur_irate = rate;
	iflp.cutoff(cur_cutoff, rate);
	normalize.samples(rate / 10);
	demod.bandwidth(150000.f / rate);
	notch.notch(19000, rate, 24);
	bandpass.bandpass(19000, rate, 24);
	aalp.lowpass(15000, rate);
}

void change_orate(int rate)
{
	cur_orate = rate;
	deemp.cutoff(cur_knee, rate);
	bias.cutoff(30, rate);
}

void change_region(int region)
{
	switch (region) {
	case 1:
		// Europe and Asia: 50µs or 3183Hz
		cur_knee = 3183;
		break;
	case 2:
		// The Americas: 75µs or 2122Hz
		cur_knee = 2122;
		break;
	default:
		cur_knee = 15000;
	}
	deemp.cutoff(cur_knee, cur_orate);
}

void change_type(int type)
{
	switch (type) {
		case 1:
			// Mono WFM
			cur_cutoff = 75000;
			mono = true;
			break;
		case 2:
			// Mono NFM
			cur_cutoff = 25000;
			mono = true;
			break;
		default:
			// Stereo WFM
			cur_cutoff = 75000;
			mono = false;
	}
	iflp.cutoff(cur_cutoff, cur_irate);
}

void toggle_rainbow()
{
	show_rainbow = !show_rainbow;
}

void toggle_phosphor()
{
	show_phosphor = !show_phosphor;
}

void change_color(int color)
{
	draw_color = 0x00ffffff & color;
}

int input_length()
{
	return INPUT;
}

unsigned char *input_pointer()
{
	return buffer;
}

int output_length()
{
	return SAMPLES;
}

float *output_pointer()
{
	return reinterpret_cast<float *>(samples);
}

int scope_length()
{
	return scope.length;
}

int scope_width()
{
	return scope.width;
}

int scope_height()
{
	return scope.height;
}

int *scope_pointer()
{
	return scope.pixels;
}

int spectrum_length()
{
	return spectrum.length;
}

int spectrum_width()
{
	return spectrum.width;
}

int spectrum_height()
{
	return spectrum.height;
}

int *spectrum_pointer()
{
	return spectrum.pixels;
}

int spectrogram_length()
{
	return spectrogram.length;
}

int spectrogram_width()
{
	return spectrogram.width;
}

int spectrogram_height()
{
	return spectrogram.height;
}

int *spectrogram_pointer()
{
	return spectrogram.pixels;
}

