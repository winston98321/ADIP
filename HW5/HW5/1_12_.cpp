/*#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <queue>
#include <complex>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace cv;
#define PI 3.141592653589793
using namespace std;
clock_t a, b;
vector<vector<int>> read_raw(const char filename[], int sizew, int sizeh) {
	FILE* input_file;
	unsigned char* image = new unsigned char[sizew * sizeh];
	input_file = fopen(filename, "rb");
	fread(image, 1, sizew * sizeh, input_file);
	vector<vector<int>> img;
	for (int i = 0; i < sizeh; i++) {
		vector<int> vec;
		for (int j = 0; j < sizew; j++) {
			vec.push_back(int(image[i * sizew + j]));
		}img.push_back(vec);
	}
	return img;
}
void write_raw(vector<vector<int>> img, const char filename[], int sizew, int sizeh) {
	FILE* output_file;
	unsigned char* image = new unsigned char[sizew * sizeh];
	output_file = fopen(filename, "wb");
	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			image[i * sizew + j] = img[i][j];
		}
	}
	fwrite(image, 1, sizeh * sizew, output_file);
	fclose(output_file);
}

vector<vector<complex<double>>> dft(vector<vector<int>> img, const char filename[], int sizew, int sizeh) {
	vector<vector<complex<double>>> f1_out(sizeh, vector<complex<double>>(sizew));
	vector<vector<complex<double>>> f2_out(sizeh, vector<complex<double>>(sizew));
	vector<vector<int>> img_out(sizeh, vector<int>(sizew, 0));
	auto M_I_2PI_DL = -(6.28318530718i / double(sizew));
	a = clock();
	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			complex<double> sum = (0,0);
			for (int k = 0; k < sizew; k++) {
				sum += double(img[i][k]) *exp(M_I_2PI_DL * double(k) * double(j)) * exp(M_I_2PI_DL * double(k) * double(128));
			}
			f1_out[i][j] = sum;
		}
	}


	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			complex<double> sum = (0,0);
			for (int k = 0; k < sizew; k++) {
				sum += f1_out[k][i] * std::exp(M_I_2PI_DL * double(k) * double(j)) * exp(M_I_2PI_DL * double(k) * double(128));//
			}
			f2_out[j][i] = sum;
		}
	}

	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			//cout << f2_out[i][j].real() << " " << f2_out[i][j].imag() << endl;
			double ans = 255*(round(sqrt(f2_out[i][j].real() * f2_out[i][j].real() + f2_out[i][j].imag() * f2_out[i][j].imag()))/(256.0*256));
			//cout << ans << endl;
			//if (ans > 255)ans = 255;
			//if (ans < 0) ans = 0;
			img_out[i][j] = ans;
		}
	}
	img = img_out;
	//log transfer for image enhancement
	double c = 255 / (log(1 + 255));
	for (int i = 0; i < sizeh; ++i) {
		for (int j = 0; j < sizew; ++j) {
			img_out[i][j] = c * log(1 + img_out[i][j]);
		}

	}
	double ma = -999;
	for (int i = 0; i < sizeh; ++i) {
		for (int j = 0; j < sizew; ++j) {
			if (img_out[i][j] > ma) {
				ma = img_out[i][j];
			}
		}

	}
	for (int i = 0; i < sizeh; ++i) {
		for (int j = 0; j < sizew; ++j) {
				img_out[i][j] = 255*(img_out[i][j]/ma);
				//cout << ma<<endl;
		}

	}
	b = clock();
	cout << "DFT Caluclate time :"<<double(b - a) / CLOCKS_PER_SEC << endl;
	write_raw(img_out, filename, sizew, sizeh);
	return f2_out;
}
void idft(vector<vector<complex<double>>> img, const char filename[], int sizew, int sizeh)
{
	std::vector<std::vector<complex<double>>> f1_out(sizeh, std::vector<complex<double>>(sizew));
	std::vector<std::vector<complex<double>>> f2_out(sizeh, std::vector<complex<double>>(sizew));
	std::vector<std::vector<int>> img_out(sizeh, std::vector<int>(sizew, 0));
	auto M_I_2PI_DL = (6.28318530718i / double(sizew));
	a = clock();
	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			complex<double> sum = (0, 0);
			for (int k = 0; k < sizew; k++) {
				sum += img[i][k] * exp(M_I_2PI_DL * double(k) * double(j));// * exp(M_I_2PI_DL * double(k) * double(128));
			}
			f1_out[i][j] = sum;
		}
	}



	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			complex<double> sum = (0, 0);
			for (int k = 0; k < sizew; k++) {
				sum += f1_out[k][i] * exp(M_I_2PI_DL * double(k) * double(j));// *exp(M_I_2PI_DL * double(k) * double(128));
			}
			f2_out[j][i] = sum;
		}
	}

	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			//cout << f2_out[i][j].real() << " " << f2_out[i][j].imag() << endl;
			int ans = round(sqrt(f2_out[i][j].real() * f2_out[i][j].real() + f2_out[i][j].imag() * f2_out[i][j].imag()))/(256*256) ;
			img_out[i][j] = ans;
		}
	}
	b = clock();
	cout << "IDFT Caluclate time :" << double(b - a) / CLOCKS_PER_SEC << endl;
	write_raw(img_out, filename, sizew, sizeh);
}
void compare(vector<vector<int>> vec, vector<vector<int>> src, string a) {
	double MSE = 0;
	int len = vec[0].size();
	for (int i = 0; i < len; i++) {
		for (int j = 0; j < len; j++) {
			MSE += pow((src[i][j] - vec[i][j]), 2);
		}
	}
	MSE = (MSE / (len * len));
	cout << a << "_MSE :" << MSE << endl;
	cout << a << "_PSNR :" << 10 * log10((255 * 255) / MSE) << endl;
}
int main() 
{
	vector<vector<int>>baboon, lena,rec_baboon,rec_lena;
	vector<vector<complex<double>>> idft_baboon(256, vector<complex<double>>(256));
	vector<vector<complex<double>>> idft_lena(256, vector<complex<double>>(256));
	baboon = read_raw("src/baboon_256.raw", 256, 256);
	lena = read_raw("src/lena_256.raw", 256, 256);
	
	idft_baboon = dft(baboon, "image_file/dft_baboon.raw", 256, 256);
	idft_lena = dft(lena, "image_file/dft_lena.raw", 256, 256);
	idft(idft_baboon, "image_file/idft_baboon.raw", 256, 256);
	idft(idft_lena, "image_file/idft_lena.raw", 256, 256);
	rec_baboon = read_raw("image_file/idft_baboon.raw", 256, 256);
	rec_lena = read_raw("image_file/idft_lena.raw", 256, 256);
	string a = "1.2 baboon: ";
	compare(baboon, rec_baboon, a);
	a = "1.2 lena: ";
	compare(lena, rec_lena, a);
	return 0;
}
*/