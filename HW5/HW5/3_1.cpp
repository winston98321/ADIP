#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <queue>
#include <complex>
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
vector<vector<complex<double>>> fdft(vector<vector<int>> img, const char filename[], int sizew, int sizeh) {
	vector<vector<complex<double>>> f1_out(sizeh, vector<complex<double>>(sizew));
	vector<vector<complex<double>>> f2_out(sizeh, vector<complex<double>>(sizew));
	vector<vector<double>> img_out(sizeh, vector<double>(sizew, 0));
	vector<vector<int>> img_out2(sizeh, vector<int>(sizew, 0));
	auto M_I_2PI_DL = -(6.28318530718i / double(sizew));
	a = clock();
	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			complex<double> sum = (0, 0);
			for (int k = 0; k < sizew; k++) {
				sum += double(img[i][k]) * exp(M_I_2PI_DL * double(k) * double(j)) * exp(M_I_2PI_DL * double(k) * double(sizew / 2));

			}
			f1_out[i][j] = sum;
		}
	}

	M_I_2PI_DL = -(6.28318530718i / double(sizeh));

	for (int i = 0; i < sizew; i++) {
		for (int j = 0; j < sizeh; j++) {
			complex<double> sum = (0, 0);
			for (int k = 0; k < sizeh; k++) {
				sum += f1_out[k][i] * std::exp(M_I_2PI_DL * double(k) * double(j)) * exp(M_I_2PI_DL * double(k) * double(sizeh / 2));//
			}
			f2_out[j][i] = sum;
		}
	}

	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			//cout << f2_out[i][j].real() << " " << f2_out[i][j].imag() << endl;
			double ans =  255*(sqrt(f2_out[i][j].real() * f2_out[i][j].real() + f2_out[i][j].imag() * f2_out[i][j].imag())) / (sizew * sizeh);

			//cout << ans << endl;
			//if (ans > 255)ans = 255;
			//if (ans < 0) ans = 0;
			img_out[i][j] = ans*255;
			cout << img_out[i][j] << endl;
		}
	}
	//img = img_out;
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
			img_out[i][j] =int(round( 255.0 * (img_out[i][j] / ma)));
			img_out2[i][j] = img_out[i][j];
			//cout << ma<<endl;
		}

	}
	
	b = clock();
	cout << "DFT Caluclate time :" << double(b - a) / CLOCKS_PER_SEC << endl;
	write_raw(img_out2, filename, sizew, sizeh);
	return f2_out;
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
				sum += double(img[i][k]) *exp(M_I_2PI_DL * double(k) * double(j)) * exp(M_I_2PI_DL * double(k) * double(sizew/2));
			}
			f1_out[i][j] = sum;
		}
	}

	M_I_2PI_DL = -(6.28318530718i / double(sizeh));

	for (int i = 0; i < sizew; i++) {
		for (int j = 0; j < sizeh; j++) {
			complex<double> sum = (0,0);
			for (int k = 0; k < sizeh; k++) {
				sum += f1_out[k][i] * std::exp(M_I_2PI_DL * double(k) * double(j)) * exp(M_I_2PI_DL * double(k) * double(sizeh/2));//
			}
			f2_out[j][i] = sum;
		}
	}

	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			//cout << f2_out[i][j].real() << " " << f2_out[i][j].imag() << endl;
			double ans = 255*(round(sqrt(f2_out[i][j].real() * f2_out[i][j].real() + f2_out[i][j].imag() * f2_out[i][j].imag()))/(sizew*sizeh));
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
	/* 
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

	}*/
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

	M_I_2PI_DL = (6.28318530718i / double(sizeh));

	for (int i = 0; i < sizew; i++) {
		for (int j = 0; j < sizeh; j++) {
			complex<double> sum = (0, 0);
			for (int k = 0; k < sizeh; k++) {
				sum += f1_out[k][i] * exp(M_I_2PI_DL * double(k) * double(j));// *exp(M_I_2PI_DL * double(k) * double(128));
			}
			f2_out[j][i] = sum;
		}
	}

	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			//cout << f2_out[i][j].real() << " " << f2_out[i][j].imag() << endl;
			int ans = round(sqrt(f2_out[i][j].real() * f2_out[i][j].real() + f2_out[i][j].imag() * f2_out[i][j].imag()))/(sizeh*sizew) ;
			img_out[i][j] = ans;
			//cout << ans << endl;
		}
	}
	b = clock();
	cout << "IDFT Caluclate time :" << double(b - a) / CLOCKS_PER_SEC << endl;
	write_raw(img_out, filename, sizew, sizeh);
}
void iidft(vector<vector<complex<double>>> img, const char filename[], int sizew, int sizeh)
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
				sum += img[i][k] * exp(M_I_2PI_DL * double(k) * double(j)) *exp(M_I_2PI_DL * double(k) * double(sizew / 2));
				
			}
			f1_out[i][j] = sum;
		}
	}

	M_I_2PI_DL = (6.28318530718i / double(sizeh));

	for (int i = 0; i < sizew; i++) {
		for (int j = 0; j < sizeh; j++) {
			complex<double> sum = (0, 0);
			for (int k = 0; k < sizeh; k++) {
				sum += f1_out[k][i] * exp(M_I_2PI_DL * double(k) * double(j)) *exp(M_I_2PI_DL * double(k) * double(sizeh / 2));
			}
			f2_out[j][i] = sum;
		}
	}

	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			//cout << f2_out[i][j].real() << " " << f2_out[i][j].imag() << endl;
			int ans = round(sqrt(f2_out[i][j].real() * f2_out[i][j].real() + f2_out[i][j].imag() * f2_out[i][j].imag()) / (sizeh * sizew));
			
			if (ans > 255)ans = 255;
			if (ans <0)ans = 0;
			img_out[i][j] = ans;
			cout << ans << endl;
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
void sobel(vector<vector<int>> img, int kernel_angle, const char filename[], int sizew, int sizeh) {
	int bund = 1;
	vector<vector<int>> out = img;
	vector<int> prewitt_kernel_0 = { 0,1,0,-1,2,0,-2,1,0,-1 };
	vector<int> prewitt_kernel_90 = { 0,-1,-2,-1,0,0,0,1,2,1 };
	vector<int> kernel = prewitt_kernel_0;
	if (kernel_angle == 90) {
		kernel = prewitt_kernel_90;
	}
	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			double total = 0;
			double combine0 = 0;
			double combine90 = 0;
			int kernel_cnt = 0;
			if (i + bund > sizeh - 1 || i - bund<0 || j + bund>sizew - 1 || j - bund < 0) {
				for (int ii = i - bund; ii <= i + bund; ii++) {
					for (int jj = j - bund; jj <= j + bund; jj++) {
						kernel_cnt += 1;
						int iii = ii, jjj = jj;
						if (ii < 0) iii = abs(ii) - 1;
						if (jj < 0) jjj = abs(jj) - 1;
						if (ii > sizeh - 1) iii = sizeh - 1 + (sizeh - 1 - ii + 1);
						if (jj > sizew - 1) jjj = sizew - 1 + (sizew - 1 - jj + 1);
						total += double(kernel[kernel_cnt]) * img[iii][jjj];
						combine0 += double(prewitt_kernel_0[kernel_cnt]) * img[iii][jjj];
						combine90 += double(prewitt_kernel_90[kernel_cnt]) * img[iii][jjj];
					}
				}
			}
			else {
				for (int ii = i - bund; ii <= i + bund; ii++) {
					for (int jj = j - bund; jj <= j + bund; jj++) {
						kernel_cnt += 1;
						total += double(kernel[kernel_cnt]) * img[ii][jj];
						combine0 += double(prewitt_kernel_0[kernel_cnt]) * img[ii][jj];
						combine90 += double(prewitt_kernel_90[kernel_cnt]) * img[ii][jj];
					}
				}
			}
			//if (abs(total > 100))total = 255;
			//else total = 0;
			int cb = pow((pow(combine0, 2) + pow(combine90, 2)), 0.5);

			out[i][j] = total;
			if (kernel_angle == 1) out[i][j] = cb;
			if (out[i][j] > 255)out[i][j] = 255;
			if (out[i][j] < 0)out[i][j] = 0;
		}
	}

	write_raw(out, filename, sizew, sizeh);
}
int main()
{
	vector<vector<int>>building,rec_building;
	vector<vector<complex<double>>> idft_building(474, vector<complex<double>>(632));
	vector<vector<complex<double>>> idft_ans(474, vector<complex<double>>(632));
	vector<vector<complex<double>>> idft_building2(474, vector<complex<double>>(632));
	vector<vector<complex<double>>> idft_sobel0(474, vector<complex<double>>(632));
	vector<vector<complex<double>>> idft_sobel90(474, vector<complex<double>>(632));

	vector<vector<int>> sobel_0(632,vector<int>(474,0));
	vector<vector<int>> sobel_90(632, vector<int>(474, 0));
	
	sobel_0[632 / 2 - 1][474 / 2 + 1] = -1;
	sobel_0[632 / 2 ][474 / 2 + 1] = -2;
	sobel_0[632 / 2 + 1][474 / 2 + 1] = -1;
	sobel_0[632 / 2 - 1][474 / 2 -1 ] = 1;
	sobel_0[632 / 2 ][474 / 2 - 1] = 2;
	sobel_0[632 / 2 + 1][474 / 2 - 1] = 1;
	building = read_raw("src/building_474x632.raw", 474, 632);
	idft_sobel0 = fdft(sobel_0, "image_file/dft_sobel0.raw", 474, 632);
	idft_building = dft(building, "image_file/dft_building.raw", 474, 632);
	idft_ans = idft_building;
	for (int i = 0; i < 632; i++) {
		for (int j = 0; j < 474; j++) {
			idft_ans[i][j] = idft_building[i][j] * idft_sobel0[i][j] ;
		}
	}

	iidft(idft_ans, "image_file/idft_sobel_0_building.raw", 474, 632);
	sobel(building, 0, "image_file/spetial_sobel_0_suilding.raw", 474, 632);
	//sobel90
	

	sobel_90[632 / 2 - 1][474 / 2 + 1] = -1;
	sobel_90[632 / 2][474 / 2 - 1] = -2;
	sobel_90[632 / 2 - 1][474 / 2 - 1] = -1;
	sobel_90[632 / 2 + 1][474 / 2 + 1] = 1;
	sobel_90[632 / 2 +1][474 / 2 ] = 2;
	sobel_90[632 / 2 + 1][474 / 2 - 1] = 1;
	building = read_raw("src/building_474x632.raw", 474, 632);
	idft_sobel90 = fdft(sobel_90, "image_file/dft_sobel90.raw", 474, 632);
	idft_building2 = idft_building;
	
	for (int i = 0; i < 632; i++) {
		for (int j = 0; j < 474; j++) {
			idft_building2[i][j] = idft_building[i][j] * idft_sobel90[i][j];
		}
	}

	iidft(idft_building2, "image_file/idft_sobel_90_building.raw", 474, 632);
	sobel(building, 1, "image_file/spetial_sobel_90_suilding.raw", 474, 632);
	return 0;
}
