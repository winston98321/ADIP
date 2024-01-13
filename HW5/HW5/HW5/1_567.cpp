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
vector<vector<double>> read_raw(const char filename[], int sizew, int sizeh) {
	FILE* input_file;
	unsigned char* image = new unsigned char[sizew * sizeh];
	input_file = fopen(filename, "rb");
	fread(image, 1, sizew * sizeh, input_file);
	vector<vector<double>> img;
	for (int i = 0; i < sizeh; i++) {
		vector<double> vec;
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
void dct1D(vector<double> img,  vector<double> &output, int size) {
	for (int k = 0; k < size; k++) {
		double sum = 0.0;
		for (int n = 0; n < size; n++) {
			sum += img[n] * cos((2.0 * n + 1.0) * k * PI / (2.0 * size));
		}

		output.push_back( sum * ((k == 0) ? 1.0 / sqrt(2.0) : 1.0));
	}
}
void idct1D(vector<double> img, vector<double>& output, int size) {
	for (int k = 0; k < size; k++) {
		double sum = 0.0;
		for (int n = 0; n < size; n++) {
			sum += ((n == 0) ? 1.0 / sqrt(2.0) : 1.0)*img[n] * cos((2.0 * k + 1.0) * n * PI / (2.0 * size));

		}

		output.push_back(sum);
	}
}
void dct2D(vector<vector<double>> input, vector<vector<double>>&output, int sizew, int sizeh, const char filename[]) {
	a = clock();
	for (int i = 0; i < sizew; i++) {
		for (int j = 0; j < sizeh; j++) {
			vector<double> rowInput, rowOutput;
			for (int k = 0; k < sizeh; k++) {
				rowInput.push_back(double(input[i][k]));
			}

			dct1D(rowInput, rowOutput, sizeh);

			for (int k = 0; k < sizeh; k++) {
				output[i][k] = rowOutput[k];
			}
		}
	}

	for (int j = 0; j < sizew; j++) {
		vector<double> colInput, colOutput;
		for (int i = 0; i < sizew; i++) {
			colInput.push_back(output[i][j]);
		}
		dct1D(colInput, colOutput, sizew);
		for (int i = 0; i < sizew; i++) {
			output[i][j] = colOutput[i];
		}
	}
	b = clock();
	cout << "DCT Caluclate time :" << double(b - a) / CLOCKS_PER_SEC << endl;
	vector<vector<int>>img_out(256, vector<int>(256, 0));
	double ma = -9999999, mi = 99999999;
	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			output[i][j] = (2 / (sqrt(sizew * sizew))) * (output[i][j]);
			img_out[i][j] = round(abs(output[i][j]));
		}
	}
	

	write_raw(img_out, filename, sizew, sizeh);
	
}
void idct2D(vector<vector<double>> input, vector<vector<double>>& output, int sizew, int sizeh, const char filename[]) {
	a = clock();
	for (int i = 0; i < sizew; i++) {
		for (int j = 0; j < sizeh; j++) {
			vector<double> rowInput, rowOutput;
			for (int k = 0; k < sizeh; k++) {
				rowInput.push_back(double(input[i][k]));
			}

			idct1D(rowInput, rowOutput, sizeh);

			for (int k = 0; k < sizeh; k++) {
				output[i][k] = rowOutput[k];
			}
		}
	}

	for (int j = 0; j < sizew; j++) {
		vector<double> colInput, colOutput;
		for (int i = 0; i < sizew; i++) {
			colInput.push_back(output[i][j]);
		}
		idct1D(colInput, colOutput, sizew);
		for (int i = 0; i < sizew; i++) {
			output[i][j] = colOutput[i];
		}
	}
	b = clock();
	cout << "IDCT Caluclate time :" << double(b - a) / CLOCKS_PER_SEC << endl;
	vector<vector<int>>img_out(256, vector<int>(256, 0));
	double ma = -9999999, mi = 99999999;
	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			output[i][j] = (2 / (sqrt(sizew * sizew))) * (output[i][j]);
			img_out[i][j] = round(abs(output[i][j]));

		}
	}
	//cout << ma << " " << mi << endl;

	write_raw(img_out, filename, sizew, sizeh);

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
void compare(vector<vector<double>> vec, vector<vector<double>> src, string a) {
	double MSE = 0;
	int len = vec[0].size();
	for (int i = 0; i < len; i++) {
		for (int j = 0; j < len; j++) {
			MSE += pow((round(src[i][j]) - round(vec[i][j])), 2);
		}
	}
	MSE = (MSE / (len * len));
	cout << a << "_MSE :" << MSE << endl;
	cout << a << "_PSNR :" << 10 * log10((255 * 255) / MSE) << endl;
}
int main()
{
	vector<vector<double>>baboon, lena,rec_baboon,rec_lena;
	vector<vector<double>> idct_baboon(256, vector<double>(256,0.0)), dct_baboon(256, vector<double>(256, 0.0));
	vector<vector<double>> idct_lena(256, vector<double>(256.0,0.0)), dct_lena(256, vector<double>(256, 0.0));
	baboon = read_raw("src/baboon_256.raw", 256, 256);
	lena = read_raw("src/lena_256.raw", 256, 256);
	int ma = -9999, mi = 9999;
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			if (baboon[i][j] > ma)ma = baboon[i][j];
			if (baboon[i][j] < mi)mi = baboon[i][j];
		}
	}
	cout << ma<<" " << mi << endl;
	dct2D(baboon,dct_baboon,256,256, "image_file/dct_baboon.raw");
	idct2D(dct_baboon,idct_baboon, 256, 256, "image_file/idct_baboon.raw");
	
	dct2D(lena, dct_lena, 256, 256, "image_file/dct_lena.raw");
	idct2D(dct_lena, idct_lena, 256, 256, "image_file/idct_lena.raw");
	string s = "1.6 lena : ";
	compare(lena, idct_lena, s);
	s = "1.6 baboon : ";
	compare(baboon, idct_baboon, s);
	return 0;
}
*/