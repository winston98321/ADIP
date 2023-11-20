#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <queue>
#include <bitset>

using namespace std;


vector<vector<int>> read_raw(const char filename[],int size) {
	FILE* input_file;
	unsigned char* image = new unsigned char[size*size];
	input_file = fopen(filename, "rb");
	fread(image, 1, size*size, input_file);
	vector<vector<int>> img;
	for (int i = 0; i < size; i++) {
		vector<int> vec;
		for (int j = 0; j < size; j++) {
			vec.push_back(int(image[i * size + j]));
		}img.push_back(vec);
	}
	return img;
}
void write_raw(vector<vector<int>> img,const char filename[],int size) {
	FILE* output_file;
	unsigned char* image = new unsigned char[size*size];
	output_file = fopen(filename, "wb");
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			image[i * size + j] = img[i][j];
		}
	}
	fwrite(image, 1, size*size, output_file);
	fclose(output_file);
}
bool myfunction(int i, int j) {
	return (i < j);
}
void histogram(vector<vector<int>> img, const char filename[], int size) {
	vector<int> plot(256, 0);
	vector<vector<int>> output(256, vector<int>(256, 255));
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			plot[img[i][j]] = plot[img[i][j]] + 1;
		}
	}

	int maxvalue = *max_element(plot.begin(), plot.end());
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			if (i >= 255-(round(255.0 * double(plot[j]) / maxvalue))) {
				output[i][j] = 0;
			}
		}
	}
	write_raw(output, filename, 256);
}
vector<int> to_Binary(int n)
{
	vector < int> d_array;
	if (n / 2 != 0) {
		to_Binary(n / 2);
	}
	d_array.push_back( n % 2);
	return d_array;
}
void Equalization( vector<vector<int>> img, const char filename1[], const char filename2[], int size) {
	vector<int> plot(256, 0);
	vector<vector<int>> output(size, vector<int>(size, 0));
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			plot[img[i][j]] = plot[img[i][j]] + 1;
		}
	}

	int maxvalue = *max_element(plot.begin(), plot.end());
	
	vector<int> vector_1 = plot;

	int sum = 0;
	int EquaHistogram[256] ;
	double cumulative[256] = { 0.0 };
	int cal_sum=0;
	for (int i = 0; i < 256; i++) {
		cal_sum += plot[i];
		cumulative[i] = cal_sum;
	}

	for (int i = 0; i < 256; i++) {
		double x = (cumulative[i] - cumulative[0])/ (size * size -cumulative[0]);
		//double x = (cumulative[i])  / (size * size );

		sum = static_cast<int>((x * 255.0));
		EquaHistogram[i] = sum;
	}
	
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			output[i][j] = round(EquaHistogram[img[i][j]]);
		}
	}
	
	histogram(output, filename1, 512);

	write_raw(output, filename2, 512);

}
void sp_equalization(vector<vector<int>> img, const char filename1[], const char filename2[], int size) {
	vector<int> plot(256, 0);
	vector<vector<int>> output(size, vector<int>(size, 0));
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			plot[img[i][j]] = plot[img[i][j]] + 1;
		}
	}

	vector<int> vector_1 = plot;

	int sum = 0;
	int EquaHistogram[256];
	double cumulative[256] = { 0.0 };
	int cal_sum = 0;

	for (int i = 0; i < 256; i++) {
		cal_sum += plot[i];
		cumulative[i] = cal_sum;
	}

	for (int i = 0; i < 256; i++) {
		double x = (cumulative[i] - cumulative[0]) / (size * size - cumulative[0]);
		sum = static_cast<int>((x * 77.0));
		EquaHistogram[i] = sum;
	}

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			output[i][j] = 77+round(EquaHistogram[img[i][j]]);
		}
	}
	

	histogram(output, filename1, 512);

	write_raw(output, filename2, 512);

}
int main() {
	vector<vector<int>> power_law_512,doreamon,doreamon_256,lena_256,log_512,inverse_log_512,neg_log512,neg_log_512;
	vector<vector<int>> dorlena(256, vector<int>(256, 0));
	doreamon = read_raw("image_file/doreamon.raw", 512);
	for (int i = 0; i < 512; i++) {
		for (int j = 0; j < 512; j++) {
			if (doreamon[i][j] >= 200) {
				doreamon[i][j] = 255;
			}
			else {
				doreamon[i][j] = 0;
			}
		}
	}
	write_raw(doreamon, "image_file/thershold_doreamon.raw", 512);
	doreamon_256 = read_raw("image_file/thershold_doreamon_256.raw", 256);
	lena_256 = read_raw("image_file/lena256.raw", 256);
	for (int cnt = 1; cnt < 9; cnt++) {
		for (int i = 0; i < 256; i++) {
			for (int j = 0; j < 256; j++) {
				int a = lena_256[i][j];
				bitset<8> bit(a);
				bit[cnt-1] = doreamon_256[i][j] / 255;//0 means the most least importance place in img LSB
				a = bit.to_ullong();
				dorlena[i][j] = a;
			}
		}
		char filename[] = "image_file/dorlena_1.raw";
		filename[19] = cnt + '0';
		write_raw(dorlena,filename , 256);
	}

	doreamon_256 = read_raw("image_file/dorlena_1.raw", 256);
	
		for (int i = 0; i < 256; i++) {
			for (int j = 0; j < 256; j++) {
				int a = doreamon_256[i][j];
				bitset<8> bit(a);
				doreamon_256[i][j] = bit[0]*255.0;
			}
		}
	char filename[] = "image_file/dorlena_1_LSB.raw";
	write_raw(doreamon_256, filename, 256);
	//2.1
	log_512 = read_raw("image_file/log512.raw", 512);
	neg_log512 = log_512;
	for (int i = 0; i < 512; ++i) {
		for (int j = 0; j < 512; ++j) {
			neg_log512[i][j] = abs(log_512[i][j] - 255);
		}
	}
	neg_log_512 = neg_log512;
	write_raw(neg_log_512, "image_file/neg_log512.raw", 512);
	double c = 255 / (log(1 + 255));
	for (int i = 0; i <512; ++i) {
		for (int j = 0; j < 512; ++j) {
			log_512[i][j] = c * log(1 + log_512[i][j]);
			neg_log_512[i][j] = c * log(1 + neg_log_512[i][j]);
		}

	}
	write_raw(log_512, "image_file/log_512.raw", 512);
	write_raw(neg_log_512, "image_file/neg_log_512.raw", 512);
	inverse_log_512 = read_raw("image_file/log512.raw", 512);
	neg_log_512 = neg_log512;

	for (int i = 0; i < 512; ++i) {
		for (int j = 0; j < 512; ++j) {
			inverse_log_512[i][j] = pow(exp(inverse_log_512[i][j]), (1 / c)) - 1;
			neg_log_512[i][j] = pow(exp(neg_log_512[i][j]),(1/c)) - 1;
		}
	}
	write_raw(inverse_log_512, "image_file/inverse_log_512.raw", 512);
	write_raw(neg_log_512, "image_file/neg_inverse_log_512.raw", 512);
	power_law_512 = read_raw("image_file/log512.raw", 512);
	neg_log_512 = neg_log512;

	for (int i = 0; i < 512; ++i) {
		for (int j = 0; j < 512; ++j) {
			power_law_512[i][j] = 255 * pow((power_law_512[i][j]/255.0 ),0.2);
			neg_log_512[i][j] = 255 * pow((neg_log_512[i][j] / 255.0), 0.2);
		}
	}
	write_raw(power_law_512, "image_file/power_law_02_512.raw", 512);
	write_raw(neg_log_512, "image_file/neg_power_law_02_512.raw", 512);

	power_law_512 = read_raw("image_file/log512.raw", 512);
	neg_log_512 = neg_log512;

	for (int i = 0; i < 512; ++i) {
		for (int j = 0; j < 512; ++j) {
			power_law_512[i][j] = 255 * pow((power_law_512[i][j] / 255.0), 200);
			neg_log_512[i][j] = 255 * pow((neg_log_512[i][j] / 255.0), 200);
		}
	}
	write_raw(power_law_512, "image_file/power_law_200_512.raw", 512);
	write_raw(neg_log_512, "image_file/neg_power_law_200_512.raw", 512);

	power_law_512 = read_raw("image_file/log512.raw", 512);
	neg_log_512 = neg_log512;

	for (int i = 0; i < 512; ++i) {
		for (int j = 0; j < 512; ++j) {
			power_law_512[i][j] = 255 * pow((power_law_512[i][j] / 255.0), 0.0002);
			neg_log_512[i][j] = 255 * pow((neg_log_512[i][j] / 255.0), 0.0002);
		}
	}
	write_raw(power_law_512, "image_file/power_law_00002_512.raw", 512);
	write_raw(neg_log_512, "image_file/neg_power_law_00002_512.raw", 512);
	power_law_512 = read_raw("image_file/log512.raw", 512);
	neg_log_512 = neg_log512;

	for (int i = 0; i < 512; ++i) {
		for (int j = 0; j < 512; ++j) {
			power_law_512[i][j] = 255 * pow((power_law_512[i][j] / 255.0), 2);
			neg_log_512[i][j] = 255 * pow((neg_log_512[i][j] / 255.0), 2);
		}
	}
	write_raw(power_law_512, "image_file/power_law_20_512.raw", 512);
	write_raw(neg_log_512, "image_file/neg_power_law_20_512.raw", 512);
	log_512 = read_raw("image_file/log512.raw", 512);
	neg_log512 = read_raw("image_file/neg_log512.raw", 512);
	histogram(log_512, "image_file/log_512_histogram.raw", 512);
	histogram(neg_log512, "image_file/neg_log_512_histogram.raw", 512);
	Equalization(log_512, "image_file/log_512_eq_histogram.raw", "image_file/log_512_eq.raw", 512);
	Equalization(neg_log512, "image_file/neg_log_512_eq_histogram.raw", "image_file/neg_log_512_eq.raw", 512);

	log_512 = read_raw("image_file/log512.raw", 512);
	sp_equalization(log_512, "image_file/log_512_speq_histogram.raw", "image_file/log_512_speq.raw", 512);
}