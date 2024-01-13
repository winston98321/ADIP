//#define _CRT_SECURE_NO_DEPRECATE
//
//#include <stdio.h>
//#include <iostream>
//#include <vector>
//#include <algorithm>
//#include <cmath>
//#include <stdlib.h>
//#include <time.h>
//#include <queue>
//#include <complex>
//
//
//#define PI 3.141592653589793
//using namespace std;
//clock_t a, b;
//vector<vector<int>> read_raw(const char filename[], int sizew, int sizeh) {
//	FILE* input_file;
//	unsigned char* image = new unsigned char[sizew * sizeh];
//	input_file = fopen(filename, "rb");
//	fread(image, 1, sizew * sizeh, input_file);
//	vector<vector<int>> img;
//	for (int i = 0; i < sizeh; i++) {
//		vector<int> vec;
//		for (int j = 0; j < sizew; j++) {
//			vec.push_back(int(image[i * sizew + j]));
//		}img.push_back(vec);
//	}
//	return img;
//}
//void write_raw(vector<vector<int>> img, const char filename[], int sizew, int sizeh) {
//	FILE* output_file;
//	unsigned char* image = new unsigned char[sizew * sizeh];
//	output_file = fopen(filename, "wb");
//	for (int i = 0; i < sizeh; i++) {
//		for (int j = 0; j < sizew; j++) {
//			image[i * sizew + j] = img[i][j];
//		}
//	}
//	fwrite(image, 1, sizeh * sizew, output_file);
//	fclose(output_file);
//}
//
//vector<vector<complex<double>>> dft(vector<vector<int>> img, const char filename[], int sizew, int sizeh) {
//	vector<vector<complex<double>>> f1_out(sizeh, vector<complex<double>>(sizew));
//	vector<vector<complex<double>>> f2_out(sizeh, vector<complex<double>>(sizew));
//	vector<vector<int>> img_out(sizeh, vector<int>(sizew, 0));
//	auto M_I_2PI_DL = -(6.28318530718i / double(sizew));
//	a = clock();
//	for (int i = 0; i < sizeh; i++) {
//		for (int j = 0; j < sizew; j++) {
//			complex<double> sum = (0,0);
//			for (int k = 0; k < sizew; k++) {
//				sum += double(img[i][k]) *exp(M_I_2PI_DL * double(k) * double(j)) * exp(M_I_2PI_DL * double(k) * double(128));
//			}
//			f1_out[i][j] = sum;
//		}
//	}
//
//
//	for (int i = 0; i < sizeh; i++) {
//		for (int j = 0; j < sizew; j++) {
//			complex<double> sum = (0,0);
//			for (int k = 0; k < sizew; k++) {
//				sum += f1_out[k][i] * std::exp(M_I_2PI_DL * double(k) * double(j)) * exp(M_I_2PI_DL * double(k) * double(128));//
//			}
//			f2_out[j][i] = sum;
//		}
//	}
//
//	for (int i = 0; i < sizeh; i++) {
//		for (int j = 0; j < sizew; j++) {
//			//cout << f2_out[i][j].real() << " " << f2_out[i][j].imag() << endl;
//			double ans = 255 * (round(sqrt(f2_out[i][j].real() * f2_out[i][j].real() + f2_out[i][j].imag() * f2_out[i][j].imag())) / (256.0 * 256));
//			//cout << ans << endl;
//			//if (ans > 255)ans = 255;
//			//if (ans < 0) ans = 0;
//			img_out[i][j] = ans;
//		}
//	}
//	img = img_out;
//	//log transfer for image enhancement
//	double c = 255 / (log(1 + 255));
//	for (int i = 0; i < sizeh; ++i) {
//		for (int j = 0; j < sizew; ++j) {
//			img_out[i][j] = c * log(1 + img_out[i][j]);
//		}
//
//	}
//	double ma = -999;
//	for (int i = 0; i < sizeh; ++i) {
//		for (int j = 0; j < sizew; ++j) {
//			if (img_out[i][j] > ma) {
//				ma = img_out[i][j];
//			}
//		}
//
//	}
//	for (int i = 0; i < sizeh; ++i) {
//		for (int j = 0; j < sizew; ++j) {
//			img_out[i][j] = 255 * (img_out[i][j] / ma);
//			//cout << ma<<endl;
//		}
//
//	}
//	b = clock();
//	cout << "DFT Caluclate time :" << double(b - a) / CLOCKS_PER_SEC << endl;
//	write_raw(img_out, filename, sizew, sizeh);
//	return f2_out;
//}
//void idft(vector<vector<complex<double>>> img, const char filename[], int sizew, int sizeh)
//{
//	std::vector<std::vector<complex<double>>> f1_out(sizeh, std::vector<complex<double>>(sizew));
//	std::vector<std::vector<complex<double>>> f2_out(sizeh, std::vector<complex<double>>(sizew));
//	std::vector<std::vector<int>> img_out(sizeh, std::vector<int>(sizew, 0));
//	auto M_I_2PI_DL = (6.28318530718i / double(sizew));
//	a = clock();
//	for (int i = 0; i < sizeh; i++) {
//		for (int j = 0; j < sizew; j++) {
//			complex<double> sum = (0, 0);
//			for (int k = 0; k < sizew; k++) {
//				sum += img[i][k] * exp(M_I_2PI_DL * double(k) * double(j));// * exp(M_I_2PI_DL * double(k) * double(128));
//			}
//			f1_out[i][j] = sum;
//		}
//	}
//
//
//
//	for (int i = 0; i < sizeh; i++) {
//		for (int j = 0; j < sizew; j++) {
//			complex<double> sum = (0, 0);
//			for (int k = 0; k < sizew; k++) {
//				sum += f1_out[k][i] * exp(M_I_2PI_DL * double(k) * double(j));// *exp(M_I_2PI_DL * double(k) * double(128));
//			}
//			f2_out[j][i] = sum;
//		}
//	}
//
//	for (int i = 0; i < sizeh; i++) {
//		for (int j = 0; j < sizew; j++) {
//			//cout << f2_out[i][j].real() << " " << f2_out[i][j].imag() << endl;
//			int ans = round(sqrt(f2_out[i][j].real() * f2_out[i][j].real() + f2_out[i][j].imag() * f2_out[i][j].imag())) / (256 * 256);
//			img_out[i][j] = ans;
//		}
//	}
//	b = clock();
//	cout << "IDFT Caluclate time :" << double(b - a) / CLOCKS_PER_SEC << endl;
//	write_raw(img_out, filename, sizew, sizeh);
//}
//void compare(vector<vector<int>> vec, vector<vector<int>> src, string a) {
//	double MSE = 0;
//	int len = vec[0].size();
//	for (int i = 0; i < len; i++) {
//		for (int j = 0; j < len; j++) {
//			MSE += pow((src[i][j] - vec[i][j]), 2);
//		}
//	}
//	MSE = (MSE / (len * len));
//	cout << a << " MSE :" << MSE << endl;
//	cout << a << " PSNR :" << 10 * log10((255 * 255) / MSE) << endl;
//}
//void ALNRF(vector<vector<int>> img, int kernelsize,const char filename[], int sizew, int sizeh) {
//	int bund = kernelsize/2;
//	vector<vector<int>> out = img;
//	double miu = 0,cov_all=0;
//	for (int i = 0; i < sizeh; i++) {
//		for (int j = 0; j < sizew; j++){
//			miu += img[i][j];
//		}
//	}miu /= (sizeh * sizew);
//	
//	for (int i = 0; i < sizeh; i++) {
//		for (int j = 0; j < sizew; j++) {
//			cov_all += pow(img[i][j] - miu, 2);
//		}
//	}cov_all /= (sizeh * sizew);
//
//	
//	for (int i = 0; i < sizeh; i++) {
//		for (int j = 0; j < sizew; j++) {
//			double total = 0;
//			double filter_miu = 0;
//			double filter_cov = 0;
//			double cov_total=0;
//			if (i + bund > sizeh - 1 || i - bund<0 || j + bund>sizew - 1 || j - bund < 0) {
//				for (int ii = i - bund; ii <= i + bund; ii++) {
//					for (int jj = j - bund; jj <= j + bund; jj++) {
//						int iii = ii, jjj = jj;
//						if (ii < 0) iii = abs(ii) - 1;
//						if (jj < 0) jjj = abs(jj) - 1;
//						if (ii > sizeh - 1) iii = sizeh - 1 + (sizeh - 1 - ii + 1);
//						if (jj > sizew - 1) jjj = sizew - 1 + (sizew - 1 - jj + 1);
//						total +=  img[iii][jjj];
//						cov_total += pow(img[iii][jjj],2);
//					}
//				}
//			}
//			else {
//				for (int ii = i - bund; ii <= i + bund; ii++) {
//					for (int jj = j - bund; jj <= j + bund; jj++) {
//						total += img[ii][jj];
//						cov_total += pow(img[ii][jj], 2);
//					}
//				}
//			}
//			filter_miu =total/ (kernelsize * kernelsize);
//			filter_cov = cov_total / (kernelsize * kernelsize) - pow(filter_miu, 2);
//			double cnt = (cov_all / filter_cov);
//			//cout << cnt << endl;
//			if (cnt > 1)cnt = 1;
//			out[i][j] = img[i][j] -cnt*(img[i][j] - filter_miu);
//			//cout << out[i][j] << endl;
//			if (out[i][j] > 255)out[i][j] = 255;
//			if (out[i][j] < 0)out[i][j] = 0;
//		}
//	}
//
//	write_raw(out, filename, sizew, sizeh);
//}
//void ATMF(vector<vector<int>> img, int kernelsize,int d, const char filename[], int sizew, int sizeh) {
//	int bund = kernelsize / 2;
//	vector<vector<int>> out = img;
//
//	for (int i = 0; i < sizeh; i++) {
//		for (int j = 0; j < sizew; j++) {
//			double total = 0;
//			vector<int> record_filter;
//			double cov_total = 0;
//			if (i + bund > sizeh - 1 || i - bund<0 || j + bund>sizew - 1 || j - bund < 0) {
//				for (int ii = i - bund; ii <= i + bund; ii++) {
//					for (int jj = j - bund; jj <= j + bund; jj++) {
//						int iii = ii, jjj = jj;
//						if (ii < 0) iii = abs(ii) - 1;
//						if (jj < 0) jjj = abs(jj) - 1;
//						if (ii > sizeh - 1) iii = sizeh - 1 + (sizeh - 1 - ii + 1);
//						if (jj > sizew - 1) jjj = sizew - 1 + (sizew - 1 - jj + 1);
//						
//						record_filter.push_back(img[iii][jjj]);
//						
//					}
//				}
//			}
//			else {
//				for (int ii = i - bund; ii <= i + bund; ii++) {
//					for (int jj = j - bund; jj <= j + bund; jj++) {
//						
//						record_filter.push_back(img[ii][jj]);
//					}
//				}
//			}
//			sort(record_filter.begin(), record_filter.end());
//			for (int q = d / 2 ; q < kernelsize * kernelsize -(d / 2 );q++) {
//				//cout << q << " " << record_filter[q] << endl;
//				total += record_filter[q];
//			}
//			//cout << total << endl;
//			total = total/((kernelsize * kernelsize) - d);
//			
//			out[i][j] = total;
//			//cout << out[i][j] << endl;
//			if (out[i][j] > 255)out[i][j] = 255;
//			if (out[i][j] < 0)out[i][j] = 0;
//		}
//	}
//
//	write_raw(out, filename, sizew, sizeh);
//}
//int main()
//{
//	vector<vector<int>>lena_saltandpepper, lena_gussuian_1, lena_gussian_2;
//	vector<vector<int>>lena_256_blur, lena_256_blur_noise;
//	lena_256_blur = read_raw("src/lena_256_blur.raw", 256, 256);
//	lena_256_blur_noise = read_raw("src/lena_256_blur_noise.raw", 256, 256);
//	lena_saltandpepper = read_raw("src/lean_256_salt&pepper.raw", 256, 256);
//	lena_gussuian_1 = read_raw("src/lena_gussuian_1_256.raw", 256, 256);
//	lena_gussian_2 = read_raw("src/lena_gussuian_2_256.raw", 256, 256);
//
//	ALNRF(lena_saltandpepper,3, "image_file/1/ALNRF_lean_256_salt&pepper_3x3.raw", 256, 256);
//	ATMF(lena_saltandpepper, 3,2 ,"image_file/1/ATMF_lean_256_salt&pepper_3x3_2.raw", 256, 256);
//	ATMF(lena_saltandpepper, 3,4, "image_file/1/ATMF_lean_256_salt&pepper_3x3_4.raw", 256, 256);
//	ATMF(lena_saltandpepper, 3,8, "image_file/1/ATMF_lean_256_salt&pepper_3x3_8.raw", 256, 256);
//
//
//
//	ALNRF(lena_gussuian_1, 3, "image_file/1/ALNRF_lena_gussuian_1.raw_3x3.raw", 256, 256);
//	ATMF(lena_gussuian_1, 3,2 ,"image_file/1/ATMF_lena_gussuian_1.raw_3x3_2.raw", 256, 256);
//	ATMF(lena_gussuian_1, 3,4, "image_file/1/ATMF_lena_gussuian_1.raw_3x3_4.raw", 256, 256);
//	ATMF(lena_gussuian_1, 3,8, "image_file/1/ATMF_lena_gussuian_1.raw_3x3_8.raw", 256, 256);
//
//	
//	ALNRF(lena_gussian_2, 3, "image_file/1/ALNRF_lena_gussuian_2.raw_3x3.raw", 256, 256);
//	ATMF(lena_gussian_2, 3,2, "image_file/1/ATMF_lena_gussuian_2.raw_3x3_2.raw", 256, 256);
//	ATMF(lena_gussian_2, 3,4, "image_file/1/ATMF_lena_gussuian_2.raw_3x3_4.raw", 256, 256);
//	ATMF(lena_gussian_2, 3,8, "image_file/1/ATMF_lena_gussuian_2.raw_3x3_8.raw", 256, 256);
//
//	
//
//	ALNRF(lena_saltandpepper, 9, "image_file/1/ALNRF_lena_256_salt&pepper_9x9.raw", 256, 256);
//	ATMF(lena_saltandpepper, 9, 8, "image_file/1/ATMF_lena_256_salt&pepper_9x9_8.raw", 256, 256);
//	ATMF(lena_saltandpepper, 9, 16, "image_file/1/ATMF_lena_256_salt&pepper_9x9_16.raw", 256, 256);
//	ATMF(lena_saltandpepper, 9, 32, "image_file/1/ATMF_lena_256_salt&pepper_9x9_32.raw", 256, 256);
//
//	ALNRF(lena_gussuian_1, 9, "image_file/1/ALNRF_lena_gussuian_1.raw_9x9.raw", 256, 256);
//	ATMF(lena_gussuian_1, 9, 8, "image_file/1/ATMF_lena_gussuian_1.raw_9x9_8.raw", 256, 256);
//	ATMF(lena_gussuian_1, 9, 16, "image_file/1/ATMF_lena_gussuian_1.raw_9x9_16.raw", 256, 256);
//	ATMF(lena_gussuian_1, 9, 32, "image_file/1/ATMF_lena_gussuian_1.raw_9x9_32.raw", 256, 256);
//
//
//	ALNRF(lena_gussian_2, 9, "image_file/1/ALNRF_lena_gussuian_2.raw_9x9.raw", 256, 256);
//	ATMF(lena_gussian_2, 9, 8, "image_file/1/ATMF_lena_gussuian_2.raw_9x9_8.raw", 256, 256);
//	ATMF(lena_gussian_2, 9, 16, "image_file/1/ATMF_lena_gussuian_2.raw_9x9_16.raw", 256, 256);
//	ATMF(lena_gussian_2, 9, 32, "image_file/1/ATMF_lena_gussuian_2.raw_9x9_32.raw", 256, 256);
//
//
//
//	return 0;
//}
