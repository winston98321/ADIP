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

vector<vector<int>> read_raw(const char filename[], int sizew,int sizeh) {
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
void write_raw(vector<vector<int>> img, const char filename[], int sizew , int sizeh) {
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
void box_smooth(vector<vector<int>> img, int filter_size, const char filename[], int sizew,int sizeh) {
	int bund = floor(filter_size / 2);
	vector<vector<int>> out = img;;
	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			int total = 0;
			if (i + bund > sizeh-1 || i - bund<0 || j + bund>sizew-1 || j - bund < 0) {
				//cout << i << " " << j << endl;
				for (int ii = i - bund; ii <= i + bund; ii++) {
					for (int jj = j - bund; jj <= j + bund; jj++) {
						int iii = ii, jjj = jj;
						if (ii < 0) iii = abs(ii)-1;
						if (jj < 0) jjj = abs(jj)-1;
						if (ii > sizeh-1) iii = sizeh-1 + (sizeh-1 - ii+1);
						if (jj > sizew-1) jjj = sizew-1 + (sizew-1 - jj+1);
						//cout << ii << " " << jj << " " << iii <<" " << jjj << endl;
						total += img[iii][jjj];
					}
				}
			}
			else {
				for (int ii = i - bund; ii <= i + bund; ii++) {
					for (int jj = j - bund; jj <= j + bund; jj++) {
						total += img[ii][jj];
					}
				}
			}
			//cout << total << " " << (filter_size * filter_size) << endl;
			out[i][j] = total / (filter_size * filter_size);
		}
	}
	write_raw(out, filename, sizew, sizeh);
}
void gaussian_smooth(vector<vector<int>> img, int filter_size, const char filename[], int sizew, int sizeh) {
	int bund = floor(filter_size / 2);
	vector<vector<int>> out = img;
	vector<int> gaussian_kernel_3x3 = { 16,1,2,1,2,4,2,1,2,1 };
	vector<int> gaussian_kernel_5x5 = { 273,1,4,7,4,1,4,16,26,16,4,7,26,41,26,7,4,16,26,16,4,1,4,7,4,1 };
	vector<int> kernel = gaussian_kernel_3x3;
	if (filter_size == 5) {
		kernel = gaussian_kernel_5x5;
	}
	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			double total = 0;
			int kernel_cnt=0;
			if (i + bund > sizeh - 1 || i - bund<0 || j + bund>sizew - 1 || j - bund < 0) {
				for (int ii = i - bund; ii <= i + bund; ii++) {
					for (int jj = j - bund; jj <= j + bund; jj++) {
						kernel_cnt += 1;
						int iii = ii, jjj = jj;
						if (ii < 0) iii = abs(ii)-1;
						if (jj < 0) jjj = abs(jj)-1;
						if (ii > sizeh - 1) iii = sizeh - 1 + (sizeh - 1 - ii+1);
						if (jj > sizew - 1) jjj = sizew - 1 + (sizew - 1 - jj+1);
						total += (double(kernel[kernel_cnt])/kernel[0])*img[iii][jjj];
					}
				}
			}
			else {
				for (int ii = i - bund; ii <= i + bund; ii++) {
					for (int jj = j - bund; jj <= j + bund; jj++) {
						kernel_cnt += 1;
						total += (double(kernel[kernel_cnt]) / kernel[0]) * img[ii][jj];
					}
				}
			}
			out[i][j] = total;
		}
	}
	write_raw(out, filename, sizew, sizeh);
}
void robert(vector<vector<int>> img, int kernel_angle, const char filename[], int sizew, int sizeh) {
	int bund = 1;
	vector<vector<int>> out = img;
	vector<int> robert_kernel_45 = { 0,0,0,0,0,-1,0,0,0,1 };
	vector<int> robert_kernel_n45 = { 0,0,0,0,0,0,-1,0,1,0};
	vector<int> kernel = robert_kernel_45;
	if (kernel_angle == -45) {
		kernel = robert_kernel_n45;
	}
	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			double total = 0;
			double combine45 = 0;
			double combinen45 = 0;
			int kernel_cnt = 0;
			if (i + bund > sizeh - 1 || i - bund<0 || j + bund>sizew - 1 || j - bund < 0) {
				for (int ii = i - bund; ii <= i + bund; ii++) {
					for (int jj = j - bund; jj <= j + bund; jj++) {
						kernel_cnt += 1;
						int iii = ii, jjj = jj;
						if (ii < 0) iii = abs(ii)-1;
						if (jj < 0) jjj = abs(jj)-1;
						if (ii > sizeh - 1) iii = sizeh - 1 + (sizeh - 1 - ii+1);
						if (jj > sizew - 1) jjj = sizew - 1 + (sizew - 1 - jj+1);
						total += double(kernel[kernel_cnt])* img[iii][jjj];
						combine45 += double(robert_kernel_45[kernel_cnt]) * img[iii][jjj];
						combinen45 += double(robert_kernel_n45[kernel_cnt]) * img[iii][jjj];
					}
				}
			}
			else {
				for (int ii = i - bund; ii <= i + bund; ii++) {
					for (int jj = j - bund; jj <= j + bund; jj++) {
						kernel_cnt += 1;
						total += double(kernel[kernel_cnt]) * img[ii][jj];
						combine45 += double(robert_kernel_45[kernel_cnt]) * img[ii][jj];
						combinen45 += double(robert_kernel_n45[kernel_cnt]) * img[ii][jj];
					}
				}
			}
			out[i][j] = total;
			if (kernel_angle == 90) out[i][j] = pow((pow(combine45,2)+pow(combinen45,2)),0.5);
			if (out[i][j] > 255)out[i][j] = 255;
			if (out[i][j] < 0)out[i][j] = 0;
		}
	}
	write_raw(out, filename, sizew, sizeh);
}
void prewitt(vector<vector<int>> img, int kernel_angle, const char filename[], int sizew, int sizeh) {
	int bund = 1;
	vector<vector<int>> out = img;
	vector<int> prewitt_kernel_0 = { 0,1,0,-1,1,0,-1,1,0,-1 };
	vector<int> prewitt_kernel_90 = { 0,-1,-1,-1,0,0,0,1,1,1 };
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
			int cb = abs(combine0+combine90);
			
			out[i][j] = total;
			if (kernel_angle == 1) out[i][j] = cb;
			if (out[i][j] > 255)out[i][j] = 255;
			if (out[i][j] < 0)out[i][j] = 0;
		}
	}
	/*
	double max = -999;
	for (int i = 0; i < 421;i++) {
		for (int j = 0; j < 700; j++) {
			if (max < out[i][j]) max = out[i][j];
		}
	}
	for (int i = 0; i < 421; i++) {
		for (int j = 0; j < 700; j++) {
			out[i][j] = ((out[i][j]) / max) * 255;
		}
	}*/
	write_raw(out, filename, sizew, sizeh);
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
	/*
	double max = -999;
	for (int i = 0; i < 421;i++) {
		for (int j = 0; j < 700; j++) {
			if (max < out[i][j]) max = out[i][j];
		}
	}
	for (int i = 0; i < 421; i++) {
		for (int j = 0; j < 700; j++) {
			out[i][j] = ((out[i][j]) / max) * 255;
		}
	}*/
	write_raw(out, filename, sizew, sizeh);
}
void sobel45(vector<vector<int>> img, int kernel_angle, const char filename[], int sizew, int sizeh) {
	int bund = 1;
	vector<vector<int>> out = img;
	vector<int> prewitt_kernel_0 = { 0,-2,-1,0,-1,0,1,0,1,2 };
	vector<int> prewitt_kernel_90 = { 0,0,-1,-2,1,0,-1,2,1,0 };
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
	/*
	double max = -999;
	for (int i = 0; i < 421;i++) {
		for (int j = 0; j < 700; j++) {
			if (max < out[i][j]) max = out[i][j];
		}
	}
	for (int i = 0; i < 421; i++) {
		for (int j = 0; j < 700; j++) {
			out[i][j] = ((out[i][j]) / max) * 255;
		}
	}*/
	write_raw(out, filename, sizew, sizeh);
}
void laplacian(vector<vector<int>> img, int kernel_angle, const char filename[], int sizew, int sizeh) {
	int bund = 1;
	vector<vector<int>> out = img;
	vector<int> prewitt_kernel_0 = { 0,0,-1,0,-1,4,-1,0,-1,0 };
	vector<int> prewitt_kernel_90 = { 0,-1,-1,-1,-1,8,-1,-1,-1,-1 };
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
	/*
	double max = -999;
	for (int i = 0; i < 421;i++) {
		for (int j = 0; j < 700; j++) {
			if (max < out[i][j]) max = out[i][j];
		}
	}
	for (int i = 0; i < 421; i++) {
		for (int j = 0; j < 700; j++) {
			out[i][j] = ((out[i][j]) / max) * 255;
		}
	}*/
	write_raw(out, filename, sizew, sizeh);
}
void high_boost_box(vector<vector<int>> img, int filter_size, const char filename[], int sizew, int sizeh) {
	int bund = floor(filter_size / 2);
	vector<vector<int>> out = img;
	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			int total = 0;
			if (i + bund > sizeh - 1 || i - bund<0 || j + bund>sizew - 1 || j - bund < 0) {
				for (int ii = i - bund; ii <= i + bund; ii++) {
					for (int jj = j - bund; jj <= j + bund; jj++) {
						int iii = ii, jjj = jj;
						if (ii < 0) iii = abs(ii) - 1;
						if (jj < 0) jjj = abs(jj) - 1;
						if (ii > sizeh - 1) iii = sizeh - 1 + (sizeh - 1 - ii + 1);
						if (jj > sizew - 1) jjj = sizew - 1 + (sizew - 1 - jj + 1);
						total += img[iii][jjj];
					}
				}
			}
			else {
				for (int ii = i - bund; ii <= i + bund; ii++) {
					for (int jj = j - bund; jj <= j + bund; jj++) {
						total += img[ii][jj];
					}
				}
			}
			out[i][j] = 2*img[i][j] - total / (filter_size * filter_size);
			//cout << out[i][j] << endl;
			if (out[i][j] < 0) { out[i][j] = 0; }
			if (out[i][j] > 255) { out[i][j] = 255; }
		}
	}
	write_raw(out, filename, sizew, sizeh);
}
void high_boost_gau(vector<vector<int>> img, int filter_size, const char filename[], int sizew, int sizeh) {
	int bund = floor(filter_size / 2);
	vector<vector<int>> out = img;
	vector<int> gaussian_kernel_3x3 = { 16,1,2,1,2,4,2,1,2,1 };
	vector<int> gaussian_kernel_5x5 = { 273,1,4,7,4,1,4,16,26,16,4,7,26,41,26,7,4,16,26,16,4,1,4,7,4,1 };
	vector<int> kernel = gaussian_kernel_3x3;
	if (filter_size == 5) {
		kernel = gaussian_kernel_5x5;
	}
	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			double total = 0;
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
						total += (double(kernel[kernel_cnt]) / kernel[0]) * img[iii][jjj];
					}
				}
			}
			else {
				for (int ii = i - bund; ii <= i + bund; ii++) {
					for (int jj = j - bund; jj <= j + bund; jj++) {
						kernel_cnt += 1;
						total += (double(kernel[kernel_cnt]) / kernel[0]) * img[ii][jj];
					}
				}
			}
			out[i][j] = 2*img[i][j] -total;
			if (out[i][j] < 0) { out[i][j] = 0; }
			if (out[i][j] > 255) { out[i][j] = 255; }
		}
	}
	write_raw(out, filename, sizew, sizeh);
}


void gaussian_threshold(vector<vector<int>> img, int filter_size, const char filename[], int sizew, int sizeh) {
	int bund = floor(filter_size / 2);
	vector<vector<int>> out = img;
	vector<int> gaussian_kernel_3x3 = { 16,1,2,1,2,4,2,1,2,1 };
	vector<int> gaussian_kernel_5x5 = { 273,1,4,7,4,1,4,16,26,16,4,7,26,41,26,7,4,16,26,16,4,1,4,7,4,1 };
	vector<int> kernel = gaussian_kernel_3x3;
	if (filter_size == 5) {
		kernel = gaussian_kernel_5x5;
	}
	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			double total = 0;
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
						total += (double(kernel[kernel_cnt]) / kernel[0]) * img[iii][jjj];
					}
				}
			}
			else {
				for (int ii = i - bund; ii <= i + bund; ii++) {
					for (int jj = j - bund; jj <= j + bund; jj++) {
						kernel_cnt += 1;
						total += (double(kernel[kernel_cnt]) / kernel[0]) * img[ii][jj];
					}
				}
			}
			if (total >190) total = 255;
			else total = 0;
			out[i][j] = total;
		}
		
	}
	write_raw(out, filename, sizew, sizeh);

	}

void box_threshold(vector<vector<int>> img, int filter_size, const char filename[], int sizew, int sizeh) {
	int bund = floor(filter_size / 2);
	vector<vector<int>> out = img;;
	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			int total = 0;
			if (i + bund > sizeh - 1 || i - bund<0 || j + bund>sizew - 1 || j - bund < 0) {
				//cout << i << " " << j << endl;
				for (int ii = i - bund; ii <= i + bund; ii++) {
					for (int jj = j - bund; jj <= j + bund; jj++) {
						int iii = ii, jjj = jj;
						if (ii < 0) iii = abs(ii) - 1;
						if (jj < 0) jjj = abs(jj) - 1;
						if (ii > sizeh - 1) iii = sizeh - 1 + (sizeh - 1 - ii + 1);
						if (jj > sizew - 1) jjj = sizew - 1 + (sizew - 1 - jj + 1);
						//cout << ii << " " << jj << " " << iii <<" " << jjj << endl;
						total += img[iii][jjj];
					}
				}
			}
			else {
				for (int ii = i - bund; ii <= i + bund; ii++) {
					for (int jj = j - bund; jj <= j + bund; jj++) {
						total += img[ii][jj];
					}
				}
			}
			//cout << total << " " << (filter_size * filter_size) << endl;
			total = total / (filter_size * filter_size);
			if (total > 190)total = 255;
			else total = 0;
			out[i][j] = total;
		}
	}
	write_raw(out, filename, sizew, sizeh);
}
void gaussian_logic(vector<vector<int>> img, vector<vector<int>> logic_filter, int filter_size, const char filename[], int sizew, int sizeh) {
	int bund = floor(filter_size / 2);
	vector<vector<int>> out = img;
	vector<int> gaussian_kernel_3x3 = { 16,1,2,1,2,4,2,1,2,1 };
	vector<int> gaussian_kernel_5x5 = { 273,1,4,7,4,1,4,16,26,16,4,7,26,41,26,7,4,16,26,16,4,1,4,7,4,1 };
	vector<int> kernel = gaussian_kernel_3x3;
	if (filter_size == 5) {
		kernel = gaussian_kernel_5x5;
	}
	for (int i = 0; i < sizeh; i++) {
		for (int j = 0; j < sizew; j++) {
			double total = 0;
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
						total += (double(kernel[kernel_cnt]) / kernel[0]) * img[iii][jjj];
					}
				}
			}
			else {
				for (int ii = i - bund; ii <= i + bund; ii++) {
					for (int jj = j - bund; jj <= j + bund; jj++) {
						kernel_cnt += 1;
						total += (double(kernel[kernel_cnt]) / kernel[0]) * img[ii][jj];
					}
				}
			}
			if (logic_filter[i][j] == 0) {
				out[i][j] = total;
			}
			else {
				out[i][j] = 2*out[i][j] - total;
				if (out[i][j] < 0) out[i][j] = 0;
				if (out[i][j] > 255) out[i][j] = 255;
			}
		}
	}
	write_raw(out, filename, sizew, sizeh);
}
void main() {
	vector<vector<vector<int>>> shape100(9, vector<vector<int> >(100, vector<int>(100, 0))), shape100_out(9, vector<vector<int> >(100, vector<int>(100, 0)));
	vector<vector<int>> slate, slate_noise, flower, logic_filter;
	//Q1
	for (int i = 1; i < 9; i++) {
		char filename[] = "src/shape100_1.raw";
		filename[13] = i + '0';
		shape100[i] = read_raw(filename, 100, 100);
	}
	for (int k = 1; k < 9; k++) {
		double x = 0, y = 0;
		double m00 = 0, m01 = 0, m10 = 0;
		double miu = 0;
		for (int i = 0; i < 100; i++) {
			for (int j = 0; j < 100; j++) {
				m00 += shape100[k][i][j];
				m01 += j * shape100[k][i][j];
				m10 += i * shape100[k][i][j];
			}
		}
		x = m10 / m00; y = m01 / m00;
		cout << "1.1 centroid of " << "shape100_" << k << " x:" << x << " y:" << y << endl;
		double p10 = 0, p01 = 0, p20 = 0, p02 = 0,p11=0, p21 = 0, p12 = 0,p30=0,p03=0;
		for (int i = 0; i < 100; i++) {
			for (int j = 0; j < 100; j++) {
				p10 += pow((i - x), 1) * pow((j - y), 0) * shape100[k][i][j];
				p01 += pow((i - x), 0) * pow((j - y), 1) * shape100[k][i][j];
				p20 += pow((i - x), 2) * pow((j - y), 0) * shape100[k][i][j];
				p02 += pow((i - x), 0) * pow((j - y), 2) * shape100[k][i][j];
				p11 += pow((i - x), 1) * pow((j - y), 1) * shape100[k][i][j];
				p21 += pow((i - x), 2) * pow((j - y), 1) * shape100[k][i][j];
				p12 += pow((i - x), 1) * pow((j - y), 2) * shape100[k][i][j];
				p03 += pow((i - x), 0) * pow((j - y), 3) * shape100[k][i][j];
				p30 += pow((i - x), 3) * pow((j - y), 0) * shape100[k][i][j];
			}
		}
		cout << "1.2 central momen order " << "p=1,q=0 :" << p10 <<
			" p=0,q=1 :" << p01 <<
			" p=2,q=0 :" << p20 <<
			" p=0,q=2 :" << p02 <<
			" p=1,q=1 :" << p11 <<
			" p=2,q=1 :" << p21 <<
			" p=1,q=2 :" << p12 << 
			" p=0,q=3 :" << p03 << 
			" p=3,q=0 :" << p30 << endl;
	}
	//Q2.1 boxsmooth and gaussiansblur

	slate = read_raw("src/slate_700x421.raw", 700, 421);
	slate_noise = read_raw("src/slate_noise_700x421.raw", 700, 421);
	box_smooth(slate, 3, "image_file/box_smooth_3x3_slate.raw", 700, 421);
	box_smooth(slate, 5, "image_file/box_smooth_5x5_slate.raw", 700, 421);
	box_smooth(slate_noise, 3, "image_file/box_smooth_3x3_slate_noise.raw", 700, 421);
	box_smooth(slate_noise, 5, "image_file/box_smooth_5x5_slate_noise.raw", 700, 421);
	gaussian_smooth(slate, 3, "image_file/gaussian_smooth_3x3_slate.raw", 700, 421);
	gaussian_smooth(slate, 5, "image_file/gaussian_smooth_5x5_slate.raw", 700, 421);
	gaussian_smooth(slate_noise, 3, "image_file/gaussian_smooth_3x3_slate_noise.raw", 700, 421);
	gaussian_smooth(slate_noise, 5, "image_file/gaussian_smooth_5x5_slate_noise.raw", 700, 421);

	//Q2.2

	robert(slate, 45, "image_file/robert_45_slate.raw", 700, 421);
	robert(slate, -45, "image_file/robert_n45_slate.raw", 700, 421);
	robert(slate, 90, "image_file/robert_combine_slate.raw", 700, 421);
	robert(slate_noise, 45, "image_file/robert_45_slate_noise.raw", 700, 421);
	robert(slate_noise, -45, "image_file/robert_n45_slate_noise.raw", 700, 421);
	robert(slate_noise, 90, "image_file/robert_combine_slate_noise.raw", 700, 421);

	//Q2.3
	prewitt(slate, 0, "image_file/prewitt_0_slate.raw", 700, 421);
	prewitt(slate, 90, "image_file/prewitt_90_slate.raw", 700, 421);
	prewitt(slate, 1, "image_file/prewitt_combine_slate.raw", 700, 421);
	prewitt(slate_noise, 0, "image_file/prewitt_0_slate_noise.raw", 700, 421);
	prewitt(slate_noise, 90, "image_file/prewitt_90_slate_noise.raw", 700, 421);
	prewitt(slate_noise, 1, "image_file/prewitt_combine_slate_noise.raw", 700, 421);

	//Q2.4
	sobel(slate, 0, "image_file/sobel_0_slate.raw", 700, 421);
	sobel(slate, 90, "image_file/sobel_90_slate.raw", 700, 421);
	sobel(slate, 1, "image_file/sobel_combine_slate.raw", 700, 421);
	sobel(slate_noise, 0, "image_file/sobel_0_slate_noise.raw", 700, 421);
	sobel(slate_noise, 90, "image_file/sobel_90_slate_noise.raw", 700, 421);
	sobel(slate_noise, 1, "image_file/sobel_combine_slate_noise.raw", 700, 421);

	//Q2.5
	sobel45(slate, 0, "image_file/sobel_45_slate.raw", 700, 421);
	sobel45(slate, 90, "image_file/sobel_n45_slate.raw", 700, 421);
	sobel45(slate, 1, "image_file/sobel45_combine_slate.raw", 700, 421);
	sobel45(slate_noise, 0, "image_file/sobel_45_slate_noise.raw", 700, 421);
	sobel45(slate_noise, 90, "image_file/sobel_n45_slate_noise.raw", 700, 421);
	sobel45(slate_noise, 1, "image_file/sobel45_combine_slate_noise.raw", 700, 421);

	//Q2.6 laplacian
	laplacian(slate, 0, "image_file/laplacian_4_slate.raw", 700, 421);
	laplacian(slate, 90, "image_file/laplacian_8_slate.raw", 700, 421);
	laplacian(slate_noise, 0, "image_file/laplacian_4_slate_noise.raw", 700, 421);
	laplacian(slate_noise, 90, "image_file/laplacian_8_slate_noise.raw", 700, 421);

	//Q2.8 high-boost = org - low pass filter
	high_boost_box(slate, 3, "image_file/highboost_box_slate_.raw", 700, 421);
	high_boost_gau(slate, 3, "image_file/highboost_gau_slate_.raw", 700, 421);
	high_boost_box(slate_noise, 3, "image_file/highboost_box_slate_noise_.raw", 700, 421);
	high_boost_gau(slate_noise, 3, "image_file/highboost_gau_slate_noise_.raw", 700, 421);

	//Q3.1 blur thresholding

	flower = read_raw("src/flower_512x384.raw", 512, 384);
	gaussian_threshold(flower, 5, "image_file/blur_thresholding_flower.raw", 512, 384);

	//Q3.2 logic blur
	logic_filter = read_raw("image_file/blur_thresholding_flower.raw", 512, 384);
	gaussian_logic(flower, logic_filter, 5, "image_file/logic_blur_flower.raw", 512, 384);

}