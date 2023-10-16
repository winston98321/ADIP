/********************************************************
* Filename    : hw2.cpp                              *
* Update date : 09/25/2023                              *
* Author      : pokai                            *
* Note        : ADIP hw1 code of image I/O			    *
*             : validation for learning                 *
*********************************************************/

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <queue>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace std;
using namespace cv;
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
void gray_down(vector<vector<int>> vec, char out[], int levle, int cho) {
	unsigned char* img_1d = new unsigned char[256 * 256];
	vector<vector<int>> vec2;
	vec2 = vec;
	for (int i = 0; i < vec.size(); i++) {
		for (int j = 0; j < vec.size(); j++) {
			double cn = (255.0 / ((pow(2, (8 - levle))) - 1));
			//cout <<levle<< cn << endl;
			//cout << vec[i][j] << " " << (vec[i][j] >> levle) << " " << round((vec[i][j] >> levle) * cn) << endl;
			//cout << (2 ^ levle) << endl;

			img_1d[i * vec.size() + j] = round((vec[i][j] >> levle) * cn);
			vec2[i][j] = round((vec[i][j] >> levle) * cn);
		}
	}
	string aa = "3_img_lena_";
	aa += out[19];
	if (cho == 2) {
		aa = "3_img_baboon_";
		aa += out[21];
	}
	aa += "_bits";
	compare(vec, vec2, aa);
	FILE* output_file22;
	output_file22 = fopen(out, "wb");
	fwrite(img_1d, 1, 256 * 256, output_file22);
	delete[] img_1d;
	fclose(output_file22);
}
int shortestPath(vector<vector<int>>grid, int control, char out[]) {

	int n = grid.size();
	if (grid[0][0] != 0 || grid[n - 1][n - 1] != 0)
		return -1;
	if (n - 1 == 0)
		return 1;

	vector<vector<int>> path(n, vector<int>(n, 0)); // Matrix to record the path
	grid[0][0] = 1;
	path[0][0] = 1;

	queue<pair<int, int>> q;
	q.push({ 0, 0 });

	int count = 0;
	int flag = 0;

	vector<int> r = { 0, 1, 0, -1 };
	vector<int> c = { 1, 0,-1, 0 };
	vector<int> rr = { 0, 1, -1, 0, 1, -1, 1, -1 };
	vector<int> cc = { 1, 0, 0, -1, 1, 1, -1, -1 };
	if (control == 2) {
		r = { 0, 1, -1, 0, 1, -1, 1, -1 };
		c = { 1, 0, 0, -1, 1, 1, -1, -1 };

	}

	while (!q.empty()) {
		int size = q.size();
		queue<pair<int, int>> temp;
		count++;

		while (size--) {
			int row = q.front().first;
			int col = q.front().second;
			q.pop();

			for (int i = 0; i < r.size(); i++) {
				int x = row + r[i], y = col + c[i];
				int nx = row + r[(i + 1) % r.size()];
				int ny = col + c[(i + 1) % r.size()];

				if (control == 3 && x >= 0 && y >= 0 && x < n && y < n && nx >= 0 && ny >= 0 && nx < n && ny < n && grid[x][y] == 1 && grid[nx][ny] == 1) {
					int sx = row + r[i] + r[(i + 1) % r.size()], sy = col + c[i] + c[(i + 1) % c.size()];
					if (grid[sx][sy] == 0) {
						x = sx; y = sy;
					}
				}


				if (x == n - 1 && y == n - 1) {
					flag = 1;
					count++;
					path[x][y] = path[row][col] + 1;
					
					break;
				}

				if (x >= 0 && y >= 0 && x < n && y < n) {
					if (grid[x][y] == 0) {
						temp.push({ x, y });
						grid[x][y] = 2;
						if (path[x][y] > 0) {
							path[x][y] = min(path[x][y], path[row][col] + 1);
						}
						else {
							path[x][y] = path[row][col] + 1;
						} // Record the path length
					}
				}
			}

			if (flag == 1)
				break;
		}

		q = temp;

		if (flag == 1)
			break;
	}

	if (flag == 0)
		return -1;

	// Backtrack to find the path
	int x = n - 1;
	int y = n - 1;
	
	
	vector<pair<int, int>> traveledPath;
	traveledPath.push_back({ x, y });
	//path[9][9] = count;
	while (x != 0 || y != 0) {
		for (int i = 0; i < 8; i++) {
			int nx = x - rr[i];
			int ny = y - cc[i];

			if (nx >= 0 && ny >= 0 && nx <= (n - 1) && ny <= (n - 1) && path[nx][ny] == (path[x][y] - 1)) {
				x = nx;
				y = ny;
				traveledPath.push_back({ x, y });
				break;
			}
		}
	}
	unsigned char* img_1d = new unsigned char[100];
	for (int i = 0; i < 100; i++) {
		img_1d[i] = 255;
	}
	
	// Print the traveled path
	cout << "Traveled path:" << endl;
	for (int i = traveledPath.size() - 1; i >= 0; i--) {
		cout << "(" << traveledPath[i].first << ", " << traveledPath[i].second << ")";
		if (i != 0) {
			cout << " -> ";
		}
	}
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			for (int n = traveledPath.size() - 1; n >= 0; n--) {
				if (i == traveledPath[n].first && j == traveledPath[n].second) {
					img_1d[i * 10 + j] = 0;
				}


			}
		}
	}


	cout << endl;

	FILE* output_file21;
	output_file21 = fopen(out, "wb");
	fwrite(img_1d, 1, 100, output_file21);
	delete[] img_1d;
	fclose(output_file21);

	return count;
}
void halftone(Mat& img, vector < vector<int>>& d2_vec) {
	for (int i = 0; i < d2_vec.size(); i += 16) {

		for (int j = 0; j < d2_vec.size(); j += 16) {
			double avg = 0;
			for (int k = 0; k < 16; k++) {
				for (int l = 0; l < 16; l++)
				{
					avg += d2_vec[i + k][j + l];
				}
			}

			avg /= 256.0;
			//cout << avg << endl;
			if (avg < 58) {
				avg = 8;
			}
			else if (avg > 201) {
				avg = 1;
			}
			else {
				avg -= 57;
				avg = -static_cast<int>(ceil(static_cast<double>(avg) / 23));
				avg += 8;

			}

			cv::Point center(j + 8, i + 8);
			cv::circle(img, center, avg, 0, -1);

		}
	}
	//cv::imshow("Circle Image", img);
	//cv::waitKey(0);
	for (int q = 0; q < d2_vec.size(); q++) {
		for (int r = 0; r < d2_vec.size(); r++) {
			d2_vec[q][r] = int(img.at<uchar>(q, r));
		}
	}
}

double cubicInterpolate(double p[4], double x) {
	return p[1] + 0.5f * x * (p[2] - p[0] + x * (2.0f * p[0] - 5.0f * p[1] + 4.0f * p[2] - p[3] + x * (3.0f * (p[1] - p[2]) + p[3] - p[0])));
}
double bicubicInterpolate(double p[4][4], double x, double y) {
	double arr[4];
	arr[0] = cubicInterpolate(p[0], y);
	arr[1] = cubicInterpolate(p[1], y);
	arr[2] = cubicInterpolate(p[2], y);
	arr[3] = cubicInterpolate(p[3], y);
	return cubicInterpolate(arr, x);
}
void neighbor(vector<vector<int>>& vec, double dim) {
	double len = vec[0].size();

	vector<vector<int>> tmp(dim, vector<int>(dim, 0));

	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			tmp[i][j] = vec[round((i / (dim - 1)) * (len - 1))][round((j / (dim - 1)) * (len - 1))];
		}
	}vec = tmp;

}
void bilinear(vector<vector<int>>& vec, double dim) {
	double len = vec[0].size();

	vector<vector<int>> tmp(dim, vector<int>(dim, 0));

	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			double srcx = ((i / (dim - 1)) * (len - 1)), srcy = ((j / (dim - 1)) * (len - 1));
			double x0 = floor(srcx); double x1 = x0 + 1;
			double y0 = floor(srcy); double y1 = y0 + 1;
			double dx = srcx - x0;
			double dy = srcy - y0;
			//cout << i << " " << j << endl;
			if (x0 == len - 1 || y0 == len - 1) {
				tmp[i][j] = vec[x0][y0];
			}
			else {
				tmp[i][j] = (1 - dx) * (1 - dy) * vec[x0][y0] +
					dx * (1 - dy) * vec[x1][y0] +
					(1 - dx) * dy * vec[x0][y1] +
					dx * dy * vec[x1][y1];
			}

		}
	}vec = tmp;

}
void bicubic(vector<vector<int>>& vec, double dim) {
	double len = vec[0].size();
	vector<vector<int>> tmp(dim, vector<int>(dim, 0));
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			// 计算原始图像坐标
			double srcx = ((i / (dim - 1)) * (len - 1)), srcy = ((j / (dim - 1)) * (len - 1));

			// 找到最接近的16个像素
			int x0 = floor(srcx);
			int y0 = floor(srcy);
			if (x0 == len - 1 || y0 == len - 1 || x0 == 0 || y0 == 0 || x0 == len - 2 || y0 == len - 2) {
				tmp[i][j] = vec[round(srcx)][round(srcy)];
			}


			else {
				// 计算插值所需的16个像素
				double values[4][4];
				for (int k = 0; k < 4; k++) {
					for (int l = 0; l < 4; l++) {
						int pixelX = x0 + k - 1;
						int pixelY = y0 + l - 1;
						values[k][l] = static_cast<double>(vec[pixelX][pixelY]);
					}
				}

				// 进行双三次插值
				double arr[4];
				arr[0] = cubicInterpolate(values[0], srcy - y0);
				arr[1] = cubicInterpolate(values[1], srcy - y0);
				arr[2] = cubicInterpolate(values[2], srcy - y0);
				arr[3] = cubicInterpolate(values[3], srcy - y0);
				double interpolatedValue = cubicInterpolate(arr, srcx - x0);

				tmp[i][j] = (interpolatedValue);
			}
		}
	}vec = tmp;
}
int main()
{
	//-----------------------1. Initial variable-----------------------//
	// Input  raw image name
	char  input_img[] = "image_file/lena256.raw";
	char  input_img2[] = "image_file/lena512.raw";
	char  input_img3[] = "image_file/lena128.raw";
	char  input_img4[] = "image_file/lena512_blur.raw";
	char  input_img5[] = "image_file/lena1024.raw";
	char input_map[] = "image_file/map10x10.raw";
	char input_baboon[] = "image_file/baboon256.raw";

	char output_img[] = "out_file/lena256_512_out.raw";

	char output_img2[] = "out_file/lena512_128_out.raw";

	char output_img3[] = "out_file/lena128_256_a_out.raw";
	char output_img4[] = "out_file/lena128_256_b_out.raw";
	char output_img5[] = "out_file/lena128_256_c_out.raw";

	char output_img6[] = "out_file/lena128_512_a_out.raw";
	char output_img7[] = "out_file/lena128_512_b_out.raw";
	char output_img8[] = "out_file/lena128_512_c_out.raw";

	char output_img9[] = "out_file/lena512_384_a_out.raw";
	char output_img10[] = "out_file/lena512_384_b_out.raw";
	char output_img11[] = "out_file/lena512_384_c_out.raw";

	char output_img12[] = "out_file/lena256_576_384_a_out.raw";
	char output_img13[] = "out_file/lena256_576_384_b_out.raw";
	char output_img14[] = "out_file/lena256_576_384_c_out.raw";

	char output_img15[] = "out_file/lena256_128_384_a_out.raw";
	char output_img16[] = "out_file/lena256_128_384_b_out.raw";
	char output_img17[] = "out_file/lena256_128_384_c_out.raw";

	char output_img18[] = "out_file/lena256_384_a_out.raw";
	char output_img19[] = "out_file/lena256_384_b_out.raw";
	char output_img20[] = "out_file/lena256_384_c_out.raw";

	char output_img21[] = "out_file/lena1024_cir_out.raw";
	FILE* input_file;
	FILE* input_file2;
	FILE* input_file3;
	FILE* input_file4;
	FILE* input_file5;
	FILE* input_gmap;
	FILE* input_bab;

	FILE* output_file;
	FILE* output_file2;

	FILE* output_file3;
	FILE* output_file4;
	FILE* output_file5;
	FILE* output_file6;
	FILE* output_file7;
	FILE* output_file8;

	FILE* output_file9;
	FILE* output_file10;
	FILE* output_file11;

	FILE* output_file12;
	FILE* output_file13;
	FILE* output_file14;
	FILE* output_file15;
	FILE* output_file16;
	FILE* output_file17;
	FILE* output_file18;
	FILE* output_file19;
	FILE* output_file20;
	FILE* output_file21;


	int width = 256;
	int height = 256;
	int size = width * height;
	int size2 = 512 * 512;
	int size3 = 384 * 384;

	unsigned char* img_lena_256 = new unsigned char[size]; // array for image data
	unsigned char* img_lena_512 = new unsigned char[size2];
	unsigned char* img_lena_512_blur = new unsigned char[size2];
	unsigned char* img_lena_128 = new unsigned char[128 * 128];
	unsigned char* img_lena_1024 = new unsigned char[1024 * 1024];
	unsigned char* grid_map = new unsigned char[100];
	unsigned char* img_bab = new unsigned char[size2];


	unsigned char* img_lena_256_a = new unsigned char[size2];

	unsigned char* img_1_1 = new unsigned char[size2];
	unsigned char* img_1_2 = new unsigned char[128 * 128];
	unsigned char* img_1_3_a = new unsigned char[size];
	unsigned char* img_1_3_b = new unsigned char[size];
	unsigned char* img_1_3_c = new unsigned char[size];

	unsigned char* img_1_3_d = new unsigned char[size2];
	unsigned char* img_1_3_e = new unsigned char[size2];
	unsigned char* img_1_3_f = new unsigned char[size2];

	unsigned char* img_1_4_a = new unsigned char[size2];
	unsigned char* img_1_4_b = new unsigned char[size2];
	unsigned char* img_1_4_c = new unsigned char[size2];

	unsigned char* img_1_5_a1 = new unsigned char[size3];
	unsigned char* img_1_5_b1 = new unsigned char[size3];
	unsigned char* img_1_5_c1 = new unsigned char[size3];
	unsigned char* img_1_5_a2 = new unsigned char[size3];
	unsigned char* img_1_5_b2 = new unsigned char[size3];
	unsigned char* img_1_5_c2 = new unsigned char[size3];
	unsigned char* img_1_5_a3 = new unsigned char[size3];
	unsigned char* img_1_5_b3 = new unsigned char[size3];
	unsigned char* img_1_5_c3 = new unsigned char[size3];

	vector<vector<int>> image_256, image_512, image_512_b, image_128, image_1_1, image_1_2, image_1_3_a, image_1_3_b, image_1_3_c, image_1_3_d, image_1_3_e, image_1_3_f, image_1_4_a, image_1_4_b, image_1_4_c;
	vector<vector<int>> image_1_5_a1, image_1_5_b1, image_1_5_c1, image_1_5_a2, image_1_5_b2, image_1_5_c2, image_1_5_a3, image_1_5_b3, image_1_5_c3, g_map, image_bab, image_1024;
	//vector<vector<int>> one(256, vector<int>(256, 0)), two(256, vector<int>(256, 0)), three(256, vector<int>(256, 0)), four(256, vector<int>(256, 0)), five(256, vector<int>(256, 0)), six(256, vector<int>(256, 0)), seven(256, vector<int>(256, 0)), eight(256, vector<int>(256, 0));
	//initialize with 0
	//-----------------------1. Initial variable-----------------------//


	//-----------------------2. Read File-----------------------//
	// using fopen as example, fstream works too
	input_file = fopen(input_img, "rb");
	input_file2 = fopen(input_img2, "rb");
	input_file3 = fopen(input_img3, "rb");
	input_file4 = fopen(input_img4, "rb");
	input_file5 = fopen(input_img5, "rb");
	input_gmap = fopen(input_map, "rb");
	input_bab = fopen(input_baboon, "rb");

	if (input_file == NULL) {
		cout << "Loading File Error!" << endl;
		system("PAUSE");
		exit(0);
	}

	fread(img_lena_256, 1, size, input_file);
	fread(img_lena_512, 1, size2, input_file2);
	fread(img_lena_512_blur, 1, size2, input_file4);
	fread(img_lena_128, 1, 128 * 128, input_file3);
	fread(img_lena_1024, 1, 1024 * 1024, input_file5);
	fread(grid_map, 1, 100, input_gmap);
	fread(img_bab, 1, size2, input_bab);
	//-----------------------2. Read File-----------------------//
	// 
	//-----------------------3. image processing-----------------------//
	for (int i = 0; i < 10; i++) {
		vector<int> vec;
		for (int j = 0; j < 10; j++) {
			vec.push_back(int(grid_map[i * 10 + j]));
		}g_map.push_back(vec);
	}
	for (int i = 0; i < 256; i++) {
		vector<int> vec, vec2;
		for (int j = 0; j < 256; j++) {
			vec.push_back(int(img_lena_256[i * 256 + j]));
			vec2.push_back(int(img_bab[i * 256 + j]));
		}image_256.push_back(vec);
		image_bab.push_back(vec2);
	}
	for (int i = 0; i < 512; i++) {
		vector<int> vec, vec2;
		for (int j = 0; j < 512; j++) {
			vec.push_back(int(img_lena_512[i * 512 + j]));
			vec2.push_back(int(img_lena_512_blur[i * 512 + j]));
		}image_512.push_back(vec);
		image_512_b.push_back(vec2);
	}
	for (int i = 0; i < 1024; i++) {
		vector<int> vec;
		for (int j = 0; j < 1024; j++) {
			vec.push_back(int(img_lena_1024[i * 1024 + j]));

		}image_1024.push_back(vec);

	}
	for (int i = 0; i < 128; i++) {
		vector<int> vec;
		for (int j = 0; j < 128; j++) {
			vec.push_back(int(img_lena_128[i * 128 + j]));
		}image_128.push_back(vec);
	}//轉成vector方便操作



	//-----------------------3. image processing-----------------------//
	for (int i = 0; i < 256; i++) {
		vector<int> vec;
		for (int j = 0; j < 256; j++) {
			vec.push_back(image_256[i][j]);
			vec.push_back(image_256[i][j]);
		}image_1_1.push_back(vec);
		image_1_1.push_back(vec);
	}
	string a = "1.1";
	compare(image_1_1, image_512, a);


	for (int i = 0; i < 128; i++) {
		vector<int> vec;
		for (int j = 0; j < 128; j++) {
			vec.push_back(image_512[i * 4][j * 4]);
		}image_1_2.push_back(vec);
	}
	image_1_3_a = image_128;
	image_1_3_d = image_128;
	double START, END;
	START = clock();
	neighbor(image_1_3_a, 256.0);
	END = clock();
	cout << endl << "1.3 256 Neighbor運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;
	START = clock();
	neighbor(image_1_3_d, 512.0);
	END = clock();
	cout << endl << "1.3 512 Neighbor運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;
	a = "1.3_256_neighbor";
	compare(image_1_3_a, image_256, a);
	a = "1.3_512_neighbor";
	compare(image_1_3_d, image_512, a);

	image_1_3_b = image_128;
	image_1_3_e = image_128;
	START = clock();
	bilinear(image_1_3_b, 256.0);
	END = clock();
	cout << endl << "1.3 256 Bilinear 運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;
	START = clock();
	bilinear(image_1_3_e, 512.0);
	END = clock();
	cout << endl << "1.3 512 Bilinear 運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;
	a = "1.3_256_bilinear";
	compare(image_1_3_b, image_256, a);
	a = "1.3_512_bilinear";
	compare(image_1_3_e, image_512, a);

	image_1_3_c = image_128;
	image_1_3_f = image_128;
	START = clock();
	bicubic(image_1_3_c, 256.0);
	END = clock();
	cout << endl << "1.3 256 Bicubic運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;
	START = clock();
	bicubic(image_1_3_f, 512.0);
	END = clock();
	cout << endl << "1.3 256 Bicubic運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;

	a = "1.3_256_bicubic";
	compare(image_1_3_c, image_256, a);
	a = "1.3_512_bicubic";
	compare(image_1_3_f, image_512, a);

	image_1_4_a = image_512_b;
	image_1_4_b = image_512_b;
	image_1_4_c = image_512_b;
	START = clock();
	neighbor(image_1_4_a, 384.0);
	END = clock();
	cout << endl << "1.4 Neighbor運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;
	START = clock();
	bilinear(image_1_4_b, 384.0);
	END = clock();
	cout << endl << "1.4 Bilinear運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;
	START = clock();
	bicubic(image_1_4_c, 384.0);
	END = clock();
	cout << endl << "1.4 Bicubic運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;

	image_1_5_a1 = image_256; image_1_5_a2 = image_256; image_1_5_a3 = image_256; image_1_5_b1 = image_256; image_1_5_b2 = image_256; image_1_5_b3 = image_256; image_1_5_c1 = image_256; image_1_5_c2 = image_256; image_1_5_c3 = image_256;
	START = clock();
	neighbor(image_1_5_a1, 576.0);
	neighbor(image_1_5_a1, 384.0);  
	END = clock();
	cout << endl << "1.5 上2.25下1.5 Neighbor運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;
	START = clock();
	bilinear(image_1_5_b1, 576.0);
	bilinear(image_1_5_b1, 384.0);
	END = clock();
	cout << endl << "1.5 上2.25下1.5 Bilinear運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;
	START = clock();
	bicubic(image_1_5_c1, 576.0);
	bicubic(image_1_5_c1, 384.0);
	END = clock();
	cout << endl << "1.5 上2.25下1.5 Bicubic運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;

	START = clock();
	neighbor(image_1_5_a2, 171.0);  
	neighbor(image_1_5_a2, 384.0);  
	END = clock();
	cout << endl << "1.5 下1.5上2.25 Neighbor運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;
	START = clock();
	bilinear(image_1_5_b2, 171.0);
	bilinear(image_1_5_b2, 384.0);
	END = clock();
	cout << endl << "1.5 下1.5上2.25 Bilinear運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;
	START = clock();
	bicubic(image_1_5_c2, 171.0);
	bicubic(image_1_5_c2, 384.0);
	END = clock();
	cout << endl << "1.5 下1.5上2.25 Bicubic運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;
	START = clock();
	neighbor(image_1_5_a3, 384.0);  
	END = clock();
	cout << endl << "1.5 上1.5 Neighbor運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;
	START = clock();
	bilinear(image_1_5_b3, 384.0);
	END = clock();
	cout << endl << "1.5 上1.5 Bilinear運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;
	START = clock();
	bicubic(image_1_5_c3, 384.0);
	END = clock();
	cout << endl << "1.5 上1.5 Bicubic運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;

	vector < vector<int>> map_72, map_72_145, map_72_145_218;
	map_72 = g_map;
	map_72_145 = g_map;
	map_72_145_218 = g_map;
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			if (map_72[i][j] == 72) {
				map_72[i][j] = 0;
				map_72_145[i][j] = 0;
				map_72_145_218[i][j] = 0;
			}
			else if (map_72_145[i][j] == 145) {
				map_72_145[i][j] = 0;
				map_72_145_218[i][j] = 0;
			}
			else if (map_72_145_218[i][j] == 218) {
				map_72_145_218[i][j] = 0;
				map_72_145[i][j] = 1;
			}
			else {
				map_72[i][j] = 1;
				map_72_145[i][j] = 1;
				map_72_145_218[i][j] = 1;
			}
		}
	}

	int cnt;
	char outmap[] = "out_file/map72_4.raw";
	cnt = shortestPath(map_72, 1, outmap);
	cout << "總步數 :" << cnt << endl;

	char outmap2[] = "out_file/map72_8.raw";
	cnt = shortestPath(map_72, 2, outmap2);
	cout << "總步數 :" << cnt << endl;

	char outmap3[] = "out_file/map72_m.raw";
	cnt = shortestPath(map_72, 3, outmap3);
	cout << "總步數 :" << cnt << endl;

	char outmap4[] = "out_file/map72145_4.raw";
	cnt = shortestPath(map_72_145, 1, outmap4);
	cout << "總步數 :" << cnt << endl;

	char outmap5[] = "out_file/map72145_8.raw";
	cnt = shortestPath(map_72_145, 2, outmap5);
	cout << "總步數 :" << cnt << endl;

	char outmap6[] = "out_file/map72145_m.raw";
	cnt = shortestPath(map_72_145, 3, outmap6);
	cout << "總步數 :" << cnt << endl;

	char outmap7[] = "out_file/map72145218_4.raw";
	cnt = shortestPath(map_72_145_218, 1, outmap7);
	cout << "總步數 :" << cnt << endl;

	char outmap8[] = "out_file/map72145218_8.raw";
	cnt = shortestPath(map_72_145_218, 2, outmap8);
	cout << "總步數 :" << cnt << endl;

	char outmap9[] = "out_file/map72145218_m.raw";
	cnt = shortestPath(map_72_145_218, 3, outmap9);
	cout << "總步數 :" << cnt << endl;

	for (int name = 0; name < 8; name++) {
		char out_gray[] = "out_file/lena_gray_0bits.raw";
		out_gray[19] = name + '0';

		gray_down(image_256, out_gray, name, 1);
		char out_gray2[] = "out_file/baboon_gray_0bits.raw";
		out_gray2[21] = name + '0';
		gray_down(image_bab, out_gray2, name, 2);

		//a = "lena_256_1";
		//compare(image_256, out_gray, a);
	}

	Mat img(Size(1024, 1024), CV_8UC1, cv::Scalar(255));
	halftone(img, image_1024);



	//-----------------------4. Save Image as raw format-----------------------//
	//image = four;
	//cout << four[128][255] << endl;
	for (int i = 0; i < 512; i++) {
		vector<int> vec;
		for (int j = 0; j < 512; j++) {
			img_lena_256_a[i * 512 + j] = image_1_1[i][j];
			img_1_3_d[i * 512 + j] = image_1_3_d[i][j];
			img_1_3_e[i * 512 + j] = image_1_3_e[i][j];
			img_1_3_f[i * 512 + j] = image_1_3_f[i][j];

		}
	}
	for (int i = 0; i < 384; i++) {
		vector<int> vec;
		for (int j = 0; j < 384; j++) {
			img_1_4_a[i * 384 + j] = image_1_4_a[i][j];
			img_1_4_b[i * 384 + j] = image_1_4_b[i][j];
			img_1_4_c[i * 384 + j] = image_1_4_c[i][j];

			img_1_5_a1[i * 384 + j] = image_1_5_a1[i][j];
			img_1_5_b1[i * 384 + j] = image_1_5_b1[i][j];
			img_1_5_c1[i * 384 + j] = image_1_5_c1[i][j];
			img_1_5_a2[i * 384 + j] = image_1_5_a2[i][j];
			img_1_5_b2[i * 384 + j] = image_1_5_b2[i][j];
			img_1_5_c2[i * 384 + j] = image_1_5_c2[i][j];
			img_1_5_a3[i * 384 + j] = image_1_5_a3[i][j];
			img_1_5_b3[i * 384 + j] = image_1_5_b3[i][j];
			img_1_5_c3[i * 384 + j] = image_1_5_c3[i][j];

		}
	}
	for (int i = 0; i < 128; i++) {
		vector<int> vec;
		for (int j = 0; j < 128; j++) {
			img_1_2[i * 128 + j] = image_1_2[i][j];
		}
	}
	for (int i = 0; i < 256; i++) {
		vector<int> vec;
		for (int j = 0; j < 256; j++) {
			img_1_3_a[i * 256 + j] = image_1_3_a[i][j];
			img_1_3_b[i * 256 + j] = image_1_3_b[i][j];
			img_1_3_c[i * 256 + j] = image_1_3_c[i][j];
		}
	}
	for (int i = 0; i < 1024; i++) {
		vector<int> vec;
		for (int j = 0; j < 1024; j++) {
			img_lena_1024[i * 1024 + j] = image_1024[i][j];
		}
	}

	output_file = fopen(output_img, "wb");
	output_file2 = fopen(output_img2, "wb");
	output_file3 = fopen(output_img3, "wb");
	output_file4 = fopen(output_img4, "wb");
	output_file5 = fopen(output_img5, "wb");
	output_file6 = fopen(output_img6, "wb");
	output_file7 = fopen(output_img7, "wb");
	output_file8 = fopen(output_img8, "wb");
	output_file9 = fopen(output_img9, "wb");
	output_file10 = fopen(output_img10, "wb");
	output_file11 = fopen(output_img11, "wb");

	output_file12 = fopen(output_img12, "wb");
	output_file13 = fopen(output_img13, "wb");
	output_file14 = fopen(output_img14, "wb");
	output_file15 = fopen(output_img15, "wb");
	output_file16 = fopen(output_img16, "wb");
	output_file17 = fopen(output_img17, "wb");
	output_file18 = fopen(output_img18, "wb");
	output_file19 = fopen(output_img19, "wb");
	output_file20 = fopen(output_img20, "wb");
	output_file21 = fopen(output_img21, "wb");

	fwrite(img_lena_256_a, 1, size2, output_file);
	fwrite(img_1_2, 1, 128 * 128, output_file2);
	fwrite(img_1_3_a, 1, size, output_file3);
	fwrite(img_1_3_b, 1, size, output_file4);
	fwrite(img_1_3_c, 1, size, output_file5);
	fwrite(img_1_3_d, 1, size2, output_file6);
	fwrite(img_1_3_e, 1, size2, output_file7);
	fwrite(img_1_3_f, 1, size2, output_file8);
	fwrite(img_1_4_a, 1, 384 * 384, output_file9);
	fwrite(img_1_4_b, 1, 384 * 384, output_file10);
	fwrite(img_1_4_c, 1, 384 * 384, output_file11);

	fwrite(img_1_5_a1, 1, 384 * 384, output_file12);
	fwrite(img_1_5_b1, 1, 384 * 384, output_file13);
	fwrite(img_1_5_c1, 1, 384 * 384, output_file14);
	fwrite(img_1_5_a2, 1, 384 * 384, output_file15);
	fwrite(img_1_5_b2, 1, 384 * 384, output_file16);
	fwrite(img_1_5_c2, 1, 384 * 384, output_file17);
	fwrite(img_1_5_a3, 1, 384 * 384, output_file18);
	fwrite(img_1_5_b3, 1, 384 * 384, output_file19);
	fwrite(img_1_5_c3, 1, 384 * 384, output_file20);
	fwrite(img_lena_1024, 1, 1024 * 1024, output_file21);


	//-----------------------5. Release memory-----------------------//
	delete[] img_lena_256;
	delete[] img_lena_256_a;
	delete[] img_lena_1024;
	delete[] img_1_2;
	delete[] img_1_3_a;
	delete[] img_1_3_b;
	delete[] img_1_3_c;
	delete[] img_1_3_d;
	delete[] img_1_3_e;
	delete[] img_1_3_f;
	delete[] img_1_4_a;
	delete[] img_1_4_b;
	delete[] img_1_4_c;
	delete[] img_1_5_a1;
	delete[] img_1_5_b1;
	delete[] img_1_5_c1;
	delete[] img_1_5_a2;
	delete[] img_1_5_b2;
	delete[] img_1_5_c2;
	delete[] img_1_5_a3;
	delete[] img_1_5_b3;
	delete[] img_1_5_c3;

	fclose(input_file);
	fclose(output_file);
	fclose(output_file2);
	fclose(output_file3);
	fclose(output_file4);
	fclose(output_file5);
	fclose(output_file6);
	fclose(output_file7);
	fclose(output_file8);
	fclose(output_file9);
	fclose(output_file10);
	fclose(output_file11);
	fclose(output_file12);
	fclose(output_file13);
	fclose(output_file14);
	fclose(output_file15);
	fclose(output_file16);
	fclose(output_file17);
	fclose(output_file18);
	fclose(output_file19);
	fclose(output_file20);
	fclose(output_file21);
}

//-----------------------5. Release memory-----------------------//

