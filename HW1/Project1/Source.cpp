/********************************************************
* Filename    : hw1.cpp                              *
* Update date : 09/19/2023                              *
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
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace cv;
using namespace std;
void rotate90(vector<vector<int>>& a, int cnt) {
	{
		vector<vector<int>> tmp(256, vector<int>(256, 0));
		int x, y, nx, ny;
		double vcos, vsin;
		double theta = 90;
		theta *= 3.14159265358979 / 180; // 轉弳度
		vsin = sin(theta * cnt), vcos = cos(theta * cnt);

		for (int i = 0; i < 256; i++) {
			for (int j = 0; j < 256; j++) {
				//2階旋轉矩陣
				int new_x =round( i * vcos - j * vsin); 
				int new_y =round( i * vsin + j * vcos);

				// 依據角度選擇平移量
				if (cnt == 1) {
					tmp[new_x + 255][new_y] = a[i][j];
				}
				if (cnt == 2) {
					tmp[new_x + 255][new_y + 255] = a[i][j];
				}
				if (cnt == 3) {
					tmp[new_x][new_y + 255] = a[i][j];
				}


			}
		}

		// Copy the rotated matrix back to the original matrix
		a = tmp;
	}
}

int main()
{
	//-----------------------1. Initial variable-----------------------//
	// Input  raw image name
	char  input_img[] = "lena256.raw";                 
	char output_img[] = "lena256out.raw";
	char output_img2[] = "lena256out_divide.raw";
	char output_img3[] = "lena256out_minpool.raw";
	char output_img4[] = "lena256out_bright.raw";
	char output_img5[] = "lena256out_bright_random.raw";

	FILE* input_file;
	FILE* output_file;
	FILE* output_file2;
	FILE* output_file3;
	FILE* output_file4;
	FILE* output_file5;

	int width = 256;
	int height = 256;
	int size = width * height;

	unsigned char* img_lena = new unsigned char[size]; // array for image data
	unsigned char* img_reverse = new unsigned char[size];
	unsigned char* img_min = new unsigned char[size];
	unsigned char* img_bright = new unsigned char[size];
	unsigned char* img_bright_r = new unsigned char[size];

	vector<vector<int>> image, image_min, image_bright, image_bright_r;
	vector<vector<int>> one(256, vector<int>(256, 0)), two(256, vector<int>(256,0)), three(256, vector<int>(256,0)), four(256, vector<int>(256,0)), five(256, vector<int>(256, 0)), six(256, vector<int>(256,0)), seven(256, vector<int>(256, 0)), eight(256, vector<int>(256, 0));
	//initialize with 0
	//-----------------------1. Initial variable-----------------------//


	//-----------------------2. Read File-----------------------//
	// using fopen as example, fstream works too
	input_file = fopen(input_img, "rb");

	if (input_file == NULL) {
		cout << "Loading File Error!" << endl;
		system("PAUSE");
		exit(0);
	}

	fread(img_lena, 1, size, input_file);
	//-----------------------2. Read File-----------------------//
	// 
	//-----------------------3. image processing-----------------------//
	for (int i = 0; i < 256; i++) {
		vector<int> vec;
		for (int j = 0; j < 256; j++) {
			vec.push_back(int(img_lena[i * 256 + j]));
		}image.push_back(vec);
	}//轉成vector方便操作


	cout << "1-2_a.1 : " << image[42][145] << endl;
	cout << "1-2_a.2 : " << int(img_lena[42018]) << endl;

	//bright
	for (int i = 0; i < 256; i++) {
		vector<int> vec;
		for (int j = 0; j < 256; j++) {
			vec.push_back(int(img_lena[i * 256 + j]) + 70);
		}image_bright.push_back(vec);
	}
	int n_max = -1000, n_min = 1000;
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			if (image_bright[i][j] < n_min) {
				n_min = image_bright[i][j];
			}
			if (image_bright[i][j] > n_max) {
				n_max = image_bright[i][j];
			}
		}
	}
	cout << "max : "<<n_max << "min : " << n_min << endl;
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			image_bright[i][j] = int(255.0 * (double((image_bright[i][j] - n_min)) / (n_max - n_min)));//rescale
		}
	}
	//bright_random
	srand(time(NULL));
	/* 指定亂數範圍 */
	double m = 60.0;
	double mi = -60.0;
	/* 產生 (min , max) 的浮點數亂數 */
	for (int i = 0; i < 256; i++) {
		vector<int> vec;
		for (int j = 0; j < 256; j++) {
			double x = (m - mi) * rand() / (RAND_MAX + 1.0) + mi;
			vec.push_back(int(img_lena[i * 256 + j]) + x);
		}image_bright_r.push_back(vec);
	}
	int r_max = -1000, r_min = 1000;
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			if (image_bright_r[i][j] < r_min) {
				r_min = image_bright_r[i][j];
			}
			if (image_bright_r[i][j] > r_max) {
				r_max = image_bright_r[i][j];
			}
		}
	}
	
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			image_bright_r[i][j] = int(255.0 * (double((image_bright_r[i][j] - r_min)) / (r_max - r_min)));//rescale
		}
	}


	//minpooling
	for (int i = 0; i < 32; i++) {
		vector<int> vec;
		for (int j = 0; j < 32; j++) {
			int min = 1000;
			for (int k = 0; k < 8; k++) {
				for (int l = 0; l < 8; l++) {
					if (image[i * 8 + k][j * 8 + l] < min) {
						min = image[i * 8 + k][j * 8 + l];
					}
				}
			}
			for (int t = 0; t < 8; t++) {
				vec.push_back(min);
			}
		}
		for (int t = 0; t < 8; t++) {
			image_min.push_back(vec);
		}

	}



	double x1 = 127, x2 = 0, x3 = 0, x4 = 0, x5 = 127, x6 = 127;
	double y1 = 0, y2 = 127, y3 = 127, y4 = 128, y5 = 0, y6 = 255;
	//象限分割
	
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			if (i <= x1) {
				if (j <= y2) {
					if (j <= y3) {
						//case 1
						one[i][j] = image[i][j];
					}
					else {
						//case 2
						two[i][j] = image[i][j];
					}
				}

				else if (j>y2){
					if (j < y4) {
						//case 3
						three[i][j] = image[i][j];
					}
					else {
						//case 4
						four[i][j] = image[i][j];

					}
				}
			}
			else {
				if (j <= y2) {
					if (j <= y5) {
						//case 5
						five[i][j] = image[i][j];
					}
					else {
						//case 6
						six[i][j] = image[i][j];
					}
				}

				else {
					if (j < y6) {
						//case 7
						seven[i][j] = image[i][j];
					}
					else {
						//case 8
						eight[i][j] = image[i][j];
					}
				}

			}
			
		}y1 += 1; x2 += 1;
		if (i < 128) {
			x3 += 1; y3 -= 1;
			x4 += 1; y4 += 1;
		}
		if (i > 127) {
			x5 += 1; y5 += 1;
			x6 += 1; y6 -= 1;
		}
	}


	x1 = 127, x2 = 0, x3 = 0, x4 = 0, x5 = 127, x6 = 127;
	 y1 = 0, y2 = 127, y3 = 127, y4 = 128, y5 = 0, y6 = 255;

	rotate90(one, 3); rotate90(three, 2); rotate90(four, 3); rotate90(six, 2); rotate90(eight, 2);
	
	// 
	//rebuild image
	for (int i = 0; i < 256; i++) {
		
		for (int j = 0; j < 256; j++) {
			if (i <= x1) {
				if (j <= y2) {
					if (j <= y3) {
						//case 1
						image[i][j] = eight[i][j];
					}
					else {
						//case 2
						image[i][j] = two[i][j];
					}
				}

				else {
					if (j < y4) {
						//case 3
						image[i][j] = six[i][j];
					}
					else {
						//case 4
						image[i][j] = one[i][j];

					}
				}
			}
			else {
				if (j <= y2) {
					if (j <= y5) {
						//case 5
						image[i][j] = five[i][j];
					}
					else {
						//case 6
						image[i][j] = three[i][j];
					}
				}

				else {
					if (j < y6) {
						//case 7
						image[i][j] = seven[i][j];
					}
					else {
						//case 8
						image[i][j] = four[i][j];
					}
				}

			}
		}
		y1 += 1; x2 += 1;
		if (i < 128) {
			x3 += 1; y3 -= 1;
			x4 += 1; y4 += 1;
		}
		if (i > 127) {
			x5 += 1; y5 += 1;
			x6 += 1; y6 -= 1;
		}
	}
	//-----------------------3. image processing-----------------------//



	//-----------------------4. Save Image as raw format-----------------------//
	//image = four;
	//cout << four[128][255] << endl;
	for (int i = 0; i < 256; i++) {
		vector<int> vec;
		for (int j = 0; j < 256; j++) {
			img_reverse[i * 256 + j] = image[i][j];
			img_min[i * 256 + j] = image_min[i][j];
			img_bright[i * 256 + j] = image_bright[i][j];
			img_bright_r[i * 256 + j] = image_bright_r[i][j];
		}
	}

	output_file = fopen(output_img, "wb");
	output_file2 = fopen(output_img2, "wb");
	output_file3 = fopen(output_img3, "wb");
	output_file4 = fopen(output_img4, "wb");
	output_file5 = fopen(output_img5, "wb");

	fwrite(img_lena, 1, size, output_file);//org
	fwrite(img_reverse, 1, size, output_file2);//reverse
	fwrite(img_min, 1, size, output_file3);//minpooling
	fwrite(img_bright, 1, size, output_file4);//birght
	fwrite(img_bright_r, 1, size, output_file5);//random bright

	//-----------------------5. Release memory-----------------------//
	delete[] img_lena;
	delete[] img_reverse;
	delete[] img_min;
	delete[] img_bright;
	delete[] img_bright_r;
	fclose(input_file);
	fclose(output_file);
	fclose(output_file2);
	fclose(output_file3);
	fclose(output_file4);
	fclose(output_file5);
	cv::Mat img(512, 512, CV_8UC1, cv::Scalar(255));
	//-----------------------5. Release memory-----------------------//
	
	//-----------------------6. opencv draw doreamon-----------------------// 
	

	//cv::circle(InputOutputArray img, Point center, int radius, const Scalar & color, int thickness = 1, int lineType = LINE_8, int shift = 0)
	//cv::line(InputOutputArray img, Point pt1, Point pt2, const Scalar &color, int thickness=1, int lineType=LINE_8, int shift=0)
	//cv::ellipse(Mat& img, Point center, Size axes, double angle, double startAngle, double endAngle, const Scalar& color, int thickness=1, int lineType=8, int shift=0)
	// Head
	cv::circle(img, cv::Point(256, 256), 200, cv::Scalar(100), -1, 16);

	
	// Eye lines
	cv::ellipse(img, cv::Point(256, 330), cv::Size(185, 120), 180, 180, 0, cv::Scalar(0), 4, 16);
	cv::ellipse(img, cv::Point(256, 330), cv::Size(185, 120), 180, 180, 0, cv::Scalar(255), -1, 16);
	cv::ellipse(img, cv::Point(256, 330), cv::Size(285, 150), 0, 180, 0, cv::Scalar(255), -1, 16);
	cv::circle(img, cv::Point(256, 256), 200, cv::Scalar(0), 4, 16);
	//eye
	cv::ellipse(img, cv::Point(200, 200), cv::Size(54, 70), 0, 0, 360, cv::Scalar(0), 4, 32);
	cv::ellipse(img, cv::Point(200, 200), cv::Size(54, 70), 0, 0, 360, cv::Scalar(255), -1);
	cv::ellipse(img, cv::Point(312, 200), cv::Size(54, 70), 0, 0, 360, cv::Scalar(0), 4, 32);
	cv::ellipse(img, cv::Point(312, 200), cv::Size(54, 70), 0, 0, 360, cv::Scalar(255), -1);
	cv::ellipse(img, cv::Point(220, 200), cv::Size(15, 20), 0, 0, 360, cv::Scalar(0), -1);
	cv::ellipse(img, cv::Point(292, 200), cv::Size(15, 20), 0, 0, 360, cv::Scalar(0), -1);
	cv::ellipse(img, cv::Point(220, 200), cv::Size(5, 10), 0, 0, 360, cv::Scalar(255), -1);
	cv::ellipse(img, cv::Point(292, 200), cv::Size(5, 10), 0, 0, 360, cv::Scalar(255), -1);
	//nose
	cv::circle(img, cv::Point(256, 265), 20, cv::Scalar(0), 4, 16);
	cv::circle(img, cv::Point(256, 265), 20, cv::Scalar(150), -1, 16);
	cv::circle(img, cv::Point(250, 265), 5, cv::Scalar(255), -1, 16);
	//mouth
	cv::ellipse(img, cv::Point(256, 300), cv::Size(120, 100), 0, 0, 180, cv::Scalar(0), 4, 32);
	cv::line(img, cv::Point(256, 287), cv::Point(256, 400), Scalar(0), 2, 32);
	//thread
	cv::line(img, cv::Point(180, 340), cv::Point(50, 370), Scalar(0), 2, 32);
	cv::line(img, cv::Point(180, 310), cv::Point(50, 310), Scalar(0), 2, 32);
	cv::line(img, cv::Point(180, 280), cv::Point(50, 250), Scalar(0), 2, 32);

	cv::line(img, cv::Point(450, 370), cv::Point(320, 340), Scalar(0), 2, 32);
	cv::line(img, cv::Point(450, 310), cv::Point(320, 310), Scalar(0), 2, 32);
	cv::line(img, cv::Point(450, 250), cv::Point(320, 280), Scalar(0), 2, 32);
	//text
	string text = "111C71007";
	int font_face = cv::FONT_HERSHEY_COMPLEX;
	cv::putText(img, text, cv::Point(160, 490), font_face, 1, cv::Scalar(180), 2, 1, 0);
	cv::imshow("hw1", img);
	cv::waitKey(5000);
	imwrite("111c71007.png", img);
	cv::destroyAllWindows();
}
//-----------------------6. opencv draw doreamon-----------------------//