/*#define _CRT_SECURE_NO_DEPRECATE
#include <iostream>
#include <vector>
#include <cmath>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace cv;
using namespace std;

const double PI = 3.14159265358979323846;

cv::Mat visualizeSpectrum(const cv::Mat& complexImage) {
    cv::Mat planes[2];
    cv::split(complexImage, planes);
    cv::Mat magnitudeImage(256, 256, CV_32S);
    cv::magnitude(planes[0], planes[1], magnitudeImage);

    // Add 1 to avoid log(0)
    magnitudeImage += Scalar::all(1);
    cv::log(magnitudeImage, magnitudeImage);

    // Normalize the magnitude spectrum for display
    cv::normalize(magnitudeImage, magnitudeImage, 0, 1, NORM_MINMAX);

    return magnitudeImage;
}
// Function to create a Gaussian filter in the frequency domain
void fftShift(Mat& image) {
    int cx = image.cols / 2;
    int cy = image.rows / 2;

    Mat q0(image, Rect(0, 0, cx, cy));   // Top-Left
    Mat q1(image, Rect(cx, 0, cx, cy));  // Top-Right
    Mat q2(image, Rect(0, cy, cx, cy));  // Bottom-Left
    Mat q3(image, Rect(cx, cy, cx, cy)); // Bottom-Right

    Mat tmp;
    q0.copyTo(tmp);
    q3.copyTo(q0);
    tmp.copyTo(q3);

    q1.copyTo(tmp);
    q2.copyTo(q1);
    tmp.copyTo(q2);
}

vector<vector<double>> read_raw(const char filename[], int sizew, int sizeh) {
    FILE* input_file;
    unsigned char* image = new unsigned char[sizew * sizeh];
    input_file = fopen(filename, "rb");
    fread(image, 1, sizew * sizeh, input_file);
    vector<vector<double>> img;
    for (int i = 0; i < sizeh; i++) {
        vector<double> vec;
        for (int j = 0; j < sizew; j++) {
            vec.push_back(static_cast<double>(image[i * sizew + j]));  // Convert to double
        }
        img.push_back(vec);
    }
    return img;
}

void write_raw(vector<vector<int>> img, const char filename[], int sizew, int sizeh) {
    FILE* output_file;
    unsigned char* image = new unsigned char[sizew * sizeh];
    output_file = fopen(filename, "wb");
    for (int i = 0; i < sizeh; i++) {
        for (int j = 0; j < sizew; j++) {
            image[i * sizew + j] = static_cast<unsigned char>(img[i][j]);  // Convert to unsigned char
        }
    }
    fwrite(image, 1, sizeh * sizew, output_file);
    fclose(output_file);
}

cv::Mat createButterworthFilter(int rows, int cols, double D0, int n, bool isLowPass) {
    cv::Mat filter(rows, cols, CV_64F);
    int centerX = cols / 2;
    int centerY = rows / 2;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double distance = sqrt((i - centerY) * (i - centerY) + (j - centerX) * (j - centerX));

            double value;
            if (isLowPass) {
                value = 1 / (1 + pow(distance / D0, 2 * n));
            }
            else {
                value = 1 - 1 / (1 + pow(distance / D0, 2 * n));
            }

            filter.at<double>(i, j) = value;
        }
    }

    return filter;
}

int main() {
    // Read the input image
    vector<vector<double>> mbaboon;
    mbaboon = read_raw("src/baboon_256.raw", 256, 256);

    cv::Mat baboon(256, 256, CV_64F);

    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            baboon.at<double>(i, j) = mbaboon[i][j];
        }
    }

    // Check if the image is loaded successfully
    if (baboon.empty()) {
        std::cerr << "Error: Unable to load the image." << std::endl;
        return -1;
    }

    // Apply DFT to the image

    int D0_values[] = { 10, 20, 40, 100 };

    for (int shi = 0; shi < 2; shi++) {
        for (int D0 : D0_values) {
            for (int n : {1, 2, 3, 4}) {
                for (int i = 0; i < 256; i++) {
                    for (int j = 0; j < 256; j++) {
                        baboon.at<double>(i, j) = mbaboon[i][j];
                    }
                }
                dft(baboon, baboon);

                if (shi == 0) {
                    fftShift(baboon);
                }

                int rows = baboon.rows;
                int cols = baboon.cols;
                // Create Gaussian low-pass filter
                cv::Mat lowPassFilter = createButterworthFilter(rows, cols, D0, n, true);


                // Create Gaussian high-pass filter
                cv::Mat highPassFilter = createButterworthFilter(rows, cols, D0, n, false);


                // Multiply the filters with the shifted DFT image
                cv::Mat lowPassResult = baboon.mul(lowPassFilter);
                cv::Mat highPassResult = baboon.mul(highPassFilter);


                if (shi == 0) {
                    fftShift(lowPassResult);
                    fftShift(highPassResult);
                }


                // Apply IDFT to get the filtered images
                cv::idft(lowPassResult, lowPassResult);
                cv::idft(highPassResult, highPassResult);

                // Normalize the results to the range [0, 255] for display
                normalize(lowPassResult, lowPassResult, 0, 255, NORM_MINMAX);
                normalize(highPassResult, highPassResult, 0, 255, NORM_MINMAX);

                lowPassFilter *= 255;  // You can adjust this scaling factor
                highPassFilter *= 255; // You can adjust this scaling factor
                if (shi == 0) {
                    imwrite("image_file/butterworth/ORIGIN_LPF_D0_" + std::to_string(D0) + "_n_" + std::to_string(n)+"_SPRCTRUM.png", lowPassFilter);
                    imwrite("image_file/butterworth/ORIGIN_HPF_D0_" + std::to_string(D0) + "_n_" + std::to_string(n)+ "_SPRCTRUM.png", highPassFilter);
                }
                else {
                    imwrite("image_file/butterworth/LPF_D0_" + std::to_string(D0) + "_n_" + std::to_string(n)+"_SPRCTRUM.png", lowPassFilter);
                    imwrite("image_file/butterworth/HPF_D0_" + std::to_string(D0) + "_n_" + std::to_string(n)+ "_SPRCTRUM.png", highPassFilter);
                }
                // Convert the results to 8-bit unsigned integer type
                lowPassResult.convertTo(lowPassResult, CV_8U);
                highPassResult.convertTo(highPassResult, CV_8U);



                // Display the filtered images
                if (shi == 0) {
                    imwrite("image_file/butterworth/ORIGIN_Butterworth_LPF_D0_" + std::to_string(D0) + "_n_" + std::to_string(n) + ".png", lowPassResult);
                    imwrite("image_file/butterworth/ORIGIN_Butterworth_HPF_D0_" + std::to_string(D0) + "_n_" + std::to_string(n) + ".png", highPassResult);
                }
                else {
                    imwrite("image_file/butterworth/Butterworth_LPF_D0_" + std::to_string(D0) + "_n_" + std::to_string(n) + ".png", lowPassResult);
                    imwrite("image_file/butterworth/Butterworth_HPF_D0_" + std::to_string(D0) + "_n_" + std::to_string(n) + ".png", highPassResult);
                }

                cv::waitKey(0);
            }
        }
    }

    return 0;
}

*/