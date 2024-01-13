/*#define _CRT_SECURE_NO_DEPRECATE

#include <iostream>
#include <vector>
#include <cmath>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace cv;
using namespace std;

#define PI 3.141592653589793

vector<vector<int>> read_raw(const char filename[], int sizew, int sizeh) {
    FILE* input_file;
    unsigned char* image = new unsigned char[sizew * sizeh];
    input_file = fopen(filename, "rb");
    fread(image, 1, sizew * sizeh, input_file);
    fclose(input_file);  // Don't forget to close the file

    vector<vector<int>> img;
    for (int i = 0; i < sizeh; i++) {
        vector<int> vec;
        for (int j = 0; j < sizew; j++) {
            vec.push_back(int(image[i * sizew + j]));
        }
        img.push_back(vec);
    }
    delete[] image;  // Free allocated memory
    return img;
}
void write_raw(vector<vector<int>> img, const char filename[], int sizew, int sizeh) {
    FILE* output_file;
    unsigned char* image = new unsigned char[sizew * sizeh];
    output_file = fopen(filename, "wb");
    for (int i = 0; i < sizeh; i++) {
        for (int j = 0; j < sizew; j++) {
            image[i * sizew + j] = static_cast<unsigned char>(img[i][j]);
        }
    }
    fwrite(image, 1, sizeh * sizew, output_file);
    fclose(output_file);
    delete[] image;  // Free allocated memory
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
clock_t a, b;
int main() {
    vector<vector<int>> baboon, rec_baboon,lena,rec_lena;
    baboon = read_raw("src/baboon_256.raw", 256, 256);
    lena = read_raw("src/lena_256.raw", 256, 256);

    a = clock();
    Mat mybaboon(256, 256, CV_32S);
    Mat mylena(256, 256, CV_32S);

    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            mybaboon.at<int>(i, j) = baboon[i][j];
            mylena.at<int>(i, j) = lena[i][j];
        }
    }

    int m = getOptimalDFTSize(mybaboon.rows);
    int n = getOptimalDFTSize(mybaboon.cols);

    Mat padded;
    copyMakeBorder(mybaboon, padded, 0, m - mybaboon.rows, 0, n - mybaboon.cols, BORDER_CONSTANT, Scalar::all(0));
    Mat planes[] = { Mat_<float>(padded), Mat_<float>::zeros(padded.size()) };  // Use Mat_<float> for zeros
    Mat complexImage;
    merge(planes, 2, complexImage);
    dft(complexImage, complexImage);
    fftShift(complexImage);
    split(complexImage, planes);
    magnitude(planes[0], planes[1], planes[0]);
    Mat magnitudeImage = planes[0];
    magnitudeImage += Scalar::all(1);
    log(magnitudeImage, magnitudeImage);
    normalize(magnitudeImage, magnitudeImage, 0, 1, NORM_MINMAX); 
    b = clock();

    cout << "DFT Caluclate time :" << double(b - a) / CLOCKS_PER_SEC << endl;

    vector<vector<int>> c_baboon(256, vector<int>(256, 0));
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            c_baboon[i][j] = 255 * (magnitudeImage.at<float>(i, j));
        }
    }

    write_raw(c_baboon, "image_file/cv_dft_baboon.raw", 256, 256);
    a = clock();
    Mat ifft(256,256, CV_32F,1);
    idft(complexImage, ifft, DFT_SCALE | DFT_REAL_OUTPUT);

    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            if (ifft.at<float>(i, j) < 0) {
                ifft.at<float>(i, j) =abs(ifft.at<float>(i, j));
            }
        }
    }
    //normalize(ifft, ifft, 0, 255, NORM_MINMAX);
    //fft *= 255;
    
    //ifft.convertTo(ifft, CV_8U);
    vector<vector<int>> c_idft_baboon(256, vector<int>(256, 0));
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
           // cout << ifft.at<float>(i, j) << endl;
            c_idft_baboon[i][j] =round((ifft.at<float>(i, j)));
        }
    }
    b = clock();
    cout << "IDFT Caluclate time :" << double(b - a) / CLOCKS_PER_SEC << endl;
    write_raw(c_idft_baboon, "image_file/cv_idft_baboon.raw", 256, 256);
    ifft.convertTo(ifft, CV_8U);
    imshow("IDFT", ifft);

    m = getOptimalDFTSize(mylena.rows);
    n = getOptimalDFTSize(mylena.cols);
    a = clock();
    Mat padded2;
    copyMakeBorder(mylena, padded2, 0, m - mylena.rows, 0, n - mylena.cols, BORDER_CONSTANT, Scalar::all(0));
    Mat planes2[] = { Mat_<float>(padded2), Mat_<float>::zeros(padded2.size()) };  // Use Mat_<float> for zeros
    Mat complexImage2;
    merge(planes2, 2, complexImage2);
    dft(complexImage2, complexImage2);
    fftShift(complexImage2);
    split(complexImage2, planes2);
    magnitude(planes2[0], planes2[1], planes2[0]);
    Mat magnitudeImage2 = planes2[0];
    magnitudeImage2 += Scalar::all(1);
    log(magnitudeImage2, magnitudeImage2);
    normalize(magnitudeImage2, magnitudeImage2, 0, 1, NORM_MINMAX);
    b = clock();
    cout << "DFT Caluclate time :" << double(b - a) / CLOCKS_PER_SEC << endl;

    vector<vector<int>> c_lena(256, vector<int>(256, 0));
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            c_lena[i][j] = 255 * (magnitudeImage2.at<float>(i, j));
        }
    }

    write_raw(c_lena, "image_file/cv_dft_lena.raw", 256, 256);
    a = clock();
    Mat ifft2(256, 256, CV_32F, 1);
    idft(complexImage2, ifft2, DFT_SCALE | DFT_REAL_OUTPUT);

    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            if (ifft2.at<float>(i, j) < 0) {
                ifft2.at<float>(i, j) = abs(ifft2.at<float>(i, j));
            }
        }
    }
    //normalize(ifft, ifft, 0, 255, NORM_MINMAX);
    //fft *= 255;

    //ifft.convertTo(ifft, CV_8U);
    vector<vector<int>> c_idft_lena(256, vector<int>(256, 0));
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            // cout << ifft.at<float>(i, j) << endl;
            c_idft_lena[i][j] = round((ifft2.at<float>(i, j)));
        }
    }
    b = clock();
    cout << "IDFT Caluclate time :" << double(b - a) / CLOCKS_PER_SEC << endl;
    write_raw(c_idft_lena, "image_file/cv_idft_lena.raw", 256, 256);
    //ifft2.convertTo(ifft2, CV_8U);
    string s = "1.4 lena : ";
    compare(lena, c_idft_lena, s);
    s = "1.4 baboon : ";
    compare(baboon, c_idft_baboon, s);
    

}*/
