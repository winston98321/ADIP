//
//#define _CRT_SECURE_NO_DEPRECATE
//#include <iostream>
//#include <vector>
//#include <cmath>
//#include <opencv2/core/core.hpp>
//#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/imgproc/imgproc.hpp>
//#include<complex>
//
//using namespace cv;
//using namespace std;
//
//const double PI = 3.14159265358979323846;
//void comparee(vector<vector<double>> vec, vector<vector<double>> src, string a) {
//    double MSE = 0;
//    int len = vec[0].size();
//    for (int i = 0; i < len; i++) {
//        for (int j = 0; j < len; j++) {
//            MSE += pow((src[i][j] - vec[i][j]), 2);
//        }
//    }
//    MSE = (MSE / (len * len));
//    cout << a << " MSE :" << MSE << endl;
//    cout << a << " PSNR :" << 10 * log10((255 * 255) / MSE) << endl;
//}
//cv::Mat visualizeSpectrum(const cv::Mat& complexImage) {
//    cv::Mat planes[2];
//    cv::split(complexImage, planes);
//    cv::Mat magnitudeImage(256, 256, CV_32S);
//    cv::magnitude(planes[0], planes[1], magnitudeImage);
//
//    magnitudeImage += Scalar::all(1);
//    cv::log(magnitudeImage, magnitudeImage);
//
//     Normalize the magnitude spectrum for display
//    cv::normalize(magnitudeImage, magnitudeImage, 0, 1, NORM_MINMAX);
//
//    return magnitudeImage;
//}
// Function to create a Gaussian filter in the frequency domain
//void fftShift(Mat& image) {
//    int cx = image.cols / 2;
//    int cy = image.rows / 2;
//
//    Mat q0(image, Rect(0, 0, cx, cy));   // Top-Left
//    Mat q1(image, Rect(cx, 0, cx, cy));  // Top-Right
//    Mat q2(image, Rect(0, cy, cx, cy));  // Bottom-Left
//    Mat q3(image, Rect(cx, cy, cx, cy)); // Bottom-Right
//
//    Mat tmp;
//    q0.copyTo(tmp);
//    q3.copyTo(q0);
//    tmp.copyTo(q3);
//
//    q1.copyTo(tmp);
//    q2.copyTo(q1);
//    tmp.copyTo(q2);
//}
//
//vector<vector<double>> read_raw(const char filename[], int sizew, int sizeh) {
//    FILE* input_file;
//    unsigned char* image = new unsigned char[sizew * sizeh];
//    input_file = fopen(filename, "rb");
//    fread(image, 1, sizew * sizeh, input_file);
//    vector<vector<double>> img;
//    for (int i = 0; i < sizeh; i++) {
//        vector<double> vec;
//        for (int j = 0; j < sizew; j++) {
//            vec.push_back(static_cast<double>(image[i * sizew + j]));  // Convert to double
//        }
//        img.push_back(vec);
//    }
//    return img;
//}
//
//void write_raw(vector<vector<int>> img, const char filename[], int sizew, int sizeh) {
//    FILE* output_file;
//    unsigned char* image = new unsigned char[sizew * sizeh];
//    output_file = fopen(filename, "wb");
//    for (int i = 0; i < sizeh; i++) {
//        for (int j = 0; j < sizew; j++) {
//            image[i * sizew + j] = static_cast<unsigned char>(img[i][j]);  // Convert to unsigned char
//        }
//    }
//    fwrite(image, 1, sizeh * sizew, output_file);
//    fclose(output_file);
//}
//
//cv::Mat createLaplacian(int rows, int cols) {
//    Mat filter = Mat::zeros(Size(rows,cols), CV_64F);
//    int centerX = cols / 2;
//    int centerY = rows / 2;
//    filter.at<double>(centerX, centerY) = 4;
//    filter.at<double>(centerX+1, centerY) = -1;
//    filter.at<double>(centerX, centerY+1) = -1;
//    filter.at<double>(centerX-1, centerY) = -1;
//    filter.at<double>(centerX, centerY-1) = -1;
//   
//    return filter;
//}
//cv::Mat createButterworthFilter(int rows, int cols, double D0, int n, bool isLowPass) {
//    cv::Mat filter(rows, cols, CV_64F);
//    int centerX = cols / 2;
//    int centerY = rows / 2;
//
//    for (int i = 0; i < rows; ++i) {
//        for (int j = 0; j < cols; ++j) {
//            double distance = sqrt((i - centerY) * (i - centerY) + (j - centerX) * (j - centerX));
//
//            double value;
//            if (isLowPass) {
//                value = 1 / (1 + pow(distance / D0, 2 * n));
//            }
//            else {
//                value = 1 - 1 / (1 + pow(distance / D0, 2 * n));
//            }
//
//            filter.at<double>(i, j) = value;
//        }
//    }
//    Mat toMerge[] = { filter, filter };
//    merge(toMerge, 2, filter);
//    return filter;
//}
//cv::Mat createShiftDegradation(int rows, int cols, double a, double b) {
//    cv::Mat degradationFunction(rows, cols, CV_64FC2);
//    cv::Mat out(rows, cols, CV_64F);
//    for (int y = 0; y < rows; y++) {
//        for (int x = 0; x < cols; x++) {
//            double phaseShift = 2 * PI * (a * x / cols + b * y / rows);
//            out.at<double>(y,x) =pow(pow(cos(phaseShift),2)+pow(sin(phaseShift),2),0.5);
//            degradationFunction.at<std::complex<double>>(y, x) = complexValue;
//        }
//    }
//
//    return out;
//}
//cv::Mat deblurInverseFilter(const cv::Mat& blurredImage,int R,int n, double shiftX, double shiftY, int flag) {
//     Apply Fourier transform to the blurred image
//    cv::Mat blurredImageFFT;
//    cv::Mat Butter_worth =createButterworthFilter(256, 256, R, n, true);
//    cv::dft(blurredImage, blurredImageFFT, cv::DFT_COMPLEX_OUTPUT);
//    fftShift(blurredImageFFT);
//     Create the inverse filter with known shifts
//    int rows = blurredImage.rows;
//    int cols = blurredImage.cols;
//    cv::Mat inverseFilter = Mat::zeros(rows, cols, CV_64FC2);
//    
//    for (int u = 0; u < rows; u++) {
//        for (int v = 0; v < cols; v++) {
//            double com = 0;
//            com = u * shiftX + v * shiftY;
//            auto M_I_PI_DL = (-3.14159265358979323846i *com);
//            complex<double> complexValue = (1/(PI*com))*sin(com*PI)*exp(M_I_PI_DL);
//            inverseFilter.at<complex<double>>(u, v) = complexValue;
//            cout << Butter_worth.at < complex<double>>(u, v) << endl;
//        }
//    }
//    mulSpectrums(inverseFilter, Butter_worth, inverseFilter, 0);
//    inverseFilter = inverseFilter.mul(Butter_worth);
//     Apply the inverse filter in the frequency domain
//    cv::Mat deblurredImageFFT(256,256,CV_64FC2);
//
//    for (int i = 0; i < 256; i++) {
//        for (int j = 0; j < 256; j++) {
//            deblurredImageFFT.at<complex<double>>(i, j) = blurredImageFFT.at<complex<double>>(i, j) / inverseFilter.at<complex<double>>(i, j);
//        }
//    }
//
//    if (flag) {
//    for (int u = 0; u < rows; u++) {
//        for (int v = 0; v < cols; v++) {
//            cout << inverseFilter.at<complex<double>>(u, v) << endl;
//            
//            if (abs(inverseFilter.at<complex<double>>(u, v).real()) < 0.0001 || isnan(inverseFilter.at<complex<double>>(u, v).real())) {
//                deblurredImageFFT.at<complex<double>>(u, v) = blurredImageFFT.at<complex<double>>(u, v);
//            }
//            
//        }
//    }
//    }
//    else {
//        for (int u = 0; u < rows; u++) {
//            for (int v = 0; v < cols; v++) {
//                cout << inverseFilter.at<complex<double>>(u, v) << endl;
//                
//                if ( isnan(inverseFilter.at<complex<double>>(u, v).real())) {
//                    deblurredImageFFT.at<complex<double>>(u, v) = blurredImageFFT.at<complex<double>>(u, v);
//                }
//
//            }
//        }
//    }
//    deblurredImageFFT = deblurredImageFFT.mul(Butter_worth);
//    
//    mulSpectrums(deblurredImageFFT, Butter_worth, deblurredImageFFT, 0);
//    
//    idft(deblurredImageFFT, deblurredImageFFT);
//    
//    return deblurredImageFFT;
//}
//cv::Mat WienerFilter(const cv::Mat& blurredImage,double k, double shiftX, double shiftY) {
//     Apply Fourier transform to the blurred image
//    cv::Mat blurredImageFFT;
//    cv::Mat Butter_worth = createButterworthFilter(256, 256, 10, 3, true);
//    cv::dft(blurredImage, blurredImageFFT, cv::DFT_COMPLEX_OUTPUT);
//    fftShift(blurredImageFFT);
//     Create the inverse filter with known shifts
//    int rows = blurredImage.rows;
//    int cols = blurredImage.cols;
//    cv::Mat inverseFilter = Mat::zeros(rows, cols, CV_64FC2);
//
//    for (int u = 0; u < rows; u++) {
//        for (int v = 0; v < cols; v++) {
//            double com = 0;
//            com = u * shiftX + v * shiftY;
//            auto M_I_PI_DL = (-3.14159265358979323846i * com);
//            complex<double> complexValue = (1 / (PI * com)) * sin(com * PI) * exp(M_I_PI_DL);
//            inverseFilter.at<complex<double>>(u, v) = complexValue;
//        }
//    }
//
//     Apply the inverse filter in the frequency domain
//    cv::Mat deblurredImageFFT(256, 256, CV_64FC2);
//   
//    for (int i = 0; i < 256; i++) {
//        for (int j = 0; j < 256; j++) {
//            double tmp =abs(pow(inverseFilter.at<complex<double>>(i, j).real(), 2)) + abs(pow(inverseFilter.at<complex<double>>(i, j).imag(), 2));
//            cout << tmp << endl;
//            deblurredImageFFT.at<complex<double>>(i, j) = blurredImageFFT.at<complex<double>>(i, j) * (1.0 / inverseFilter.at<complex<double>>(i, j))*(tmp/(tmp+k));
//        }
//    }
//  
//    
//    for (int u = 0; u < rows; u++) {
//        for (int v = 0; v < cols; v++) {
//            abs(inverseFilter.at<complex<double>>(u, v).real()) < 0.0001 || 
//            if (isnan(inverseFilter.at<complex<double>>(u, v).real())) {
//                deblurredImageFFT.at<complex<double>>(u, v) = blurredImageFFT.at<complex<double>>(u, v);
//            }
//
//        }
//    }
//    
//    mulSpectrums(deblurredImageFFT, Butter_worth, deblurredImageFFT, 0);
//
//    idft(deblurredImageFFT, deblurredImageFFT);
//    
//    
//    return deblurredImageFFT;
//}
//cv::Mat constrained_lsfilter(const cv::Mat& blurredImage, double gamma, double shiftX, double shiftY) {
//     Apply Fourier transform to the blurred image
//    cv::Mat blurredImageFFT;
//    cv::Mat Butter_worth = createButterworthFilter(256, 256, 10, 3, true);
//    dft(blurredImage, blurredImageFFT, cv::DFT_COMPLEX_OUTPUT);
//    fftShift(blurredImageFFT);
//     Create the inverse filter with known shifts
//    int rows = blurredImage.rows;
//    int cols = blurredImage.cols;
//    Mat inverseFilter = Mat::zeros(rows, cols, CV_64FC2);
//    Mat Lapacian = createLaplacian(256, 256);
//    Mat LapacianDFT;
//    dft(Lapacian, LapacianDFT, cv::DFT_COMPLEX_OUTPUT);
//
//    for (int u = 0; u < rows; u++) {
//        for (int v = 0; v < cols; v++) {
//            double com = 0;
//            com = u * shiftX + v * shiftY;
//            auto M_I_PI_DL = (-3.14159265358979323846i * com);
//            complex<double> complexValue = (1 / (PI * com)) * sin(com * PI) * exp(M_I_PI_DL);
//            inverseFilter.at<complex<double>>(u, v) = complexValue;
//        }
//    }
//
//     Apply the inverse filter in the frequency domain
//    Mat deblurredImageFFT(256, 256, CV_64FC2);
//    
//
//    for (int i = 0; i < 256; i++) {
//        for (int j = 0; j < 256; j++) {
//            
//            double tmp = abs(pow(inverseFilter.at<complex<double>>(i, j).real(), 2)) + abs(pow(inverseFilter.at<complex<double>>(i, j).imag(), 2));
//            double p = abs(pow(LapacianDFT.at<complex<double>>(i, j).real(), 2)) + abs(pow(LapacianDFT.at<complex<double>>(i, j).imag(), 2));
//            cout << conj(inverseFilter.at<complex<double>>(i, j)) << endl;
//            blurredImageFFT.at<complex<double>>(i, j) /= 65526.0;
//            deblurredImageFFT.at<complex<double>>(i, j) = blurredImageFFT.at<complex<double>>(i, j) * conj(inverseFilter.at<complex<double>>(i, j)) / (tmp + gamma * p);
//        }
//    }
//
//
//    for (int u = 0; u < rows; u++) {
//        for (int v = 0; v < cols; v++) {
//            (abs(inverseFilter.at<complex<double>>(u, v).real()) < 0.0001 ||
//            if  (isnan(inverseFilter.at<complex<double>>(u, v).real())) {
//                deblurredImageFFT.at<complex<double>>(u, v) = blurredImageFFT.at<complex<double>>(u, v);
//            }
//
//        }
//    }
//
//    mulSpectrums(deblurredImageFFT, Butter_worth, deblurredImageFFT, 0);
//
//    idft(deblurredImageFFT, deblurredImageFFT);
//
//    return deblurredImageFFT;
//}
//void my_main(vector<vector<double>>src,int filter,const char filename[] ,double R,double n,double k ,double gamma,int flag) {
//    Mat lena_256_blur(256, 256, CV_64F);
//    vector<vector<double>>target(256, vector<double>(256));
//
//    for (int i = 0; i < 256; i++) {
//        for (int j = 0; j < 256; j++) {
//            lena_256_blur.at<double>(i, j) = src[i][j];
//        }
//    }
//
//    Mat deblurredImage;
//    if (filter == 1) {
//        deblurredImage = deblurInverseFilter(lena_256_blur, R, n, -0.07, 0.07,flag);
//    }
//    else if (filter == 2) {
//        deblurredImage = WienerFilter(lena_256_blur, k, -0.07, 0.07);
//    }
//    else if (filter == 3) {
//        deblurredImage = constrained_lsfilter(lena_256_blur, gamma, -0.07, 0.07);
//    }
//
//    lena_256_blur.convertTo(lena_256_blur, CV_8U);
//    Mat Final_image(256, 256, CV_64F);
//
//    for (int i = 0; i < 256; i++) {
//        for (int j = 0; j < 256; j++) {
//            Final_image.at<double>(i, j) = round(sqrt(pow(deblurredImage.at<complex<double>>(i, j).real(), 2) + pow(deblurredImage.at<complex<double>>(i, j).imag(), 2)));
//        }
//    }
//    normalize(Final_image, Final_image, 0, 255, NORM_MINMAX);
//    Final_image /= 65525.0;
//    Final_image.convertTo(Final_image, CV_8U);
//    for (int i = 0; i < 256; i++) {
//        for (int j = 0; j < 256; j++) {
//            cout << Final_image.at<double>(i, j) << endl;
//            target[i][j] = double(Final_image.at<uchar>(i, j));
//
//        }
//    }
//     Display the original, blurred, and deblurred images
//    comparee(target, src, filename);
//    cv::imwrite("image_file" + string(filename) +".png", Final_image);
//    cv::imshow("Deblurred Image", Final_image);
//    cv::waitKey(0);
//}
//int main() 
//{
//     Read the input image
//    vector<vector<double>>mlena_256_blur, mlena_256_blur_noise;
//    
//
//    mlena_256_blur = read_raw("src/lena_256_blur.raw", 256, 256);
//    mlena_256_blur_noise = read_raw("src/lena_256_blur_noise.raw", 256, 256);
//    int filterselect[] = { 1,2,3 };
//    int Radius[] = { 30,40,50,60,70,80 };
//    double ka[] = { 3,0.3,0.03,0.003,0.0003,0.00003 };\
//    double Gamma[] = { 1.5,0.15,0.015,0.0015,0.00015,0.000015 };
//    for (int R : Radius) {
//        const char filename[] = "/2.1/2_1_limit_";
//        my_main(mlena_256_blur, 1, (filename + to_string(R)).c_str(), R, 2, 0.00003, 0.00003,1);
//        const char filename2[] = "/2.1/2_1_nolimit_";
//        my_main(mlena_256_blur, 1, (filename2 + to_string(R)).c_str(), R, 2, 0.00003, 0.00003, 0);
//    }
//    for (int R : Radius) {
//        const char filename[] = "/2.1/2_1_limit_noise_";
//        my_main(mlena_256_blur_noise, 1, (filename + to_string(R)).c_str(), R, 2, 0.00003, 0.00003, 1);
//        const char filename2[] = "/2.1/2_1_nolimit_noise_";
//        my_main(mlena_256_blur_noise, 1, (filename2 + to_string(R)).c_str(), R, 2, 0.00003, 0.00003, 0);
//    }
//    for (double k : ka) {
//        const char filename[] = "/2.2/2_2_";
//        my_main(mlena_256_blur, 2, (filename + to_string(k)).c_str(), 1, 2, k, 0.00003, 1);
//        const char filename2[] = "/2.2/2_2_noise_";
//        my_main(mlena_256_blur_noise, 2, (filename2 + to_string(k)).c_str(), 1, 2, k, 0.00003, 1);
//        
//    }
//    
//    for (double gamma :Gamma) {
//        const char filename[] = "/2.3/2_3_";
//        my_main(mlena_256_blur, 3, (filename + to_string(gamma)).c_str(), 1, 2, 0.00003, gamma, 1);
//        const char filename2[] = "/2.3/2_3_noise_";
//        my_main(mlena_256_blur_noise, 3, (filename2 + to_string(gamma)).c_str(), 1, 2, 0.000003, gamma, 1);
//    }
//    return 0;
//}
//
