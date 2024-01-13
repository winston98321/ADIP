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
//void write_raw(vector<vector<double>> img, const char filename[], int sizew, int sizeh) {
//    FILE* output_file;
//    unsigned char* image = new unsigned char[sizew * sizeh];
//    output_file = fopen(filename, "wb");
//    for (int i = 0; i < sizeh; i++) {
//        for (int j = 0; j < sizew; j++) {
//            image[i * sizew + j] = static_cast<unsigned char>(int(img[i][j]));  // Convert to unsigned char
//        }
//    }
//    fwrite(image, 1, sizeh * sizew, output_file);
//    fclose(output_file);
//}
//vector<std::vector<double>> radonTransform(const std::vector<std::vector<double>>& image) {
//    
//    const int numAngles = 180; // Number of projection angles
//
//    
//    vector<vector<double>> radonResult(numAngles, vector<double>(256, 0.0));
//
//    // Iterate over each projection angle
//    for (int angle = 0; angle < numAngles; ++angle) {
//        double theta = angle * PI / numAngles; // Convert angle to radians
//
//        // Perform Radon transform for the current angle
//        for (int x = 0; x < 256; x++) {
//            double sum = 0.0;
//
//            // Sum along the projection line
//            for (int y = 0; y < 256; y++) {
//                int projectionX = static_cast<int>(x * std::cos(theta) + y * std::sin(theta));
//                if (projectionX >= 0 && projectionX < 256) {
//                    sum += image[y][projectionX];
//                }
//            }
//
//            radonResult[angle][x] = sum;
//        }
//    }
//
//    return radonResult;
//}
//struct m_Point {
//    double x;
//    double y;
//};
//void rotatePoint(m_Point& p, double angleDegrees) {
//    // Convert angle to radians
//    double angleRadians = angleDegrees * PI / 180.0;
//
//    // Perform rotation using trigonometric functions
//    double newX = p.x * cos(angleRadians) - p.y * sin(angleRadians);
//    double newY = p.x * sin(angleRadians) + p.y * cos(angleRadians);
//
//    // Update the point with the rotated coordinates
//    p.x = newX;
//    p.y = newY;
//}
//
//int main(){
//    // Read the input image
//    vector<vector<double>>chessboard_distored, cat_distored;
//    vector<vector<double>>re_chessboard(256, vector<double>(256));
//    vector<vector<double>>re_cat(256, vector<double>(256));
//    chessboard_distored = read_raw("src/chessboard_distorted_256.raw", 256, 256);
//    cat_distored = read_raw("src/cat_distorted_256.raw", 256, 256);
//    
//    for (int i = 0; i < 256; i++) {
//        for (int j = 0; j < 256; j++) {
//            double x_pos = i / 255.0, y_pos = j / 255.0;
//            int left_point_x = round(x_pos * 104.0 + 1), right_point_x = round(x_pos * 206.0 + 49);
//            int left_min_y = 999, right_max_y = -999;
//            for (int k = 0; k < 256; k++) {
//                if (k < left_min_y && chessboard_distored[left_point_x][k] != 128)left_min_y = k;
//                if (k > right_max_y && chessboard_distored[right_point_x][k] != 128)right_max_y = k;
//            }
//            
//            int final_x = round(left_point_x + y_pos * (right_point_x - left_point_x));
//            int final_y = round(left_min_y + y_pos * (right_max_y - left_min_y));
//            re_chessboard[i][j] = chessboard_distored[final_x][final_y];
//            re_cat[i][j] = cat_distored[final_x][final_y];
//        }
//        }   
//    write_raw(re_chessboard, "image_file/mychessboard.raw", 256, 256);
//    write_raw(re_cat, "image_file/mycat.raw", 256, 256);
//       // cout << chessboard_distored[105][1] << endl;
//    //cout << max_x << " " << min_x << " " << max_y << " " << min_y << endl;
//}