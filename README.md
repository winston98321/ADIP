# ADIP
NTUT 2023_Fall_ADIP_hw
# Advanced Digital Image Processing Project

## Introduction
This project is part of the Advanced Digital Image Processing (ADIP) course. I used pure C++ to implement the algorithms, but in some instances, I utilized OpenCV's built-in cod in professor's require. 
Please give a star if you think it's helpful

## Note to Beginners
I hope that open-sourcing this project can help beginners in C++ to get an introduction and overcome the most challenging times. I understand the frustrations that can occur while taking this course. However, I do not condone mindless copying, as it will not help you in the long run. It may lead to difficulties in completing your final project.



## Project Contents

## Homework 1: Raw Image I/O and OpenCV

### 1. Raw Image I/O 
1.1 Raw image file format understanding
   - Download and view raw images using a preferred viewer (e.g., Xnview)

1.2 Raw image file input/output (without OpenCV)
   - Read and manipulate raw image files
   - Tasks include:
     - Reading specific pixel values
     - Dividing the image into triangular sub-blocks and rotating them
     - Partitioning the image and applying pixel value modifications

1.3 Brightness adjustment
   - Increase image brightness uniformly and randomly
   - Discussion on handling overflow/underflow issues

### 2. OpenCV Image I/O
- Set up OpenCV environment
- Create a program to draw Doraemon using OpenCV
- Add student ID to the image
- Save the result as a PNG file

## Important Notes
- Raw image manipulation should be done without using OpenCV
- Proper handling of data types and potential overflow/underflow is crucial

## Homework 2: Image Processing Fundamentals

### 1. Zooming and Shrinking
- Implement various image resizing techniques without using OpenCV
- Tasks include:
  - Zoom lena256.raw with 2:1 ratio using row-column replication
  - Shrink lena512.raw with 1:4 ratio using row-column deletion
  - Resize lena128.raw using nearest neighbor, bilinear, and bicubic interpolation
  - Resize lena512.raw to 384x384 with pre-blurring
  - Compare ↑2.25↓1.5, ↓1.5↑2.25, and ↑1.5 resizing on lena256.raw
- Compare results, calculate MSE and PSNR, discuss image quality and execution time

### 2. Distance and Path
- Find shortest paths on map10x10.raw using D4, D8, and Dm distances
- Consider different gray-value roads

### 3. Gray-level Resolution
- Quantize lena256.raw and baboon256.raw from 8 bits to 1 bit
- Calculate MSE and PSNR, discuss bit rate saving

### 4. Halftone Image
- Convert lena1024.raw to a binary halftone image
- Partition into 16x16 blocks and plot circles based on average intensity
- OpenCV circle functions allowed for this problem

## Homework 3: Image Enhancement and Histograms

### 1. Bit Plane
- Binarize and resize Doraemon image from HW1
- Replace bit planes in lena256.raw with binarized Doraemon image
- Discuss visual impact of replacing different bit planes

### 2. Grey Level Transformation
- Perform log, inverse log, and power-law transformations on log512.raw
- Apply transformations to the negative image as well
- Compare and discuss results

### 3. Histogram
- Plot histograms for log512.raw and its negative
- Implement histogram equalization
- Perform histogram matching with a specified histogram

## Important Notes
- OpenCV usage is generally not allowed unless specified
- Focus on implementing algorithms from scratch

## Homework 4: Image Moments and Spatial Filtering

### 1. Central Moments 
- Calculate centroids of shapes
- Apply central moments of orders 1 to 3

### 2. Spatial Filtering
- Perform various filtering techniques:
  - Smoothing (Box and Gaussian)
  - Roberts
  - Prewitt
  - Sobel (including custom -45° & +45° version)
  - Laplacian
  - High-boost filters

### 3. Image Smoothing & Sharpening
- Extract flowers as blob
- Perform selective filtering using ROI mask

## Homework 5: Frequency Domain Processing

### 1. 2D-DFT
- Implement DFT and IDFT
- Compare with OpenCV built-in functions
- Implement DCT and IDCT
- Discuss DCT vs DFT for image compression

### 2. Filtering in Frequency Domain
- Apply Gaussian, Butterworth, and Ideal filters (LPF and HPF)

### 3. Sobel Filtering in Frequency Domain
- Implement frequency-domain Sobel filter
- Compare with spatial-domain Sobel filter

### 4. Homomorphic Filtering)
- Implement homomorphic filter
- Compare with histogram equalization

## General Notes
- Avoid using problem-related OpenCV API unless specified
- Use origin shifting and contrast enhancement for spectrum display
- Pad image borders by mirroring or replication in mask processing

## Homework 6: Image Denoising and Deblurring

### 1. Image Denoising
- Perform Adaptive Local Noise Reduction filtering
- Perform Alpha-Trimmed Mean filtering
- Compare results on various noisy images

### 2. Image Deblurring
- Perform Inverse Filtering
- Perform Wiener Filtering
- Perform Constrained Least-Square Filter
- Compare results and discuss noise sensitivity

### 3. Geometric Transformation
- Restore distorted image using warping relationship



