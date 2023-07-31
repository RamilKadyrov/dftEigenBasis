#pragma once
// This is a library to calculate Eigenbasis of the Discrete Fourier Transform
// based and Fortran code from the article
// @article{Kuznetsov2017MinimalHE,
// title = { Minimal Hermite - Type Eigenbasis of the Discrete Fourier Transform },
// author = { Alexey Kuznetsov and Mateusz Kwasnicki },
// journal = { Journal of Fourier Analysis and Applications },
// year = { 2017 },
// volume = { 25 },
// pages = { 1053 - 1079 },
// url = { https://api.semanticscholar.org/CorpusID:119165157}
// }

#ifdef __cplusplus
extern "C" {
#endif
	int getDFTEigenVectors(double* eigenVectors, const unsigned n, const unsigned precision);
#ifdef __cplusplus
}
#endif