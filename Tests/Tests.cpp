#include "pch.h"
#include "../dftEigen/dft_eigen.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace dftTest
{
	constexpr double m_pi = 3.14159265358979323846;
	constexpr unsigned minN = 10;
	constexpr unsigned maxN = 33;
	constexpr unsigned Precision = 300;

	TEST_CLASS(dftTest)
	{
	public:

		/// <summary>
		/// test a vector from the DFT eigen basis
		/// </summary>
		/// <param name="N"></param>
		/// <param name="index"></param>
		/// <param name="lambda_index"></param>
		/// <param name="eigenMatrix"></param>
		void testRow(unsigned N, unsigned index, unsigned lambda_index, const double* eigenVector)
		{
			auto dft = calculateDFT(eigenVector, N);

			for (unsigned i = 0; i < N; ++i)
			{
				double delta = 0.0;

				// lambda = pow(-i, index) 
				switch (lambda_index)
				{
				case 0:
					// 1.0
					delta = abs(dft[i].first - eigenVector[i]);
					break;
				case 1:
					// -i
					delta = abs(dft[i].second + eigenVector[i]);
					break;
				case 2:
					// -1.0
					delta = abs(dft[i].first + eigenVector[i]);
					break;
				case 3:
					// i
					delta = abs(dft[i].second - eigenVector[i]);
					break;
				}

				if (delta >= FLT_EPSILON)
				{
					Logger::WriteMessage(("index: " + std::to_string(index)).c_str());
					Logger::WriteMessage((", i: " + std::to_string(i)).c_str());
					Logger::WriteMessage((", Delta: " + std::to_string(delta)).c_str());
				}

				Assert::IsTrue(delta < FLT_EPSILON);
			}
		}

		/// <summary>
		/// test a vector from the DFT eigen basis scaled by a scalar 
		/// </summary>
		/// <param name="N"></param>
		/// <param name="index"></param>
		/// <param name="lambda_index"></param>
		/// <param name="eigenMatrix"></param>
		void testRowScaled(const unsigned N, unsigned index, unsigned lambda_index, double scale, const double* eigenVector)
		{
			double scaledRow[maxN];

			for (unsigned i = 0; i < N; ++i)
			{
				scaledRow[i] = scale * eigenVector[i];
			}

			testRow(N, 0, lambda_index, scaledRow);
		}

		/// <summary>
		/// test the DFT eigen basis
		/// </summary>
		TEST_METHOD(DFT_EigenBasisTest)
		{

			ScaledDFT_EigenBasisTest(1.0);

			// just a random nonzero constant
			double scale = 17.2;

			ScaledDFT_EigenBasisTest(scale);
		}

		/// <summary>
		/// test the DFT eigen basis scaled by a scalar
		/// </summary>
		void ScaledDFT_EigenBasisTest(double scale)
		{
			double eigenBasis[maxN * maxN];

			for (unsigned N = minN; N <= maxN; ++N)
			{
				if (getDFTEigenVectors(eigenBasis, N, Precision) != 0)
				{
					return;
				}

				for (unsigned index = 0; index < N - 1; ++index)
				{
					unsigned lambda_index = index % 4;
					unsigned row = index * N;
					testRowScaled(N, index, lambda_index, scale, &eigenBasis[row]);
				}

				unsigned index = N - 1;
				unsigned row = index * N;
				unsigned last_index = N % 2 == 0 ? N : N - 1;
				unsigned lambda_index = last_index % 4;
				testRowScaled(N, index, lambda_index, scale, &eigenBasis[row]);
			}
		}

		/// <summary>
		/// Simple function to calculate the Discrete Fourier transform
		/// The zero frequency in the center
		/// </summary>
		/// <param name="x"></param>
		/// <param name="n"></param>
		/// <returns></returns>
		std::vector<std::pair<double, double> > calculateDFT(const double* x, const unsigned n)
		{
			std::vector<std::pair<double, double> > X(n, std::make_pair(0.0, 0.0));

			long N0 = (long)ceil(n / 2.0) - 1;
			long N1 = (long)floor(n / 2.0);
			const double pi2 = 2.0 * m_pi / n;
			const double scale = 1.0 / sqrt((double)n);

			for (long k = -N0; k <= N1; ++k)
			{
				long k1 = k + N0;
				double w = pi2 * k;
				for (long i = -N0; i <= N1; ++i)
				{
					long i1 = i + N0;
					X[k1].first += x[i1] * cos(w * i);
					X[k1].second -= x[i1] * sin(w * i);
				}

				X[k1].first *= scale;
				X[k1].second *= scale;
			}

			return X;
		}
	};
}
