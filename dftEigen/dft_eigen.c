// dftEigen.cpp : Defines the functions for the static library.

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>

static mpfr_t pi;
static mpfr_t w;
static mpfr_t w2;
static mpfr_t* S; // (2 * N + 1);
static mpfr_t one;
static unsigned N;

static void init(unsigned precision)
{
	unsigned j;

	mpfr_t s, wj, m;

	mpf_set_default_prec(precision);

	mpfr_init(pi);
	mpfr_const_pi(pi, MPFR_RNDN);

	mpfr_init(w);
	mpfr_div_ui(w, pi, N, MPFR_RNDN);

	mpfr_init(w2);
	mpfr_mul_ui(w2, w, 2, MPFR_RNDN);

	mpfr_init(one);
	mpfr_set_d(one, 1.0, MPFR_RNDN);

	mpfr_init(s);
	mpfr_init(wj);
	mpfr_init(m);

	S = malloc(sizeof(mpfr_t) * (2ull * N + 1));

	//create S - sequence
	mpfr_init(S[0]);
	mpfr_set_d(S[0], 1.0, MPFR_RNDN);
	for (j = 1; j < N; ++j)
	{
		// S[j] = S[j - 1] * 2 * sin(w * j);
		mpfr_init(S[j]);
		mpfr_mul_ui(wj, w, j, MPFR_RNDN);
		mpfr_sin(s, wj, MPFR_RNDN);
		mpfr_mul_ui(m, s, 2, MPFR_RNDN);
		mpfr_mul(S[j], S[j - 1], m, MPFR_RNDN);
	}

	for (j = N; j < 2 * N + 1; ++j)
	{
		mpfr_init(S[j]);
		mpfr_set_d(S[j], 0.0, MPFR_RNDN);
	}
}

static void clear()
{
	mpfr_clear(pi);
	mpfr_clear(w);
	mpfr_clear(one);

	for (unsigned j = 0; j < 2 * N + 1; ++j)
	{
		mpfr_clear(S[j]);
	}

	free(S);

	mpfr_free_cache();
}

void create_u(mpfr_t* u_n, unsigned m)
{
	mpfr_t alpha_n, s, wm, d;

	mpfr_init(alpha_n);
	mpfr_init(s);
	mpfr_init(wm);
	mpfr_init(d);

	//compute alpha_n
	if (N % 2 == 0)
	{
		if (m == 0)
		{
			mpfr_set_d(alpha_n, 0.5, MPFR_RNDN);
		}
		else
		{
			//alpha_n = sqrt(S[2 * m - 1] * sin(w * m)) / (S[m] * S[m]);
			mpfr_mul_ui(wm, w, m, MPFR_RNDN);
			mpfr_sin(s, wm, MPFR_RNDN);
			mpfr_sqr(d, S[m], MPFR_RNDN);
			mpfr_mul(wm, S[2 * m - 1], s, MPFR_RNDN);
			mpfr_sqrt(s, wm, MPFR_RNDN);
			mpfr_div(alpha_n, s, d, MPFR_RNDN);
		}
	}
	else
	{
		if (m == 0)
		{
			mpfr_set_d(alpha_n, 1.0, MPFR_RNDN);
		}
		else
		{
			//alpha_n = sqrt(S[2 * m]) / (S[m] * S[m]);
			mpfr_sqr(d, S[m], MPFR_RNDN);
			mpfr_sqrt(s, S[2 * m], MPFR_RNDN);
			mpfr_div(alpha_n, s, d, MPFR_RNDN);
		}
	}

	//compute vector u_n
	long N0 = (long)ceil(N / 2.0) - 1;
	long N1 = (long)floor(N / 2.0);

	if (N % 2 == 0)
	{
		if (m < N / 2)
		{
			unsigned long N2 = N * N;
			for (long k = -N0; k <= N1; ++k)
			{
				long j = k + N0 + 1 - 1;
				//u_n[j] = (alpha_n * S[m] * S[m] / (N * N)) * cos(w * k) * S[N - m - 1 - k] * S[N - m - 1 + k];
				mpfr_sqr(d, S[m], MPFR_RNDN);
				mpfr_mul(s, d, alpha_n, MPFR_RNDN);
				mpfr_div_ui(d, s, N2, MPFR_RNDN);
				mpfr_mul_si(wm, w, k, MPFR_RNDN);
				mpfr_cos(s, wm, MPFR_RNDN);
				mpfr_mul(wm, d, s, MPFR_RNDN);
				mpfr_mul(d, wm, S[N - m - 1 - k], MPFR_RNDN);
				mpfr_mul(u_n[j], d, S[N - m - 1 + k], MPFR_RNDN);
			}
		}
		else
		{
			for (long k = (-N0); k < N1; ++k)
			{
				long j = k + N0 + 1 - 1;
				//u_n[j] = 1 / (2 * sqrt(N));
				mpfr_set_ui(d, N, MPFR_RNDN);
				mpfr_sqrt(s, d, MPFR_RNDN);
				mpfr_mul_ui(d, s, 2, MPFR_RNDN);
				mpfr_div(u_n[j], one, d, MPFR_RNDN);
			}
		}
	}
	else
	{
		unsigned long N2 = N * N;
		for (long k = -N0; k <= N1; ++k)
		{
			long j = k + N0 + 1 - 1;
			//u_n[j] = (alpha_n * S[m] * S[m] / (N * N)) * S[N - m - 1 - k] * S[N - m - 1 + k];
			mpfr_sqr(d, S[m], MPFR_RNDN);
			mpfr_mul(s, d, alpha_n, MPFR_RNDN);
			mpfr_div_ui(d, s, N2, MPFR_RNDN);
			mpfr_mul(s, d, S[N - m - 1 - k], MPFR_RNDN);
			mpfr_mul(u_n[j], s, S[N - m - 1 + k], MPFR_RNDN);
		}
	}

	mpfr_clear(alpha_n);
	mpfr_clear(s);
	mpfr_clear(wm);
	mpfr_clear(d);
}

void create_v(mpfr_t* v_n, unsigned m)
{
	mpfr_t beta_n, s, wm, d;

	mpfr_init(beta_n);
	mpfr_init(s);
	mpfr_init(wm);
	mpfr_init(d);

	//compute beta_n
	if (N % 2 == 0)
	{
		//beta_n = sqrt(S[2 * m - 1] * cos(w * m)) / (S[m] * S[m]);
		mpfr_sqr(d, S[m], MPFR_RNDN);
		mpfr_mul_ui(wm, w, m, MPFR_RNDN);
		mpfr_cos(s, wm, MPFR_RNDN);
		mpfr_mul(wm, S[2 * m - 1], s, MPFR_RNDN);
		mpfr_sqrt(s, wm, MPFR_RNDN);
		mpfr_div(beta_n, s, d, MPFR_RNDN);
	}
	else
	{
		//beta_n = sqrt(S[2 * m - 1]) / (S[m] * S[m]);
		mpfr_sqr(d, S[m], MPFR_RNDN);
		mpfr_sqrt(s, S[2 * m - 1], MPFR_RNDN);
		mpfr_div(beta_n, s, d, MPFR_RNDN);
	}

	//compute vector v_n
	long N0 = (long)ceil(N / 2.0) - 1;
	long N1 = (long)floor(N / 2.0);
	unsigned long N2 = N * N;

	if (N % 2 == 0)
		for (long k = -N0; k <= N1; ++k)
		{
			long j = k + N0 + 1 - 1;
			//v_n[j] = (2 * beta_n * S[m] * S[m] / (N * N)) * sin(w * k) * S[N - m - 1 - k] * S[N - m - 1 + k];
			mpfr_sqr(d, S[m], MPFR_RNDN);
			mpfr_mul(s, d, beta_n, MPFR_RNDN);
			mpfr_div_ui(d, s, N2, MPFR_RNDN);
			mpfr_mul_ui(d, d, 2, MPFR_RNDN);
			mpfr_mul_si(wm, w, k, MPFR_RNDN);
			mpfr_sin(s, wm, MPFR_RNDN);
			mpfr_mul(wm, d, s, MPFR_RNDN);
			mpfr_mul(d, wm, S[N - m - 1 - k], MPFR_RNDN);
			mpfr_mul(v_n[j], d, S[N - m - 1 + k], MPFR_RNDN);
		}
	else
		for (long k = -N0; k <= N1; ++k)
		{
			long j = k + N0 + 1 - 1;
			
			//v_n[j] = (beta_n * S[m] * S[m] / (N * N)) * sin(w * (k * 2)) * S[N - m - 1 - k] * S[N - m - 1 + k];
			mpfr_sqr(d, S[m], MPFR_RNDN);
			mpfr_mul(s, d, beta_n, MPFR_RNDN);
			mpfr_div_ui(d, s, N2, MPFR_RNDN);
			mpfr_mul_si(wm, w2, k, MPFR_RNDN);
			mpfr_sin(s, wm, MPFR_RNDN);
			mpfr_mul(wm, d, s, MPFR_RNDN);
			mpfr_mul(d, wm, S[N - m - 1 - k], MPFR_RNDN);
			mpfr_mul(v_n[j], d, S[N - m - 1 + k], MPFR_RNDN);
		}

	mpfr_clear(beta_n);
	mpfr_clear(s);
	mpfr_clear(wm);
	mpfr_clear(d);
}

void norm(mpfr_t* f_result, mpfr_t* v)
{
	mpfr_t s;

	mpfr_init(s);

	mpfr_set_d(*f_result, 0.0, MPFR_RNDN);
	for (unsigned k = 0; k < N; ++k)
	{
		//f_result = f_result + v[k] * v[k];
		mpfr_sqr(s, v[k], MPFR_RNDN);
		mpfr_add(*f_result, *f_result, s, MPFR_RNDN);
	}

	//f_result = sqrt(f_result);
	mpfr_sqrt(*f_result, *f_result, MPFR_RNDN);

	mpfr_clear(s);
}

void create_first_four_T(mpfr_t* T0, mpfr_t* T1, mpfr_t* T2, mpfr_t* T3)
{
	//u1(n), u2(n), v1(n), v2(n);
	mpfr_t* u1 = malloc(sizeof(mpfr_t) * N);
	mpfr_t* u2 = malloc(sizeof(mpfr_t) * N);
	mpfr_t* v1 = malloc(sizeof(mpfr_t) * N);
	mpfr_t* v2 = malloc(sizeof(mpfr_t) * N);
	for (unsigned k = 0; k < N; ++k)
	{
		mpfr_init(u1[k]);
		mpfr_init(u2[k]);
		mpfr_init(v1[k]);
		mpfr_init(v2[k]);
	}

	mpfr_t c;
	mpfr_init(c);

	unsigned K0 = (long)floor((N + 2) / 4.0);
	unsigned K1 = (long)floor((N + 3) / 4.0);
	unsigned K2 = (long)floor((N + 4) / 4.0);
	unsigned K3 = (long)floor((N + 5) / 4.0);
	unsigned Nc = (long)ceil(N / 2.0);
	unsigned Nf = (long)floor(N / 2.0);

	//create T_0
	create_u(u1, K0);
	create_u(u2, Nf - K0);
	for (unsigned k = 0; k < N; ++k)
	{
		//T0[k] = u1[k] + u2[k];
		mpfr_add(T0[k], u1[k], u2[k], MPFR_RNDN);
	}

	norm(&c, T0);
	for (unsigned k = 0; k < N; ++k)
	{
		//T0[k] = T0[k] / c;
		mpfr_div(T0[k], T0[k], c, MPFR_RNDN);
	}

	//create T_1
	create_v(v1, K1);
	create_v(v2, Nc - K1);
	for (unsigned k = 0; k < N; ++k)
	{
		//T1[k] = v1[k] + v2[k];
		mpfr_add(T1[k], v1[k], v2[k], MPFR_RNDN);
	}

	norm(&c, T1);
	for (unsigned k = 0; k < N; ++k)
	{
		//T1[k] = T1[k] / c;
		mpfr_div(T1[k], T1[k], c, MPFR_RNDN);
	}

	//create T_2
	create_u(u1, K2);
	create_u(u2, Nf - K2);
	for (unsigned k = 0; k < N; ++k)
	{
		//T2[k] = u1[k] - u2[k];
		mpfr_sub(T2[k], u1[k], u2[k], MPFR_RNDN);
	}

	norm(&c, T2);
	for (unsigned k = 0; k < N; ++k)
	{
		//T2[k] = T2[k] / c;
		mpfr_div(T2[k], T2[k], c, MPFR_RNDN);
	}

	//create T_3
	create_v(v1, K3);
	create_v(v2, Nc - K3);
	for (unsigned k = 0; k < N; ++k)
	{
		//T3[k] = v1[k] - v2[k];
		mpfr_sub(T3[k], v1[k], v2[k], MPFR_RNDN);
	}
	norm(&c, T3);
	for (unsigned k = 0; k < N; ++k)
	{
		//T3[k] = T3[k] / c;
		mpfr_div(T3[k], T3[k], c, MPFR_RNDN);
	}

	free(u1);
	free(u2);
	free(v1);
	free(v2);

	mpfr_clear(c);
}

void L_operator(mpfr_t* Lv, mpfr_t* v)
{
	mpfr_t tmp1, tmp2;
	mpfr_init(tmp1);
	mpfr_init(tmp2);

	long N0 = (long)ceil(N / 2.0) - 1;

	//Lv[1] = v[2] + v[N] + 2 * cos(w * k) * v[1];
	mpfr_mul_si(tmp1, w2, -N0, MPFR_RNDN);
	mpfr_cos(tmp2, tmp1, MPFR_RNDN);
	mpfr_mul_ui(tmp1, tmp2, 2, MPFR_RNDN);
	mpfr_mul(tmp2, tmp1, v[1 - 1], MPFR_RNDN);
	mpfr_add(tmp1, tmp2, v[N - 1], MPFR_RNDN);
	// Lv[0]
	mpfr_add(Lv[1 - 1], tmp1, v[2 - 1], MPFR_RNDN);

	long N1 = (long)floor(N / 2.0);

	//Lv[N] = v[N - 1] + v[1] + 2 * cos(w * k) * v[N];
	mpfr_mul_si(tmp1, w2, N1, MPFR_RNDN);
	mpfr_cos(tmp2, tmp1, MPFR_RNDN);
	mpfr_mul_ui(tmp1, tmp2, 2, MPFR_RNDN);
	mpfr_mul(tmp2, tmp1, v[N - 1], MPFR_RNDN);
	mpfr_add(tmp1, tmp2, v[1 - 1], MPFR_RNDN);
	// Lv[N - 1]
	mpfr_add(Lv[N - 1], tmp1, v[N - 1 - 1], MPFR_RNDN);

	for (long k = -N0 + 1; k <= N1 - 1; ++k)
	{
		long j = k + N0 + 1 - 1;
		//Lv(j) = v(j + 1) + v(j - 1) + 2 * cos(w * k) * v(j);
		mpfr_mul_si(tmp1, w2, k, MPFR_RNDN);
		mpfr_cos(tmp2, tmp1, MPFR_RNDN);
		mpfr_mul_ui(tmp1, tmp2, 2, MPFR_RNDN);
		mpfr_mul(tmp2, tmp1, v[j], MPFR_RNDN);
		mpfr_add(tmp1, tmp2, v[j - 1], MPFR_RNDN);
		// Lv[j]
		mpfr_add(Lv[j], tmp1, v[j + 1], MPFR_RNDN);
	}

	mpfr_clear(tmp2);
	mpfr_clear(tmp1);
}

#define vectorMinimumSize 8
#define vectorMaximumSize 2048
#define precisionMinimum 32
#define precisionMaximum 1024

int getDFTEigenVectors(double* eigenMatrix, const unsigned n, const unsigned precision)
{
	if (n < vectorMinimumSize || n > vectorMaximumSize || precision > precisionMaximum || precision < precisionMinimum)
	{
		return -1;
	}

	N = n;

	init(precision);

	// compute vectors T_0, T_1, T_2, T_3
	mpfr_t* T0 = malloc(sizeof(mpfr_t) * N);
	mpfr_t* T1 = malloc(sizeof(mpfr_t) * N);
	mpfr_t* T2 = malloc(sizeof(mpfr_t) * N);
	mpfr_t* T3 = malloc(sizeof(mpfr_t) * N);

	for (unsigned i = 0; i < N; ++i)
	{
		mpfr_init(T0[i]);
		mpfr_init(T1[i]);
		mpfr_init(T2[i]);
		mpfr_init(T3[i]);
	}

	create_first_four_T(T0, T1, T2, T3);

	// result T0,..,T3

	for (unsigned i = 0; i < n; ++i)
	{
		eigenMatrix[i] = mpfr_get_d1(T0[i]);
		eigenMatrix[i + N] = mpfr_get_d1(T1[i]);
		eigenMatrix[i + 2 * N] = mpfr_get_d1(T2[i]);
		eigenMatrix[i + 3 * N] = mpfr_get_d1(T3[i]);
	}

	// compute vectors T_4, T_5, ...., T_{ N - 2 }
	mpfr_t a0, b0, b_4, b_3, b_2, b_1;
	mpfr_t tmp1, tmp2;
	mpfr_init(tmp1);
	mpfr_init(tmp2);

	mpfr_init(b0);
	mpfr_set_d(b0, 0.0, MPFR_RNDN);
	mpfr_init(b_1);
	mpfr_set_d(b_1, 0.0, MPFR_RNDN);
	mpfr_init(b_2);
	mpfr_set_d(b_2, 0.0, MPFR_RNDN);
	mpfr_init(b_3);
	mpfr_set_d(b_3, 0.0, MPFR_RNDN);
	mpfr_init(b_4);
	mpfr_set_d(b_4, 0.0, MPFR_RNDN);
	mpfr_init(a0);
	mpfr_set_d(a0, 0.0, MPFR_RNDN);

	// v(n)
	mpfr_t* v = malloc(sizeof(mpfr_t) * N);

	// Tm1(n), Tm2(n), Tm3(n), Tm4(n), Tnew(n);
	mpfr_t* Tm1 = malloc(sizeof(mpfr_t) * N);
	mpfr_t* Tm2 = malloc(sizeof(mpfr_t) * N);
	mpfr_t* Tm3 = malloc(sizeof(mpfr_t) * N);
	mpfr_t* Tm4 = malloc(sizeof(mpfr_t) * N);
	mpfr_t* Tnew = malloc(sizeof(mpfr_t) * N);
	for (unsigned i = 0; i < N; ++i)
	{
		mpfr_init(v[i]);
		mpfr_init(Tm1[i]);
		mpfr_set_d(Tm1[i], 0.0, MPFR_RNDN);
		mpfr_init(Tm2[i]);
		mpfr_set_d(Tm2[i], 0.0, MPFR_RNDN);
		mpfr_init(Tm3[i]);
		mpfr_set_d(Tm3[i], 0.0, MPFR_RNDN);
		mpfr_init(Tm4[i]);
		mpfr_set_d(Tm4[i], 0.0, MPFR_RNDN);
		mpfr_init(Tnew[i]);
	}

	for (unsigned l = 4; l < N - 1; ++l)
	{
		L_operator(v, T0);
		mpfr_set_d(a0, 0.0, MPFR_RNDN);
		for (unsigned k = 0; k < N; ++k)
		{
			// a0 = a0 + T0(k) * v(k);
			mpfr_mul(tmp1, T0[k], v[k], MPFR_RNDN);
			mpfr_add(a0, a0, tmp1, MPFR_RNDN);
		}
		mpfr_set_d(b0, 0.0, MPFR_RNDN);
		for (unsigned k = 0; k < N; ++k)
		{
			//Tnew(k) = v(k) - a0 * T0(k) - b(-4) * Tm4(k)
			mpfr_mul(tmp1, b_4, Tm4[k], MPFR_RNDN);
			mpfr_sub(tmp2, v[k], tmp1, MPFR_RNDN);
			mpfr_mul(tmp1, a0, T0[k], MPFR_RNDN);
			mpfr_sub(Tnew[k], tmp2, tmp1, MPFR_RNDN);

			//b0 = b0 + Tnew(k) * *2
			mpfr_sqr(tmp1, Tnew[k], MPFR_RNDN);
			mpfr_add(b0, b0, tmp1, MPFR_RNDN);
		}
		mpfr_sqrt(b0, b0, MPFR_RNDN);
		for (unsigned k = 0; k < N; ++k)
		{
			//Tnew(k) = Tnew(k) / b0
			mpfr_div(Tnew[k], Tnew[k], b0, MPFR_RNDN);
		}
		// result Tnew(j)
		
		for (unsigned i = 0; i < n; ++i)
		{
			eigenMatrix[i + l * N] = mpfr_get_d1(Tnew[i]);
		}

		// do i = -4, -2
		// b(i) = b(i + 1)
		// end do

		mpfr_set(b_4, b_3, MPFR_RNDN);
		mpfr_set(b_3, b_2, MPFR_RNDN);
		mpfr_set(b_2, b_1, MPFR_RNDN);
		mpfr_set(b_1, b0, MPFR_RNDN);

		//Tm4 = Tm3
		//Tm3 = Tm2
		//Tm2 = Tm1
		//Tm1 = T0
		//T0 = T1
		//T1 = T2
		//T2 = T3
		//T3 = Tnew

		for (unsigned k = 0; k < N; ++k)
		{
			mpfr_set(Tm4[k], Tm3[k], MPFR_RNDN);
			mpfr_set(Tm3[k], Tm2[k], MPFR_RNDN);
			mpfr_set(Tm2[k], Tm1[k], MPFR_RNDN);
			mpfr_set(Tm1[k], T0[k], MPFR_RNDN);
			mpfr_set(T0[k], T1[k], MPFR_RNDN);
			mpfr_set(T1[k], T2[k], MPFR_RNDN);
			mpfr_set(T2[k], T3[k], MPFR_RNDN);
			mpfr_set(T3[k], Tnew[k], MPFR_RNDN);
		}
	}

	// compute the remaining vector T_{ N - 1 }
	if (N % 2 == 1)
	{
		L_operator(v, T0);
		mpfr_set_d(a0, 0.0, MPFR_RNDN);
		for (unsigned k = 0; k < N; ++k)
		{
			// a0 = a0 + T0(k) * v(k);
			mpfr_mul(tmp1, T0[k], v[k], MPFR_RNDN);
			mpfr_add(a0, a0, tmp1, MPFR_RNDN);
		}
		mpfr_set_d(b0, 0.0, MPFR_RNDN);
		for (unsigned k = 0; k < N; ++k)
		{
			//Tnew(k) = v(k) - a0 * T0(k) - b(-4) * Tm4(k)
			mpfr_mul(tmp1, b_4, Tm4[k], MPFR_RNDN);
			mpfr_sub(tmp2, v[k], tmp1, MPFR_RNDN);
			mpfr_mul(tmp1, a0, T0[k], MPFR_RNDN);
			mpfr_sub(Tnew[k], tmp2, tmp1, MPFR_RNDN);

			//b0 = b0 + Tnew(k) * *2
			mpfr_sqr(tmp1, Tnew[k], MPFR_RNDN);
			mpfr_add(b0, b0, tmp1, MPFR_RNDN);
		}
		mpfr_sqrt(b0, b0, MPFR_RNDN);
		for (unsigned k = 0; k < N; ++k)
		{
			//Tnew(k) = Tnew(k) / b0
			mpfr_div(Tnew[k], Tnew[k], b0, MPFR_RNDN);
		}

		// result Tnew[j]
		for (unsigned i = 0; i < n; ++i)
		{
			eigenMatrix[i + (N - 1) * N] = mpfr_get_d1(Tnew[i]);
		}
	}
	else
	{
		L_operator(v, T1);
		mpfr_set_d(a0, 0.0, MPFR_RNDN);
		for (unsigned k = 0; k < N; ++k)
		{
			// a0 = a0 + T1(k) * v(k);
			mpfr_mul(tmp1, T1[k], v[k], MPFR_RNDN);
			mpfr_add(a0, a0, tmp1, MPFR_RNDN);
		}
		mpfr_set_d(b0, 0.0, MPFR_RNDN);
		for (unsigned k = 0; k < N; ++k)
		{
			//Tnew(k) = v(k) - a0 * T1(k) - b(-4) * Tm3(k)
			mpfr_mul(tmp1, b_4, Tm3[k], MPFR_RNDN);
			mpfr_sub(tmp2, v[k], tmp1, MPFR_RNDN);
			mpfr_mul(tmp1, a0, T1[k], MPFR_RNDN);
			mpfr_sub(Tnew[k], tmp2, tmp1, MPFR_RNDN);

			//b0 = b0 + Tnew(k) * *2
			mpfr_sqr(tmp1, Tnew[k], MPFR_RNDN);
			mpfr_add(b0, b0, tmp1, MPFR_RNDN);
		}
		mpfr_sqrt(b0, b0, MPFR_RNDN);
		for (unsigned k = 0; k < N; ++k)
		{
			//Tnew(k) = Tnew(k) / b0
			mpfr_div(Tnew[k], Tnew[k], b0, MPFR_RNDN);
		}

		// result Tnew[j]
		for (unsigned i = 0; i < n; ++i)
		{
			eigenMatrix[i + (N - 1) * N] = mpfr_get_d1(Tnew[i]);
		}
	}

	for (unsigned i = 0; i < N; ++i)
	{
		mpfr_clear(v[i]);
		mpfr_clear(Tm1[i]);
		mpfr_clear(Tm2[i]);
		mpfr_clear(Tm3[i]);
		mpfr_clear(Tm4[i]);
		mpfr_clear(Tnew[i]);
		mpfr_clear(T0[i]);
		mpfr_clear(T1[i]);
		mpfr_clear(T2[i]);
		mpfr_clear(T3[i]);
	}

	free(Tm4);
	free(Tm3);
	free(Tm2);
	free(Tm1);
	free(Tnew);
	free(v);

	mpfr_clear(b0);
	mpfr_clear(b_1);
	mpfr_clear(b_2);
	mpfr_clear(b_3);
	mpfr_clear(b_4);
	mpfr_clear(a0);
	mpfr_clear(tmp1);
	mpfr_clear(tmp2);

	free(T3);
	free(T2);
	free(T1);
	free(T0);

	clear();

	return 0;
}

