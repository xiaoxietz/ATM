#include "FHEW.h"

int FHEW_new_Pk_a0[Polynomial_def_N] = { 0 };
int FHEW_calced_Pk_b0_b1[2 * Polynomial_def_N] = { 0 };

int Sample(const Distrib& Chi) 
{
	if (Chi.max) 
	{
		double r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		for (int i = 0; i < Chi.max; ++i)
		{
			if (r <= Chi.table[i])
			{
				return i - Chi.offset;
			}
		}
		cout << "Sampling Error: distribution table ending before (double) 1.0" << endl;
		exit(1);
	}

	double r, s = Chi.std_dev;
	if (s < 500)
	{
		int x, maxx = ceil(s * 8);
		while (true) 
		{
			x = rand() % (2 * maxx + 1) - maxx;
			r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			if (r < exp(-x*x / (2 * s*s))) return x;
		}
	}
	// For some reason unknown to us, the previous implementation provides a bad distribution for large s...
	// We switch from "discrete gaussian" to rounded gaussian when s gets larger
	double x;

	while (true) 
	{
		x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		x = 16 * x - 8;
		r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		if (r < exp(-x*x / 2)) return floor(.5 + x*s);
	}
}

void FHEW_Key_gen(int *FHEW_sk)
{
	for (int i = 0; i < Polynomial_def_N; ++i)
	{
		FHEW_sk[i] = Sample(Chi1);
	}
}

void FHEW_Encrypt(int plain_text, int *FHEW_sk, bool Pk_a0_save_flag, bool Pk_b0_b1_save_flag, int(*FHEW_ct)[FHEW_columns], int plain_text2)
{
	int Plain_text_Mod_q = (((plain_text % FHEW_Mod_q) + FHEW_Mod_q) % FHEW_Mod_q); /* * (2 * Polynomial_def_N / FHEW_Mod_q);*/// Reduce mod q (dealing with negative number as well)
	int sign = 1;
	int FHEW_result_ciphertext[Polynomial_def_N] = { 0 };

	if (Plain_text_Mod_q >= Polynomial_def_N)
	{ 
		Plain_text_Mod_q -= Polynomial_def_N;
		sign = -1; 
	}

	for (int i = 0; i < FHEW_rows; ++i)
	{
		for (int j = 0; j < Polynomial_def_N; ++j)
		{
			FHEW_ct[i * Polynomial_def_N + j][0] = rand() % FHEW_Mod_q;
			if (Pk_a0_save_flag)
			{
				FHEW_new_Pk_a0[j] = FHEW_ct[i * Polynomial_def_N + j][0];
			}
			if (Pk_b0_b1_save_flag)
			{
				FHEW_calced_Pk_b0_b1[i * Polynomial_def_N + j] = FHEW_ct[i * Polynomial_def_N + j][0];
			}
		}

		Pk_a0_save_flag = false;
		if (i == 1)
		{
			Pk_b0_b1_save_flag = false;
		}

		Polynomial_Mult(&FHEW_ct[i * Polynomial_def_N], 0, FHEW_sk, FHEW_result_ciphertext);//as
		for (int k = 0; k < Polynomial_def_N; ++k)
		{
			FHEW_ct[i * Polynomial_def_N + k][1] = FHEW_result_ciphertext[k] + Sample(Chi1);//as+e
		}
		memset(FHEW_result_ciphertext, 0, sizeof(FHEW_result_ciphertext));
	}

	for (int i = 0; i < FHEW_rows / 2; ++i)
	{
		FHEW_ct[2 * i * Polynomial_def_N + Plain_text_Mod_q][0] =
			((((FHEW_ct[2 * i * Polynomial_def_N + Plain_text_Mod_q][0] + sign * vgprime[i]) % FHEW_Mod_q) + FHEW_Mod_q) % FHEW_Mod_q); // Add G Multiple
		FHEW_ct[(2 * i + 1) * Polynomial_def_N + Plain_text_Mod_q][1] =
			((((FHEW_ct[(2 * i + 1) * Polynomial_def_N + Plain_text_Mod_q][1] + sign * vgprime[i]) % FHEW_Mod_q) + FHEW_Mod_q) % FHEW_Mod_q); // [a,as+e] + X^m *G

		if (plain_text2 >=0)
		{
			FHEW_ct[2 * i * Polynomial_def_N + plain_text2][0] =
				((((FHEW_ct[2 * i * Polynomial_def_N + plain_text2][0] + sign * vgprime[i]) % FHEW_Mod_q) + FHEW_Mod_q) % FHEW_Mod_q);
			FHEW_ct[(2 * i + 1) * Polynomial_def_N + plain_text2][1] =
				((((FHEW_ct[(2 * i + 1) * Polynomial_def_N + plain_text2][1] + sign * vgprime[i]) % FHEW_Mod_q) + FHEW_Mod_q) % FHEW_Mod_q);
		}
	}
}

void Polynomial_Mult(int(*FHEW_ct1)[FHEW_columns], int ct_columns, int *FHEW_ct2, int *FHEW_result_ct)
{
	int Temporary_array[Polynomial_def_N][Polynomial_def_N] = { 0 };
	int postion_displacement = 0;

	for (int i = 0; i < Polynomial_def_N; ++i)
	{
		for (int j = 0; j < Polynomial_def_N; ++j)
		{
			postion_displacement = j + i;
			if (postion_displacement > Polynomial_def_N - 1)
			{
				Temporary_array[postion_displacement % Polynomial_def_N][i] = -FHEW_ct1[j][ct_columns];
			}
			else
			{
				Temporary_array[postion_displacement][i] = FHEW_ct1[j][ct_columns];
			}
		}
	}

	for (int m = 0; m < Polynomial_def_N; m++)
	{
		for (int n = 0; n < Polynomial_def_N; n++)
		{
			FHEW_result_ct[m] += Temporary_array[m][n] * FHEW_ct2[n];
		}
		FHEW_result_ct[m] = (((FHEW_result_ct[m] % FHEW_Mod_q) + FHEW_Mod_q) % FHEW_Mod_q);
	}
}

void FHEW_Ciphertext_Mult(int(*FHEW_ct1)[FHEW_columns], int(*FHEW_ct2)[FHEW_columns], int(*FHEW_ct_result)[FHEW_columns])
{
	int Polynomial_Mult_Result1[Polynomial_def_N] = { 0 };
	int Polynomial_Mult_Result2[Polynomial_def_N] = { 0 };
	int Polynomial_Mult_Result3[Polynomial_def_N] = { 0 };
	int Polynomial_Mult_Result4[Polynomial_def_N] = { 0 };
	int Temporary_FHEW_ct2_0_0[Polynomial_def_N] = { 0 };
	int Temporary_FHEW_ct2_1_0[Polynomial_def_N] = { 0 };
	int Temporary_FHEW_ct2_0_1[Polynomial_def_N] = { 0 };
	int Temporary_FHEW_ct2_1_1[Polynomial_def_N] = { 0 };
	int k = 0;

	for (k = 0; k < Polynomial_def_N; ++k)
	{
		Temporary_FHEW_ct2_0_0[k] = FHEW_ct2[0 * Polynomial_def_N + k][0];
		Temporary_FHEW_ct2_1_0[k] = FHEW_ct2[1 * Polynomial_def_N + k][0];
		Temporary_FHEW_ct2_0_1[k] = FHEW_ct2[0 * Polynomial_def_N + k][1];
		Temporary_FHEW_ct2_1_1[k] = FHEW_ct2[1 * Polynomial_def_N + k][1];
	}

	for (int i = 0; i < FHEW_rows; ++i)
	{
		Polynomial_Mult(&FHEW_ct1[i * Polynomial_def_N], 0, Temporary_FHEW_ct2_0_0, Polynomial_Mult_Result1);
		Polynomial_Mult(&FHEW_ct1[i * Polynomial_def_N], 1, Temporary_FHEW_ct2_1_0, Polynomial_Mult_Result2);
		Polynomial_Mult(&FHEW_ct1[i * Polynomial_def_N], 0, Temporary_FHEW_ct2_0_1, Polynomial_Mult_Result3);
		Polynomial_Mult(&FHEW_ct1[i * Polynomial_def_N], 1, Temporary_FHEW_ct2_1_1, Polynomial_Mult_Result4);
		
		for (int j = 0; j < Polynomial_def_N; ++j)
		{
			FHEW_ct_result[i * Polynomial_def_N + j][0] = (Polynomial_Mult_Result1[j] + Polynomial_Mult_Result2[j]) % FHEW_Mod_q;
			FHEW_ct_result[i * Polynomial_def_N + j][1] = (Polynomial_Mult_Result3[j] + Polynomial_Mult_Result4[j]) % FHEW_Mod_q;
		}
		memset(Polynomial_Mult_Result1, 0, sizeof(Polynomial_Mult_Result1));
		memset(Polynomial_Mult_Result2, 0, sizeof(Polynomial_Mult_Result2));
		memset(Polynomial_Mult_Result3, 0, sizeof(Polynomial_Mult_Result3));
		memset(Polynomial_Mult_Result4, 0, sizeof(Polynomial_Mult_Result4));
	}
}

void FHEW_SimulCom_1_1_Poly(int(*FHEW_ct1)[FHEW_columns], int(*FHEW_ct2)[FHEW_columns], int plainText, int *FHEW_SimulCom_result, int(*FHEW_ct3)[FHEW_columns])
{
	int Polynomial_Mult_Result1[Polynomial_def_N] = { 0 };
	int Polynomial_Mult_Result2[Polynomial_def_N] = { 0 };
	int temp = 0;

	Polynomial_Mult(&FHEW_ct1[0], 0, &FHEW_calced_Pk_b0_b1[0 * Polynomial_def_N], Polynomial_Mult_Result1);//D0,0 * b0
	Polynomial_Mult(&FHEW_ct1[0], 1, &FHEW_calced_Pk_b0_b1[1 * Polynomial_def_N], Polynomial_Mult_Result2);//D0,1 * b1
	for (int i = 0; i < Polynomial_def_N; ++i)
	{
		FHEW_SimulCom_result[i] = FHEW_new_Pk_a0[i];
	}

	for (int i = 0; i < plainText; i++)
	{
		temp = Polynomial_def_N - FHEW_SimulCom_result[Polynomial_def_N - 1];
		for (int j = Polynomial_def_N - 1; j > 0; j--)
		{
			FHEW_SimulCom_result[j] = FHEW_SimulCom_result[j - 1];
		}
		FHEW_SimulCom_result[0] = temp;
	}

	for (int l = 0; l < Polynomial_def_N; ++l)
	{
		FHEW_SimulCom_result[l] = (FHEW_SimulCom_result[l] + Polynomial_Mult_Result1[l] + Polynomial_Mult_Result2[l]) % FHEW_Mod_q;
		FHEW_calced_Pk_b0_b1[l] = FHEW_SimulCom_result[l];
		FHEW_calced_Pk_b0_b1[Polynomial_def_N + l] = FHEW_ct3[1 * Polynomial_def_N + l][0];
	}
}

void test()
{
	int FHEW_sk[Polynomial_def_N] = { 0 };
	int FHEW_Ciphertext1[Polynomial_def_N * FHEW_rows][FHEW_columns] = { 0 };
	int FHEW_Ciphertext2[Polynomial_def_N * FHEW_rows][FHEW_columns] = { 0 };
	int FHEW_Ciphertext3[Polynomial_def_N * FHEW_rows][FHEW_columns] = { 0 };
	int FHEW_SimulCom_Result[Polynomial_def_N] = { 0 };

	FHEW_Key_gen(FHEW_sk);
	FHEW_Encrypt(32, FHEW_sk, true, false, FHEW_Ciphertext1, 96);
	FHEW_Encrypt(17, FHEW_sk, false, true, FHEW_Ciphertext2);

	/*fstream fs;
	fs.open("d:\\1.csv", ios_base::out);
	fs << "a£º";
	for (int i = 0; i < 128; i++)
	{
		fs << FHEW_new_Pk_a0[i] << ',';
		if ((i+1) %50 == 0)
		{
			fs << '\n';
		}
	}
	fs << '\n';
	fs << '\n';
	fs << "b£º";
	for (int i = 0; i < 128; i++)
	{
		fs << FHEW_calced_Pk_b0_b1[i] << ',';
		if ((i + 1) % 50 == 0)
		{
			fs << '\n';
		}
	}*/
	//clock_t start, end;
	//start = clock();
	DWORD start, end;
	start = GetTickCount();
	FHEW_Ciphertext_Mult(FHEW_Ciphertext1, FHEW_Ciphertext2, FHEW_Ciphertext3);
	cout << GetTickCount() - start << endl;
	//cout << clock() - start << endl;

	/*fs << '\n';
	fs << '\n';
	fs << "ct result£º";
	for (int i = 0; i < 128; i++)
	{
		fs << FHEW_Ciphertext3[i][0] << ',';
		if ((i + 1) % 50 == 0)
		{
			fs << '\n';
		}
	}*/
	start = GetTickCount();
	FHEW_SimulCom_1_1_Poly(FHEW_Ciphertext1, FHEW_Ciphertext2, 17, FHEW_SimulCom_Result, FHEW_Ciphertext3);
	cout << GetTickCount() - start << endl;
	/*fs << '\n';
	fs << '\n';
	fs << "SimulCom£º";
	for (int i = 32; i < 96; i++)
	{
		fs << FHEW_SimulCom_Result[i] << ',';
	}
	fs.close();*/

	//Part decryption
	int compare_len = 64;
	int add_result = -1;
	for (int i = 0; i < compare_len; i++)
	{
		if (FHEW_SimulCom_Result[i + 32] != FHEW_Ciphertext3[i + 32][0])
		{
			FHEW_Ciphertext3[i + 32][0] = FHEW_SimulCom_Result[i + 32];
			if (i < compare_len / 2 - 1)
			{
				add_result = 0;
			}
			else
			{
				add_result = 1;
			}
		}
	}

	int FHEW_Ciphertext4[Polynomial_def_N * FHEW_rows][FHEW_columns] = { 0 };
	int FHEW_Ciphertext5[Polynomial_def_N * FHEW_rows][FHEW_columns] = { 0 };
	FHEW_Encrypt(32, FHEW_sk, true, false, FHEW_Ciphertext4, 96);
	start = GetTickCount();
	FHEW_Ciphertext_Mult(FHEW_Ciphertext4, FHEW_Ciphertext3, FHEW_Ciphertext5);
	cout << GetTickCount() - start << endl;
	start = GetTickCount();
	FHEW_SimulCom_1_1_Poly(FHEW_Ciphertext4, FHEW_Ciphertext3, 113, FHEW_SimulCom_Result, FHEW_Ciphertext5);
	cout << GetTickCount() - start << endl;


	int b = 0;
}