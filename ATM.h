#ifndef FHEW_H
#define FHEW_H
#include <iostream>
#include <cstdlib> 
#include <ctime> 
#include<fstream>
#include <Windows.h>
using namespace std;

const int Polynomial_def_N = 1024;
const int FHEW_rows = 6;
const int FHEW_columns = 2;
const int FHEW_Mod_q = 512;
const int vgprime[3] = { 1, 1 << 11, 1 << 22 };

extern int FHEW_new_Pk_a0[Polynomial_def_N];
extern int FHEW_calced_Pk_b0_b1[2 * Polynomial_def_N];


typedef struct 
{
	double std_dev;
	int max; // Size of table. Set to 0 to ignore table and use rejection sampling
	int offset;
	const float* table; // CDF of Gaussian of standard deviation std_dev centered around offset
} Distrib;

const float Chi1_Table[23] = 
{
	1.12011750313263e-14, 2.38717233762211e-12, 3.04966971020178e-10,
	2.34394541319773e-8, 1.08538196465647e-6, 0.0000303513786306856,
	0.000514575939439740, 0.00532464699317562, 0.0340111330223921,
	0.136723892128727, 0.357520614142345, 0.642479385857655,
	0.863276107871273, 0.965988866977608, 0.994675353006824,
	0.999485424060560, 0.999969648621369, 0.999998914618035,
	0.999999976560546, 0.999999999695033, 0.999999999997613,
	0.999999999999989, 1.00000000000000 
};

const Distrib Chi1 = 
{
	1.4,
	23,
	11,
	Chi1_Table
};

void FHEW_Key_gen(int *FHEW_sk);
int Sample(const Distrib& Chi);

void FHEW_Encrypt(int plain_text, int *FHEW_sk, bool Pk_a0_save_flag, bool Pk_b0_b1_save_flag, int(*FHEW_ct)[FHEW_columns], int plain_text2 = -1);
void Polynomial_Mult(int(*FHEW_ct)[FHEW_columns], int ct_columns, int *FHEW_ct2, int *FHEW_result_ct);
void FHEW_Ciphertext_Mult(int(*FHEW_ct1)[FHEW_columns], int(*FHEW_ct2)[FHEW_columns], int(*FHEW_ct_result)[FHEW_columns]);

void test();

void FHEW_SimulCom_1_1_Poly(int(*FHEW_ct1)[FHEW_columns], int(*FHEW_ct2)[FHEW_columns], int plainText, int *FHEW_SimulCom_result, int(*FHEW_ct3)[FHEW_columns]);



#endif FHEW_H