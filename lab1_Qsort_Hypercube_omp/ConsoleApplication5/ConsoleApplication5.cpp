// ConsoleApplication5.cpp: ���������� ����� ����� ��� ����������� ����������.
//

// ���� ��������, �� ������� ���������
// SelectJ - ����� �� ��� ��� ���� �������� � ����� ��������� �������� ������� �������
// ������� ��� ��� ���������� � ���� �� ���� ����� tbb, ������� � ����� "2-� ������������ ������". � ��� ��� ���� ���������.
// ��� ��� ��� ���������� �� �����, ����� ������:)

#include "stdafx.h"
#include <iostream>
#include <omp.h>
#include <ctime>
#include <math.h>

using namespace std;

// �������� ������� �� ������������
bool PRKK(double *a, int n)
{
	for (int i = 0; i < n - 1; i++)
		if (a[i] > a[i + 1]) return false;

	return true;
}

//void swap(double &a, double &b)
//{
//	double c = a;
//	a = b;
//	b = c;
//}

void SerialQuickSort(double* pData, int first, int last) {
	if (first >= last)
		return;

	int PivotPos = first;
	double Pivot = pData[first];
	//Pivot = (pData[first] + pData[last]) / 2; //!
	//int tmp = pData[first];
	//pData[first] = Pivot;

	for (int i = first + 1; i <= last; i++) {
		if (pData[i] < Pivot) {
			if (i != PivotPos + 1)
				swap(pData[i], pData[PivotPos + 1]);
			PivotPos++;
		}
	}
	swap(pData[first], pData[PivotPos]);
	//pData[PivotPos] = tmp;

	SerialQuickSort(pData, first, PivotPos - 1);
	SerialQuickSort(pData, PivotPos + 1, last);
}

// ������� ���������� ������ � ������� "step"
int Step2(int step)
{
	int res = 1;
	for (int i = 0; i < step; i++)
		res *= 2;
	return res;
}

// ���������� ����� ����� �� �������� ��������� ������� �������
int SelectMedium(int _ProcNum, int _i, int _j)
{
	//(_j - _j % ( _ProcNum / 2^(_i - 1) ));

	if (_i == 0)
		return 0;
	else
	{
		return (_j - _j % (_ProcNum / Step2(_i - 1)));
	}
}
// ������� ��� ������ ��� ������� ���������
int SelectJ(int _ProcNum, int _i, int _j)
{
	if (_i == 0)
		return _j;
	else
	{
		return _j + (_j / (_ProcNum / Step2(_i))) * (_ProcNum / Step2(_i));
	}
}

// ������� �������� � ����� ������� �������� "c_size" � ������� "position" �� "shift" ������� ������/�����
void ShiftMas(double* _pData, int position, int shift, int c_size)
{
	if (shift >= 0)
		for (int i = c_size + position - 1; i >= position; i--)
		{
			_pData[i + shift] = _pData[i];
		}
	else
		for (int i = position; i <= c_size + position - 1; i++)
		{
			_pData[i + shift] = _pData[i];
		}
}

int logCalc(int twoT)
{
	int res = 0;
	while (twoT != 1)
	{
		twoT = twoT / 2;
		res++;
	}

	return res;
}

void ParallelQuickSotr(double* pData, int size, int ProcNum, int ThreadNum)
{
	int part = size / (ProcNum * 2); // ������ �����

									 // ����������� ������� �����. ������ �� 0-�� ����� �������� ������ �����
	double** Bloks = new double*[ProcNum * 2];
	for (int k = 0; k < ProcNum * 2; k++)
		Bloks[k] = new double[size + 1];
	//-------------------------------------------
	for (int i = 0; i < ProcNum * 2 - 1; i++)
	{
		Bloks[i][0] = part;
		for (int j = 1; j < part + 1; j++)
		{
			Bloks[i][j] = pData[i*part + j - 1];
		}
	}

	int r = 1;
	for (int j = (ProcNum * 2 - 1) * part; j < size; j++)
	{
		Bloks[ProcNum * 2 - 1][r] = pData[j];
		r++;
	}
	Bloks[ProcNum * 2 - 1][0] = r - 1;

	/*for (int i1 = 0; i1 < ProcNum * 2; i1++)
	{
	cout << i1 << " block: ";
	for (int j1 = 0; j1 <= Bloks[i1][0]; j1++)
	{
	cout << Bloks[i1][j1] << " ";
	}
	cout << endl;
	}*/
	//-------------------------------------------

	for (int i = 0; i <= logCalc(ProcNum); i++) // ���������� �������� ����� �������
	{
		/*int countLife = ProcNum / Step2(i); // ��� ����������� ��������
		double pivot2;*/

		omp_set_num_threads(ProcNum);
		//#pragma omp parallel num_threads(ProcNum)
#pragma omp parallel for //shared(size, ProcNum, Bloks, i)
		for (int jj = 0; jj < ProcNum; jj++) // ���������� �������� ����������� �� ���
		{
			double* temp = new double[size + 1]; // �������� tmp �������

			temp[0] = Bloks[SelectMedium(ProcNum, i, jj)][1]; // ����� ��������
															  /*cout << endl << "Pivot[" << jj << "] : " << temp[0] << endl;
															  cout << endl;*/

			int j = SelectJ(ProcNum, i, jj); // ���������� ����� ������ ���� j(new) ������������ �� i � j


											 // ���������� ������� tmp
											 // ����� ��� ���������� �������� ��������
											 /*if ((countLife < ProcNum / (i + 1)) && (countLife > 0))
											 {
											 temp[0] = pivot2;
											 countLife--;
											 }
											 else
											 {*/
											 /*pivot2 = temp[0];
											 countLife = ProcNum / Step2(i) - 1;
											 }*/

			int l = 1;
			for (int k = 1; k < Bloks[j][0] + 1; k++)
			{
				temp[l] = Bloks[j][k];
				l++;
			}
			for (int k = 1; k < Bloks[j + Step2(logCalc(ProcNum) - i)][0] + 1; k++)
			{
				temp[l] = Bloks[j + Step2(logCalc(ProcNum) - i)][k];
				l++;
			}
			// "l" - � �����, ������ ������� tmp

			// ������������� ��������� � ����������� �� ��������
			int PivotPos = 0;
			double Pivot = temp[0];
			for (int t = 1; t < l; t++)
			{
				if (temp[t] <= Pivot)
				{
					if (t != (PivotPos + 1))
						swap(temp[t], temp[PivotPos + 1]);
					PivotPos++;
				}
			}
			//swap(temp[0], temp[PivotPos]);

			int h = 1;
			// ������� ������ ������
			for (int k = 0; k < PivotPos; k++)
			{
				Bloks[j][k + 1] = temp[h];
				h++;
			}
			Bloks[j][0] = PivotPos;
			// ������� ������ ������
			for (int k = 0; k < l - (PivotPos + 1); k++)
			{
				Bloks[j + Step2(logCalc(ProcNum) - i)][k + 1] = temp[h];
				h++;
			}
			Bloks[j + Step2(logCalc(ProcNum) - i)][0] = l - (PivotPos + 1);
		}

		// �������� ������
		/*for (int i1 = 0; i1 < ProcNum * 2; i1++)
		{
		cout << i1 << " block: ";
		for (int j1 = 0; j1 <= Bloks[i1][0]; j1++)
		{
		cout << Bloks[i1][j1] << " ";
		}
		cout << endl;
		}
		cout << endl;*/
	}

#pragma omp parallel for
	for (int i1 = 0; i1 < ProcNum * 2; i1++)
	{
		SerialQuickSort(Bloks[i1], 1, Bloks[i1][0]); // ���������� ��������� ������

		int h = 0;

		for (int count = 0; count < i1; count++)
		{
			h += Bloks[i1][count];
		}

		for (int j1 = 0; j1 < Bloks[i1][0]; j1++) // ��������� ������ � �������� ������
		{
			pData[h] = Bloks[i1][j1 + 1];
			h++;
		}
	}


	// ����� ������
	/*for (int i1 = 0; i1 < ProcNum * 2; i1++)
	{
	cout << i1 << " block: ";
	for (int j1 = 0; j1 <= Bloks[i1][0]; j1++)
	{
	cout << Bloks[i1][j1] << " ";
	}
	cout << endl;
	}*/
}

int main()
{
	cout << "--My QuickSort--" << endl;

	// ��������� �������
	//------------------------------------
	int n = 1000000;
	int ProcNum = 2;
	int ThreadNum = Step2(ProcNum);
	//------------------------------------
	double *mas = new double[n];
	double *mascopy = new double[n];
	//------------------------------------
	//mas[0] = 0; mas[1] = 2; mas[2] = 5; mas[3] = -6; mas[4] = -3; mas[5] = 1; mas[6] = 8; mas[7] = -5;
	for (int i = 0; i < n; i++) // ������� �������
	{
		mas[i] = rand() % n - n/2;
		mascopy[i] = mas[i];
		//cout << mas[i] << " "; //!
	}
	//------------------------------------
	double start_time = clock();

	//-----------------------------------����� QSort-----------------------------------	
	SerialQuickSort(mas, 0, n);
	//--------------------------------------------------------------------------------

	double end_time = clock();

	cout << endl << endl << "Time: " << (end_time - start_time) / 1000 << " sec" << endl << endl;

	//--------------------------------------------------------------------------------
	if (n > 100) {
		for (int i = 0; i < 100; i++) // ����� ������ 100 ��������� � ��������� 100 ���������
		{
			cout << mas[i] << " ";
		}

		cout << endl << endl;

		for (int i = n - 100; i < n; i++)
		{
			cout << mas[i] << " ";
		}
	}
	else
		for (int i = 0; i < n; i++)
		{
			cout << mas[i] << " ";
		}

	cout << endl << endl;

	// �������� �� ������������
	if (PRKK(mas, n)) cout << "Correctly";
	else cout << "Not correctly";

	cout << endl << endl;
	//--------------------------------------------------------------------------------

	start_time = clock();
	ParallelQuickSotr(mascopy, n, ProcNum, ThreadNum);
	end_time = clock();

	cout << endl << endl << "Time: " << abs((end_time - start_time) / 1000) << " sec" << endl << endl;

	if (PRKK(mascopy, n)) cout << "Correctly" << endl << endl;
	else cout << "Not correctly" << endl << endl;

	/*for (int i = 0; i < n; i++)
	{
	cout << mascopy[i] << " ";
	}*/

	cout << endl << endl;

	return 0;
}
