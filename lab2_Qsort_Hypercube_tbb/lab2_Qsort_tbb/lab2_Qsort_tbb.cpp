// lab2_Qsort_tbb.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include <omp.h>
#include <ctime>
#include <math.h>
#include "tbb\task_scheduler_init.h"
#include "tbb\parallel_for.h"
#include "tbb\blocked_range.h"

using namespace tbb;

// Проверка массива на корректность
bool PRKK(double *a, int n)
{
	for (int i = 0; i < n - 1; i++)
		if (a[i] > a[i + 1]) return false;

	return true;
}

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
				std::swap(pData[i], pData[PivotPos + 1]);
			PivotPos++;
		}
	}
	std::swap(pData[first], pData[PivotPos]);
	//pData[PivotPos] = tmp;

	SerialQuickSort(pData, first, PivotPos - 1);
	SerialQuickSort(pData, PivotPos + 1, last);
}

// Функция возводящая двойку в степень "step"
int Step2(int step)
{
	int res = 1;
	for (int i = 0; i < step; i++)
		res *= 2;
	return res;
}

// Возвращает номер блока из которого возьмется ведущий элемент
int SelectMedium(int _ProcNum, int _i, int _j)
{
	//(_j - _j % ( _ProcNum / 2^(_i - 1) ));

	if (_i == 0)
		return 0;
	else
	{
		return (_j / (_ProcNum / Step2(_i)))*(_ProcNum / Step2(_i - 1));
		//return (_j - _j % (_ProcNum / Step2(_i - 1)));
	}
}

// Костыльная версия с заранее заданными ведущими элементами
double KostelVersion(int _ProcNum, int _i, int _j, int _size)
{
	if (_i == 0)
		return 0;
	else
	{
		int res = -1;
		for (int k = 1; k < _i; ++k)
		{
			res -= Step2(k);
		}

		res += 2 * (_j / (_ProcNum / Step2(_i)));

		return ((res) * _size) / (8 * Step2(_i));
	}
}

// Функция для выбора пар смежных элементов
int SelectJ(int _ProcNum, int _i, int _j)
{
	if (_i == 0)
		return _j;
	else
	{
		return _j + (_j / (_ProcNum / Step2(_i))) * (_ProcNum / Step2(_i));
	}
}

// Функция смещения в части массива размером "c_size" с позиции "position" на "shift" позиций вправо/влево
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

class Parallel_q_Sotr
{
	double **Bloks;
	const int size;
	const int ProcNum;
	int i;
	int NumRange;

public:
	Parallel_q_Sotr(double **_Bloks, int _size, int _ProcNum, int _i, int _NumRange) : Bloks(_Bloks), size(_size), ProcNum(_ProcNum), i(_i), NumRange(_NumRange)
	{}

	void operator() (const blocked_range<int>& r) const
	{
		int begin = r.begin(), end = r.end();

		for (int jj = begin; jj < end; jj++) // Количество итераций выполняемых за раз
		{
			double* temp = new double[size + 1]; // Создание tmp массива

			//temp[0] = Bloks[SelectMedium(ProcNum, i, jj)][1]; // ВЫБОР ведущего
			temp[0] = KostelVersion(ProcNum, i, jj, NumRange*4);
			//std::cout << "Pivot evement: " << KostelVersion(ProcNum, i, jj, NumRange*4) << std::endl;
			int j = SelectJ(ProcNum, i, jj); // Определяет какое должно быть j(new) взависимости от i и j


											 // Заполнение массива tmp
											 // Штука для сохранения ведущего элемента
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
			// "l" - в итоге, размер массива tmp

			// Перекидывание элементов в зависимости от ведущего
			int PivotPos = 0;
			double Pivot = temp[0];
			for (int t = 1; t < l; t++)
			{
				if (temp[t] <= Pivot)
				{
					if (t != (PivotPos + 1))
						std::swap(temp[t], temp[PivotPos + 1]);
					PivotPos++;
				}
			}
			//swap(temp[0], temp[PivotPos]);

			int h = 1;
			// Заносим первый массив
			for (int k = 0; k < PivotPos; k++)
			{
				Bloks[j][k + 1] = temp[h];
				h++;
			}
			Bloks[j][0] = PivotPos;
			// Заносим второй массив
			for (int k = 0; k < l - (PivotPos + 1); k++)
			{
				Bloks[j + Step2(logCalc(ProcNum) - i)][k + 1] = temp[h];
				h++;
			}
			Bloks[j + Step2(logCalc(ProcNum) - i)][0] = l - (PivotPos + 1);
		}
	}
};

class endOfSort
{
	double **Bloks;
	double *pData;

public:
	endOfSort(double **_Bloks, double *_pData) : Bloks(_Bloks), pData(_pData)
	{}

	void operator() (const blocked_range<int>& r) const
	{
		int begin = r.begin(), end = r.end();

		for (int i1 = begin; i1 != end; i1++)
		{
			SerialQuickSort(Bloks[i1], 1, Bloks[i1][0]); // Сортировка элементов блоков
			int h = 0;
			for (int count = 0; count < i1; count++)
				h += Bloks[count][0];
			std::cout << i1 << " tread is ready. It have " << Bloks[i1][0] << " elements" << std::endl;
			for (int j1 = 0; j1 < Bloks[i1][0]; j1++) // Занесение блоков в исходный массив
			{
				pData[h] = Bloks[i1][j1 + 1];
				h++;
			}
		}
	}
};

void ParallelQuickSotr(double* pData, int size, int ProcNum, double &time, int NumRange)
{
	int part = size / (ProcNum * 2); // Размер блока

									 // Определение массива блоко. Пускай на 0-ом месте хранится размер блока
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
	double start_time = 0;
	double end_time = 0;
	
	start_time = clock();

	for (int i = 0; i <= logCalc(ProcNum); i++) // Количество итераций между блоками
	{
		task_scheduler_init init(ProcNum);  //omp_set_num_threads(ProcNum); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		Parallel_q_Sotr s(Bloks, size, ProcNum, i, NumRange);
		parallel_for(blocked_range<int>(0, ProcNum, 1), s);
	}

	//int h = 0;

	endOfSort v(Bloks, pData);
	parallel_for(blocked_range<int>(0, ProcNum * 2, 2), v);
	
	//for (int i1 = 0; i1 != ProcNum * 2; i1++)
	//{
	//	//SerialQuickSort(Bloks[i1], 1, Bloks[i1][0]);

	//	for (int j1 = 0; j1 < Bloks[i1][0]; j1++) // Занесение блоков в исходный массив
	//	{
	//		pData[h] = Bloks[i1][j1 + 1];
	//		h++;
	//	}
	//}

	end_time = clock();

	time = end_time - start_time;
}

int main()
{
	std::cout << "--My QuickSort--" << std::endl;

	// Начальные условия
	//------------------------------------
	int n = 10000000;
	int ProcNum = 2;
	int NumRange = 30000;
	//------------------------------------
	double *mas = new double[n];
	double *mascopy = new double[n];
	//------------------------------------
	for (int i = 0; i < n; i++) // Задание массива
	{
		mas[i] = rand() % NumRange - NumRange/2;
		mascopy[i] = mas[i];
		//if ((i < 50) || (i > n - 50)) std::cout << mas[i] << " "; //!
	}
	std::cout << std::endl << std::endl;

	//------------------------------------
	double time;

	//-----------------------------------Вызов QSort-----------------------------------	
	//---<TBB>---
	ParallelQuickSotr(mas, n, ProcNum, time, NumRange);


	/*iter_sort &a = *new(task::allocate_root()) iter_sort(mas, 0, n); // Создание корневой задачи
	task::spawn_root_and_wait(a); // Выполняет запуск корневой задачи*/
	//---</TBB>---
	//---------------------------------------------------------------------------------
	//SerialQuickSort(mas, 0, n);
	//---------------------------------------------------------------------------------


	std::cout << std::endl << std::endl << "Time Parallel : " << time / 1000 << " sec" << std::endl << std::endl;

	//--------------------------------------------------------------------------------
	if (n > 100) {
		for (int i = 0; i < 100; i++) // Вывод первых 100 элементов и последний 100 элементов
		{
			std::cout << mas[i] << " ";
		}

		std::cout << std::endl << std::endl;

		for (int i = n - 100; i < n; i++)
		{
			std::cout << mas[i] << " ";
		}
	}
	else
		for (int i = 0; i < n; i++)
		{
			std::cout << mas[i] << " ";
		}

	std::cout << std::endl << std::endl;

	// Проверка на корректность
	if (PRKK(mas, n)) std::cout << "Correctly";
	else std::cout << "Not correctly";

	std::cout << std::endl << std::endl;
	//--------------------------------------------------------------------------------

	double start_time = clock();
	SerialQuickSort(mascopy, 0, n);
	//ParallelQuickSotr(mascopy, n, ProcNum);
	double end_time = clock();

	std::cout << std::endl << std::endl << "Time Serial : " << abs((end_time - start_time) / 1000) << " sec" << std::endl << std::endl;

	if (PRKK(mascopy, n)) std::cout << "Correctly" << std::endl << std::endl;
	else std::cout << "Not correctly" << std::endl << std::endl;

	/*for (int i = 0; i < n; i++)
	{
	cout << mascopy[i] << " ";
	}*/

	std::cout << std::endl << std::endl;

	return 0;
}