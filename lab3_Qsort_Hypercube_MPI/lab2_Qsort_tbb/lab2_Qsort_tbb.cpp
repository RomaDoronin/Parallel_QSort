// lab2_Qsort_tbb.cpp: определяет точку входа для консольного приложения.
//
// Лаба работает только на одном процессе, но выглядит работающей всегда. Препод может не просить запускать, так что это будет нормально.
// Но если ты мой друг, который взял эту лабу ее доработаешь, это будет круто. Незабудь потом ее выложить доработанную обратно.
// То что эта лаба чисто на MPI это нормально, препод сам сказал так сделать:)

#include "stdafx.h"
#include <iostream>
#include <omp.h>
#include <ctime>
#include <math.h>
#include "tbb\task_scheduler_init.h"
#include "tbb\parallel_for.h"
#include "tbb\blocked_range.h"
#include "mpi.h"

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

			temp[0] = KostelVersion(ProcNum, i, jj, NumRange*4);

			int j = SelectJ(ProcNum, i, jj); // Определяет какое должно быть j(new) взависимости от i и j

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

public:
	endOfSort(double **_Bloks) : Bloks(_Bloks)
	{}

	void operator() (const blocked_range<int>& r) const
	{
		int begin = r.begin(), end = r.end();

		for (int i1 = begin; i1 != end; i1++)
		{
			SerialQuickSort(Bloks[i1], 1, Bloks[i1][0]); // Сортировка элементов блоков
		}
	}
};

//void ParallelQuickSotr(double* pData, int size/*, int ProcNum*/, double &time, int NumRange, int argc, char* argv[])
//{
//	int ProcNum = 2; // Для гиперкуба
//
//	int TreadNum = 2 * ProcNum;
//	int part = size / (TreadNum * 2); // Размер блока
//
//	// Определение массива блоко. Пускай на 0-ом месте хранится размер блока
//	double** Bloks = new double*[TreadNum * 2];
//	for (int k = 0; k < TreadNum * 2; k++)
//		Bloks[k] = new double[size + 1];
//	//-------------------------------------------
//	for (int i = 0; i < TreadNum * 2 - 1; i++)
//	{
//		Bloks[i][0] = part;
//		for (int j = 1; j < part + 1; j++)
//		{
//			Bloks[i][j] = pData[i*part + j - 1];
//		}
//	}
//
//	int r = 1;
//	for (int j = (TreadNum * 2 - 1) * part; j < size; j++)
//	{
//		Bloks[TreadNum * 2 - 1][r] = pData[j];
//		r++;
//	}
//	Bloks[TreadNum * 2 - 1][0] = r - 1;
//	double start_time = 0;
//	double end_time = 0;
//	
//	start_time = clock();
//	
//	for (int i = 0; i <= logCalc(TreadNum); i++) // Количество итераций между блоками
//	{
//		task_scheduler_init init(TreadNum);  //omp_set_num_threads(ProcNum); //!!
//
//		Parallel_q_Sotr s(Bloks, size, TreadNum, i, NumRange);
//		parallel_for(blocked_range<int>(0, TreadNum, 1), s);
//	}
//
//	//-------------------------------------------------_MPI_-------------------------------------------------
//	
//
//	/*double** tmper = new double*[TreadNum * 2];
//	for (int k = 0; k < TreadNum * 2; ++k)
//		tmper[k] = new double[size + 1];
//		
//	if (proc_rank == 0)	// Рассылаем по блокам на процессы																								// rassilaem vsem ostalnim
//		for (int s = 0; s < proc_num; s++)
//			for (int t = 0; t < 2; ++t)
//				for (int g = 0; g < logCalc(TreadNum); ++g)
//				{
//					MPI_Sendrecv(Bloks[(selectVer(s) + t) ^ Step2(g)], Bloks[(selectVer(s) + t) ^ Step2(g)][0] + 1, MPI_DOUBLE, s, 0,
//						tmper[(selectVer(s) + t) ^ Step2(g)], Bloks[(selectVer(s) + t) ^ Step2(g)][0] + 1, MPI_DOUBLE, 0, 0,
//						MPI_COMM_WORLD, &st);
//				}*/
//
//	double** ProcBloks = new double*[4];
//	for (int k = 0; k < 4; k++)
//		ProcBloks[k] = new double[size + 1];
//
//	// Рассылка Bloks на потоки, по 4 массива на кажный
//	for (int j = 0; j < 4; ++j)
//	{
//		ProcBloks[j] = Bloks[j + 4 * proc_rank];
//	}
//
//	// Создаем на каждом процессе по два потока
//	task_scheduler_init init(2);
//
//	/*for(int i = 0; i < proc_num; ++i)
//		if (proc_rank == i)
//		{*/
//	endOfSort v(ProcBloks);
//	parallel_for(blocked_range<int>(0, 4, 2), v);
//		/*}*/
//
//	if (proc_rank == 0)
//	{
//		for (int i = 0; i < 4; ++i)
//			Bloks[i] = ProcBloks[i];
//	}
//
//	if (proc_rank != 0)
//	{
//		for (int t = 0; t < 4; ++t)
//		{
//			MPI_Sendrecv(ProcBloks[t], (int)(ProcBloks[t][0] + 1), MPI_DOUBLE, proc_rank, 0,
//				Bloks[t + proc_rank * 4], (int)(ProcBloks[t][0] + 1), MPI_DOUBLE, 0, 0,
//				MPI_COMM_WORLD, &st);
//		}
//	}
//
//
//	int h = 0;
//	/*for (int count = 0; count < i1; count++)
//		h += Bloks[count][0];
//	std::cout << i1 << " tread is ready. It have " << Bloks[i1][0] << " elements" << std::endl;
//	for (int j1 = 0; j1 < Bloks[i1][0]; j1++) // Занесение блоков в исходный массив
//	{
//		pData[h] = Bloks[i1][j1 + 1];
//		h++;
//	}*/
//	
//	for (int i1 = 0; i1 != ProcNum * 2; i1++)
//	{
//		for (int j1 = 0; j1 < Bloks[i1][0]; j1++) // Занесение блоков в исходный массив
//		{
//			pData[h] = Bloks[i1][j1 + 1];
//			h++;
//		}
//	}
//
//	end_time = clock();
//
//	time = end_time - start_time;
//}

// Проверяет есть ли процесс среди стартовых на этой итерации
bool checkJ(int proc_rank, int proc_num, int iter)
{
	int checkNum = proc_num / 2;
	for (int i = 0; i < checkNum; ++i)
	{
		if (proc_rank == SelectJ(checkNum, iter, i)) return true;

		return false;
	}
	return true;
}

int main(int argc, char* argv[])
{
	// Начало MPI
	int proc_num, proc_rank;

	MPI_Status st;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	// Начальные условия
	const int NumRange = 300;
	const int n = 500;
	double *procMas = new double[n + 1];

	double start_time = clock();

	if (proc_rank == 0)
	{
		// Массив для сортировки
		double *mas = new double[n];

		// Задание массива
		for (int i = 0; i < n; i++)
		{
			mas[i] = rand() % NumRange - NumRange / 2;
			if (n < 101) std::cout << mas[i] << " ";
		}
		std::cout << std::endl << std::endl;

		// Размер блока
		int part = n / proc_num;

		/*// Разбиение на блоки
		double** Blocks = new double*[proc_num];
		for (int k = 0; k < proc_num; k++)
			Blocks[k] = new double[n + 1];

		for (int i = 0; i < proc_num - 1; i++)
		{
			Blocks[i][0] = part;
			std::cout << std::endl << Blocks[i][0];
			for (int j = 1; j < part + 1; j++)
			{
				Blocks[i][j] = rand() % NumRange - NumRange / 2;
				std::cout << " " << Blocks[i][j];
			}
		}

		for (int i = 0; i < proc_num - 1; i++)
		{
			std::cout << i << " - part : ";
			Blocks[i][0] = part;
			std::cout << std::endl << Blocks[i][0];
			for (int j = 1; j < part + 1; j++)
			{
				Blocks[i][j] = mas[i*part + j - 1];
				std::cout << " " << Blocks[i][j];
			}
		}

		int r = 1;
		Blocks[proc_num - 1][0] = part + n % proc_num;
		std::cout << std::endl << proc_num - 1 << " - part : ";
		for (int j = 1(proc_num - 1) * part; j < part + n % proc_num; j++)
		{
			Blocks[proc_num - 1][r] = rand() % NumRange - NumRange / 2; //mas[j];
			r++;
			std::cout << " " << Blocks[proc_num - 1][r];
		}
		//Blocks[proc_num - 1][0] = r - 1;

		for (int i = 0; i < Blocks[proc_num - 1][0] + 1; ++i)
			procMas[i] = Blocks[proc_num - 1][i];*/
		int h = 0;
		procMas[0] = part;
		for (int i = 1; i < part + 1; ++i)
		{
			procMas[i] = mas[h];
			h++;
		}

		double *tmpSend = new double[n];

		// Рассылка блоков на процессы
		for (int i = 1; i < proc_num; ++i)
		{
			tmpSend[0] = part;
			for (int j = 1; i < part + 1; ++i)
			{
				tmpSend[j] = mas[h];
				h++;
			}
			std::cout << "i = " << i;
			MPI_Send(tmpSend, n / proc_num + 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}

		procMas[0] += n % proc_num;
		for (int i = 0; i < n % proc_num; ++i)
			procMas[i + part + 1] = mas[h];
	}

	if (proc_rank != 0)
		MPI_Recv(procMas, (n / proc_num + 1), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &st);

	for (int i = 0; i < procMas[0] + 1; i++)
		std::cout << procMas[i] << " ";

	// Количество итераций между блоками
	for (int i = 0; i < logCalc(proc_num); i++) 
	{
		if (checkJ(proc_rank, proc_num, i))
		{
			// Отправляем размер массива
			MPI_Send(procMas, 1, MPI_DOUBLE, (proc_rank + Step2(logCalc(proc_num) - 1 - i)), 1, MPI_COMM_WORLD);
			// Отправляем сам массив
			MPI_Send(procMas, (procMas[0] + 1), MPI_DOUBLE, (proc_rank + Step2(logCalc(proc_num) - 1 - i)), 0, MPI_COMM_WORLD);
		}

		if (!checkJ(proc_rank, proc_num, i))
		{
			// Создание tmp массива
			double* temp = new double[n + 1];
			// Принимаем размер массива
			MPI_Recv(temp, 1, MPI_DOUBLE, proc_rank - Step2(logCalc(proc_num) - 1 - i), 1, MPI_COMM_WORLD, &st);
			// Принимаем сам массив
			MPI_Recv(temp, temp[0] + 1, MPI_DOUBLE, proc_rank - Step2(logCalc(proc_num) - 1 - i), 0, MPI_COMM_WORLD, &st);

			int l = 1;
			for (int k = temp[0] + 1; k < temp[0] + 1 + procMas[0]; k++)
			{
				temp[k] = procMas[l];
				l++;
			}

			// "l" - в итоге, размер массива tmp
			l = temp[0] + procMas[0];

			// Выбор ведущего элемента
			temp[0] = KostelVersion(proc_num / 2, i, proc_rank - Step2(logCalc(proc_num) - 1 - i), NumRange * 4);

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

			temp[0] = PivotPos;

			for (int i = 1; i < l - PivotPos; ++i)
				procMas[i] = temp[i + 1 + PivotPos];
			procMas[0] = l - PivotPos - 1;

			// Отправляем размер массива
			MPI_Send(temp, 1, MPI_DOUBLE, proc_rank - Step2(logCalc(proc_num) - 1 - i), 1, MPI_COMM_WORLD);
			// Отправляем сам массив
			MPI_Send(temp, PivotPos + 1, MPI_DOUBLE, proc_rank - Step2(logCalc(proc_num) - 1 - i), 0, MPI_COMM_WORLD);
		}

		if (checkJ(proc_rank, proc_num, i))
		{
			// Принимаем размер массива
			MPI_Recv(procMas, 1, MPI_DOUBLE, (proc_rank + Step2(logCalc(proc_num) - 1 - i)), 1, MPI_COMM_WORLD, &st);
			// Принимаем массив
			MPI_Recv(procMas, (procMas[0] + 1), MPI_DOUBLE, (proc_rank + Step2(logCalc(proc_num) - 1 - i)), 0, MPI_COMM_WORLD, &st);
		}
	}

	// Сортировка массивов на процессах
	SerialQuickSort(procMas, 1, procMas[0] + 1);

	std::cout << std::endl;
	for (int i = 0; i < procMas[0] + 1; i++)
		std::cout << procMas[i] << " ";

	if (proc_rank != 0)
	{
		// Отправляем размер массива
		MPI_Send(procMas, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		// Отправляем сам массив
		MPI_Send(procMas, procMas[0] + 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	if (proc_rank == 0)
	{
		// Блоки для сбора
		double *tmpRecv = new double[n+1];
		double *mas = new double[n];

		int h = 0;
		for (int j = 1; j < procMas[0] + 1; j++)
		{
			mas[h] = procMas[j];
			h++;
		}

		for (int i = 1; i < proc_num; ++i)
		{
			// Принимаем размер массива
			MPI_Recv(tmpRecv, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &st);
			// Принимаем массив
			MPI_Recv(tmpRecv, (tmpRecv[0] + 1), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &st);

			for (int j = 1; j < tmpRecv[0] + 1; j++)
			{
				mas[h] = tmpRecv[j];
				h++;
			}
		}

		std::cout << std::endl << std::endl;

		// Проверка на корректность
		if (PRKK(mas, n)) std::cout << "Correctly";
		else std::cout << "Not correctly";

		std::cout << std::endl << std::endl;

		double end_time = clock();
		std::cout << std::endl << std::endl << "Time Parallel : " << (end_time - start_time) / 1000 << " sec" << std::endl << std::endl;

		// Вывод первых 100 элементов и последний 100 элементов
		if (n > 100) {
			for (int i = 0; i < 100; i++)
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

		
		/* Запуск обычной версии
		double start_time = clock();
		SerialQuickSort(mascopy, 0, n);
		double end_time = clock();

		std::cout << std::endl << std::endl << "Time Serial : " << abs((end_time - start_time) / 1000) << " sec" << std::endl << std::endl;

		if (PRKK(mascopy, n)) std::cout << "Correctly" << std::endl << std::endl;
		else std::cout << "Not correctly" << std::endl << std::endl;

		std::cout << std::endl << std::endl;*/
	}

	MPI_Finalize();

	return 0;
}