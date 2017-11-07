---
layout: concise
title: Параллельное решение уравнения переноса
---

# Постановка задачи и разностная схема

Рассмотрим уравнение переноса с периодическими граничными условиями

$$
\begin{aligned}
&\frac{\partial u}{\partial t} + \frac{\partial u}{\partial x} = 0, \qquad x \in [0, 1]\\
&u\big|_{t=0} = u_0(x)\\
&u\big|_{x=0} = u\big|_{x=1}
\end{aligned}
$$

Для дискретизации этой задачи используем схему «явный левый уголок»

$$
\begin{aligned}
&\frac{u^{n+1}_m - u^n_m}{\tau} + \frac{u^n_m - u^n_{m-1}}{h} = 0, &\quad m = 1, 2, \dots, M - 1\\
&\frac{u^{n+1}_0 - u^n_0}{\tau} + \frac{u^n_0 - u^n_{M-1}}{h} = 0, &\quad m = 0\\
&u^0_m = u_0(x_m), &\qquad m = 0, \dots, M-1
\end{aligned}
$$

Запись разностной схемы можно упростить, если ввести число Куранта $$\sigma = \frac{\tau}{h}$$:

$$
\begin{aligned}
&u^{n+1}_m = (1 - \sigma) u^n_m + \sigma u^n_{m-1}, &\quad m = 1, 2, \dots, M-1\\
&u^{n+1}_0 = (1 - \sigma) u^n_0 + \sigma u^n_{M-1}, &\quad m = 0
\end{aligned}
$$

# Последовательная программа

Будем решать задачу, двигаясь шагами по времени. На каждом шаге необходимо
по известным знаениям $$u^n_m, m = 0, 1, \dots M-1$$ вычислить значения $$u^{n+1}_m, m = 0, 1, \dots, M-1$$.
Таким образом, для решения задачи достаточно двух массивов --- один для известных $$u^n_m$$, другой для
вычисляемых $$u^{n+1}_m$$. В конце шага по времени эти массивы меняются местами: старые значения $$u^n_m$$
уже не нужны, а значения $$u^{n+1}_m$$ становятся значениями $$u^{n}_m$$ для следующего шага по времени.

Шаг по времени подберем таким, что число Куранта было примерно равно 0.5.

```c++
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <fstream>

// Начальное условие - ступенька
double u0(double x) {
	if (x > 0.4 && x < 0.6)
		return 1;
	return 0;
}

// Сохраняем содержимое массива в файл в формате gnuplot
void save(const int M, const double h, const double *u, const std::string &filename) {
	std::ofstream f(filename);

	for (int m = 0; m < M; m++) {
		double x = m*h;
		f << x << " " << u[m] << "\n";
	}
}

// Число узлов M задается аргументом программы при запуске
int main(int argc, char **argv) {
	if (argc < 2) {
		std::cerr << "USAGE: ./main <M>" << std::endl;
		return 1;
	}

	const int M = atoi(argv[1]);
	const double h = 1.0 / M; // Шаг по пространству

	// Массив u хранит значение u^n на текущем слое по времени
	double *u     = new double[M];
	// Массив unext хранит значение u^{n+1}
	double *unext = new double[M];

	// Заполняем начальное условие
	for (int m = 0; m < M; m++) {
		double x = m*h;
		u[m] = u0(x);
	}

	const double tmax = 1;    // Конечный момент времени
	double sigma = 0.5;
	double dt = sigma * h;    // Этот шаг может не делить tmax нацело!
	const int N = tmax / dt;  // Округляем tmax/dt вниз до целого

	dt = tmax / N;
	sigma = dt / h;           // Корректируем dt и sigma


	// Делаем ровно N шагов
	for (int n = 0; n < N; n++) {
		for (int m = 1; m < M; m++)
			unext[m] = (1 - sigma) * u[m] + sigma * u[m-1];
		unext[0] = (1 - sigma) * u[0] + sigma * u[M-1];

		std::swap(u, unext); // Меняем местами массивы
	}

	save(M, h, u, "serial.csv");

	delete[] u;
	delete[] unext;

	return 0;
}
```
# Разбиение области

Разобьем работу по вычислению $$M$$ значений с $$n+1$$ слоя по времени между процессами. Условимся разделить данные «честно»,
то есть так, чтобы число элементов у разных процессов отличалось не более, чем на один.

Пусть процессы занумерованы от $$0$$ до $$P-1$$. Следующая формула задает требуемое разбиение
$$
M_i = \left\lfloor\frac{M}{P}\right\rfloor + \begin{cases}
1, &i < \operatorname{mod}(M, P)\\
0, &i \geqslant \operatorname{mod}(M, P)
\end{cases}
$$

Введем на каждом процессе переменные `mb` и `me` означающие начало и конец его блока соответственно. Добавим необходимый MPI код:

```c++
...

#include <mpi.h>  // К заголовочным файлам добавляем mpi.h

...

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv); // Включаем MPI

	if (argc < 2) {
		std::cerr << "USAGE: ./main <M>" << std::endl;
		return 1;
	}

	const int M = atoi(argv[1]);

	// Узнаем число процессов и номер каждого процесса
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Эти вычисления продублированы всеми процессами
	const int Mi[size];
	// Сначала раздаем всем поровну по M / size (округление вниз!)
	for (int i = 0; i < size; i++)
		Mi[i] = M / size;
	// Первым M % size процессам добавляем по одному элементу
	for (int i = 0; i < M % size; i++)
		Mi[i]++;

	// Каждый процесс вычисляет собственные mb и me
	int mb, me;
	mb = 0;
	for (int i = 0; i < rank; i++)
		mb += Mi[i];
	me = mb + Mi[rank] - 1;

	std::cout << "Процесс " << rank << " обрабатывает элементы с "
		<< mb << " до " << me << std::endl;

	// Вычисления сейчас закомментированы

	MPI_Finalize(); // Выключаем MPI

	return 0;
}
```

Компилируем этот код с помощью компилятора `mpicxx`, запускаем его на трех процессах и получаем вывод
```
Процесс 0 обрабатывает элементы с 0 до 66
Процесс 1 обрабатывает элементы с 67 до 133
Процесс 2 обрабатывает элементы с 134 до 199
```
