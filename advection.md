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

	for (int m = 0; m < M; i++) {
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

	// Массив u хранит значение u^n на текущем слое по времени
	double *u     = new double[M];
	// Массив unext хранит значение u^{n+1}
	double *unext = new double[M];

	// Заполняем начальное условие
	for (int m = 0; m < M; m++) {
		double x = m*h;
		u[m] = u0(x);
	}

	const double tmax = 1;   // Конечный момент времени
	double sigma = 0.5;
	double dt = sigma * h;   // Этот шаг может не делить tmax нацело!
	const int N = tmax / dt; // Округляем tmax/dt вниз до целого

	dt = tmax / N;
	sigma = dt / h;          // Корректируем dt и sigma


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
