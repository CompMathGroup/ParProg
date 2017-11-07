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

Решение можно посмотреть с помощью программы `gnuplot`. Для этого нужно её запустить и ввести команду
```
gnuplot> plot 'serial.csv'
```
Для выхода из программы можно нажать `Ctrl+D`.

# Разбиение области

Разобьем работу по вычислению $$M$$ значений с $$n+1$$ слоя по времени между процессами. Условимся разделить данные «честно»,
то есть так, чтобы число элементов у разных процессов отличалось не более, чем на один. При этом аналогичное разбиение применим
и к массиву, сожржащему $$n$$ слой по времени.

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
    int Mi[size];
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

# Действия, не требующие взаимодействия

Разбив область на зоны ответственности процессов, можно выполнить некоторые действия, не требующие обмена данными с
соседями. Например, можно заполнить начальные значения в свою часть массива `u`. Для простоты оставим массивам
`u` и `unext` тот же размер, что и в последовательной программе. Это слегка расточительно, так как все процессы будут
иметь свою копию массива, и большая часть массива просто не используется ими. Зато такая схема позволяет сравнительно
легко переписать последовательный код в параллельный код.

```c++
...

#include <mpi.h>  // К заголовочным файлам добавляем mpi.h

...

// Поправим save так, чтобы в файл записывалась только часть массива,
// хранимая на данном процессе
void save(const int mb, const int me, const double h, const double *u, const std::string &filename) {
    std::ofstream f(filename);

    for (int m = mb; m <= me; m++) {
        double x = m*h;
        f << x << " " << u[m] << "\n";
    }
}

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
    int Mi[size];
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

    double *u = new double[M];
    double *unext = new double[M];

    for (int m = mb; m <= me; m++) {
        double x = m*h;
        u[m] = u0(x);
    }

    // Вычисления сейчас закомментированы

    save(mb, me, h, u, "par." + std::to_string(rank) + ".csv");

    delete[] u;
    delete[] unext;

    MPI_Finalize(); // Выключаем MPI

    return 0;
}
```

Каждый процесс записывает свою часть результата в файл с именем `par.<rank>.csv`. Для этого используется
функция `to_string` из стандарта c++11, так что компилировать придется с указанием этого компилятору:
```bash
mpicxx -std=c++11 parallel.cpp -o parallel
```

# Обмен данными

Каждый процесс не может вычислить свою часть $$n+1$$ слоя по времени без взаимодействия с другими процессами, 
так как предыдущий слой по времени также разбит между процессами (каждый процесс «отвечает» за свою часть данных
с `mb` до `me`). Единственное значение, при вычислении которого возникает проблема --- это элемент
`mb` на $$n+1$$ слое по времени, так как для его вычисления требуется значение `mb-1` с $$n$$-го слоя, но это
значение хранится на другом процессе!

Для того, чтобы запросить данные у соседнего процесса, необходимо, чтобы процесс-отправитель вызвал функцию
`MPI_Send`, а процесс-приемник вызвал `MPI_Recv`. Учтем, что мы хотим написать универсальный код, который будет исполнен
всеми процессами. Фиксированный процесс с номером `rank` должен поучаствовать в двух обменах: получить от левого соседа
значение элемента `mb-1` и отправить правому соседу свое значение `me` (которое для правого соседа как-раз и будет его
`mb-1` элементом). Заметим, что в такой схеме периодическое граничное условие естественным образом получается,
если считать последний процесс левым соседом нулевого.

Функции `MPI_Send/MPI_Recv` ориентированы на передачу целых массивов, так что для передачи одного числа необходмо
использовать `&` (можно считать, что `&x` --- массив из одного элемента, содержащий значение `x`). При использовании
функций MPI можно сообщениям присваивать дополнительный идентификатор `tag`. В основном, он служит для проверки
корректности пересылки (принято именно то сообщение, которое ожидалось).

```c++
    // Вычисление номеров левого и правого соседа
    int left = (rank > 0) ? rank - 1 : size - 1;
    int right = (rank < size - 1) ? rank + 1 : 0;

    // Переменная, в которую придет значение от левого процесса
    double uleft; 
    MPI_Recv(&uleft, 1, MPI_DOUBLE, left, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Отправляем правому процессу значение u[me]
    MPI_Send(&u[me], 1, MPI_DOUBLE, right, 1, MPI_COMM_WORLD);
```

Добавим код пересылок и вернем код вычислений:
```c++
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <mpi.h>

double u0(double x) {
    if (x > 0.4 && x < 0.6)
        return 1;
    return 0;
}

void save(const int mb, const int me, const double h, const double *u, const std::string &filename) {
    std::ofstream f(filename);

    for (int m = mb; m <= me; m++) {
        double x = m*h;
        f << x << " " << u[m] << "\n";
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    if (argc < 2) {
        std::cerr << "USAGE: ./main <M>" << std::endl;
        return 1;
    }

    const int M = atoi(argv[1]);
    const double h = 1. / M;

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int Mi[size];
    for (int i = 0; i < size; i++)
        Mi[i] = M / size;
    for (int i = 0; i < M % size; i++)
        Mi[i]++;

    int mb, me;
    mb = 0;
    for (int i = 0; i < rank; i++)
        mb += Mi[i];
    me = mb + Mi[rank] - 1;

    std::cout << "Процесс " << rank << " обрабатывает элементы с "
        << mb << " до " << me << std::endl;

    double *u = new double[M];
    double *unext = new double[M];

    for (int m = mb; m <= me; m++) {
        double x = m*h;
        u[m] = u0(x);
    }

    const double tmax = 1;
    double sigma = 0.5;
    double dt = sigma * h;
    const int N = tmax / dt;
    dt = tmax / N;
    sigma = dt / h;

    for (int n = 0; n < N; n++) {
        // Теперь цикл не от 1 до M-1, а от mb+1 до me
        for (int m = mb+1; m <= me; m++)
            unext[m] = (1 - sigma) * u[m] + sigma * u[m-1];

        int left = (rank > 0) ? rank - 1 : size - 1;
        int right = (rank < size - 1) ? rank + 1 : 0;

        double uleft;
        MPI_Recv(&uleft, 1, MPI_DOUBLE, left, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&u[me], 1, MPI_DOUBLE, right, 1, MPI_COMM_WORLD);

        unext[mb] = (1 - sigma) * u[mb] + sigma * uleft;

        std::swap(u, unext);
    }

    save(mb, me, h, u, "par." + std::to_string(rank) + ".csv");

    delete[] u;
    delete[] unext;

    MPI_Finalize(); // Выключаем MPI

    return 0;
}
```

Компилируем этот код, запускаем и ... обнаруживаем, что он просто висит. Если вы запускали задачу на кластере,
она будет либо убита по истечении времени, либо необходимо удалить её из очереди счета вручную.

Дело в том, что все процессы одновременно начинают ожидать получения сообщения от соседа. В данной ситуации процессы
будут бесконечно ждать, и ни один из них так и не дойдет до того, чтобы отправить сообщение. Чтобы снять эту блокировку,
можно поменять последовательность отправки и приема: пусть все процессы сначала отправляют данные, а потом принимают:
```c++
    ...
    double uleft;
    MPI_Send(&u[me], 1, MPI_DOUBLE, right, 1, MPI_COMM_WORLD);
    MPI_Recv(&uleft, 1, MPI_DOUBLE, left, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    unext[mb] = (1 - sigma) * u[mb] + sigma * uleft;
    ...
```
После этого исправления код действительно работает и даже считает то, что нужно (проверьте!). Допустим, что задача
запускалась на трех процессах. В результате запуска должно было образоваться три файла
```bash
[v.pupkin@head ~/advection] ls par.*.csv
par.0.csv  par.1.csv  par.2.csv
```
Склеим эти файлы в один файл `par.csv` и откроем его в gnuplot:
```bash
[v.pupkin@head ~/advection] cat par.*.csv > par.csv
[v.pupkin@head ~/advection] gnuplot
gnuplot> plot 'par.csv'
```
Однако, хотя такой код и работает, он не очень надежен. Функция `MPI_Send`, в отличие от `MPI_Recv`, может закончить свою работу
еще до того, как сообщение будет принято на той стороне. Само сообщение при этом хранится в буфферах библиотеки MPI.
Если этого буффера не хватит (размер зависит от используемой конкретной библиотеки), функция `MPI_Send` будет ожидать в точности так же,
как и `MPI_Recv`.

Для решения этой проблемы необходимо гарантировать, что часть процессов начнут обмен данными с отправки, а часть --- с приёма сообщений.
Это можно реализовать самостоятельно, чередуя `MPI_Send/MPI_Recv` в зависимости от номера процесса `rank`, но гораздо лучше поручить
эту задачу библиотеке MPI: функция `MPI_Sendrecv` предназначена именно для этих целей: она позволяет объединить последовательнные
вызовы `MPI_Send` и `MPI_Recv` так, чтобы блокировки гарантированно не произошло. Эта функция принимает совокупность аргументов, переданных
в отдельности `MPI_Send` и `MPI_Recv`:
```c++
    double uleft;
    MPI_Sendrecv(
        &u[me], 1, MPI_DOUBLE, right, 1, // параметры Send
        &uleft, 1, MPI_DOUBLE, left, 1,  // параметры Recv
        MPI_COMM_WORLD, MPI_STATUS_IGNORE // общие параметры
    );

    unext[mb] = (1 - sigma) * u[mb] + sigma * uleft;
```
