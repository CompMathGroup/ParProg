Создадим директорию для кода и программы:

```bash
[v.pupkin@head ~]$ mkdir simple
[v.pupkin@head ~]$ cd simple
[v.pupkin@head ~/simple]$
```

Отредактируем файл `main.cpp`, который будет содержать весь исполняемый код программы. Мы собираемся использовать библиотеку MPI, так что необходимо подключить заголовочный файл, в котором описаны все имеющиеся в библиотеке MPI функции. Далее идет исходный код программы `main.cpp`.

```c++
#include <mpi.h>
```

Весь код программы будет состоять из единственной функции `main`, принимающей стандартные аргументы `argc` и `argv`:

```c++
#include <mpi.h>

int main(int argc, char *argv[]) {
    // код функции main
}
```

Библиотеку MPI необходимо инициализировать (включить) в начале и финализировать (выключить) в конце. Делается это стандартным образом, вызовами функций `MPI_Init` и `MPI_Finalize`:

```c++
#include <mpi.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    // код, который может использовать библиотеку MPI
    MPI_Finalize();
}
```

MPI может при необходимости менять параметры `argc` и `argv` у функции `main`, извлекая лишние аргументы (предназначенные для самой библиотеки MPI), которые, возможно, были добавлены при запуске нескольких задач с помощью `mpirun`. Конкретные значения сильно зависят от конкретной реализации MPI, но инициализация всегда производится именно таким образом.

При работе MPI задачи запускается несколько процессов, каждый из которых является копией одной и той же программы. Эти процессы решают задачу совместно, общаясь посредством передачи сообщений (отсюда название библиотеки MPI — Message Passing Interface). Каждый процесс получает уникальный идентификатор, позволяющий коду отличать между собой копии программы. Процессы могут общаться в пределах так называемых коммуникаторов. Для многих задач достаточно всего одного коммуникатора, существующего в любой программе — `MPI_COMM_WORLD`. Этот коммуникатор объединяет все процессы в задаче. MPI позволяет узнать, сколько процессов объединяет данный конкретный коммуникатор, а также для каждого процесса узнать его уникальный номер в рамках данного коммуникатора. Для этого служат функции `MPI_Comm_size` и `MPI_Comm_rank`.

```c++
#include <mpi.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Записать в size количество процессов в коммуникаторе MPI_COMM_WORLD
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Записать в rank идентификатор процесса в коммуникаторе MPI_COMM_WORLD

    MPI_Finalize();
}
```

Необходимо понимать, что каждая запущенная копия программы содержит свои переменные `size` и `rank`. Переменная 
`size` у всех процессов будет иметь одинаковые значения, а `rank` будет у всех разная, принимающая значения от `0` до `size-1`. Распечатаем эти значения. Для этого подключим дополнительно заголовочный файл `iostream`, а выведем с помощью `std::cout`:

```c++
#include <mpi.h>
#include <iostream>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Записать в size количество процессов в коммуникаторе MPI_COMM_WORLD
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Записать в rank идентификатор процесса в коммуникаторе MPI_COMM_WORLD

    std::cout << "I'm process #" << rank << " of " << size << " processes" << std::endl;

    MPI_Finalize();
}
```

Данный код необходимо собрать с помощью компилятора `mpic++`. Этот компилятор использует стандартный компилятор `c++`, просто добавляет опции, указывающие где искать заголовочный файл `mpi.h` а также присоединяет библиотеку `libmpi.so` к коду. Для компиляции выполним

    mpic++ main.cpp -o main

Здесь мы указали с помощью ключа `-o` имя выходного файла `main`. По историческим причинам именем по умолчанию является `a.out`. Задачу уже можно попробовать запустить, однако управляющий узел не предназначен для запуска MPI задач и в выводе будет содержаться некоторое число ошибок:

```bash
[v.pupkin@head ~/simple]$ ./main 
[head.control:14092] mca: base: component_find: unable to open /usr/local/lib/openmpi/mca_btl_usnic: libpsm_infinipath.so.1: cannot open shared object file: No such file or directory (ignored)
[head.control:14092] mca: base: component_find: unable to open /usr/local/lib/openmpi/mca_mtl_psm: libpsm_infinipath.so.1: cannot open shared object file: No such file or directory (ignored)
[head.control:14092] mca: base: component_find: unable to open /usr/local/lib/openmpi/mca_mtl_ofi: libpsm_infinipath.so.1: cannot open shared object file: No such file or directory (ignored)
--------------------------------------------------------------------------
[[35854,1],0]: A high-performance Open MPI point-to-point messaging module
was unable to find any relevant network interfaces:

Module: OpenFabrics (openib)
  Host: head

Another transport will be used instead, although this may result in
lower performance.
--------------------------------------------------------------------------
I'm process #0 of 1 processes
```

В последней строчке находится интересующий нас вывод. Была запущена одна копия программы, получившей номер `0`.
Для того, чтобы запустить несколько копий программы, составим pbs файл для задачи. Создадим `simple.pbs` с содержанием

```bash
#PBS -N simple
#PBS -l nodes=2:ppn=4
#PBS -l walltime=00:01:00

cd simple
/usr/local/bin/mpirun ./main
```

Здесь мы описываем задачу с именем `simple`, запрашиваем 2 узла по 4 ядра на узел, обещаем не более 1 минуты выполнения программы. Задача запускается на вычислительных узлах в домашней директории, так что необходимо прямо в pbs файле переходить в директорию simple, содержащую программу main. 

Запускаем задачу командой `qsub simple.pbs` и ждем, пока она завершится (обычно на запуск задачи уходит несколько секунд) запуская `qstat`. Как только задача выполнится, в файле `simple.o<номер задачи>` будет содержаться вывод от всех 8 процессов:

```
[v.pupkin@head ~/simple]$ cat simple.o4851
I'm process #0 of 8 processes
I'm process #1 of 8 processes
I'm process #2 of 8 processes
I'm process #3 of 8 processes
I'm process #6 of 8 processes
I'm process #7 of 8 processes
I'm process #4 of 8 processes
I'm process #5 of 8 processes
```

Каждый процесс сообщил свой номер, но делают это они одновременно, так что последовательность строк может отличаться.