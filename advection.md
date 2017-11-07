---
layout: concise
title: Параллельное решение уравнения переноса
---

# Постановка задачи и разностная схема

Рассмотрим уравнение переноса с периодическими граничными условиями

$$
\begin{aligned} &\frac{\partial u}{\partial t} + \frac{\partial u}{\partial x} = 0, \qquad x \in [0, 1]\\ &u\big|_{t=0} = u_0(x)\\ &u\big|_{x=0} = u\big|_{x=1} \end{aligned}
$$

Для дискретизации этой задачи используем схему «явный левый уголок»

$$
\begin{aligned} &\frac{u^{n+1}_m - u^n_m}{\tau} + \frac{u^n_m - u^n_{m-1}}{h} = 0, &\quad m = 1, 2, \dots, M - 1\\
&\frac{u^{n+1}_0 - u^n_0}{\tau} + \frac{u^n_0 - u^n_{M-1}}{h} = 0, &\quad m = 0\\ &u^0_m = u_0(x_m), &\qquad m = 0, \dots, M-1 \end{aligned}
$$

# Последовательная программа
