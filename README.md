# Numerical Solutions to SDEs

This repository contains the numerical implementations and convergence analysis for the course **Student Seminar 2** at the **University of Luxembourg (Academic Year 2025-2026)**.

The project focuses on approximating solutions to SDEs using 

**Euler Scheme:** 
$$X_{i+1} = X_{i} + a(X_{i})\Delta t + b(X_{i})\Delta W_{i}$$, and

**Milstein Scheme:**
$$X_{i+1} = X_{i} + a(X_{i})\Delta t + b(X_{i})\Delta W_{i} + \frac{1}{2}b(X_{i})b'(X_{i})[(\Delta W_{i})^{2} - \Delta t]$$.

The table below summarizes the theoretical and experimental order od convergence (o.o.c.).
| Scheme | th. strong o.o.c. | exp. strong o.o.c. | th. weak o.o.c. | exp. weak o.o.c. |
| :--- | :---: | :---: | :---: | :---: |
| **Euler** | $0.5$ | $0.496$  | $1.0$ | $0.980$ |
| **Milstein** |$1.0$  |$0.992$  | $1.0$ | $0.995$ |
