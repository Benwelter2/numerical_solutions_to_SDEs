# Numerical Solutions to SDEs

This repository contains the numerical implementations and convergence analysis for the course **Student Seminar 2** at the **University of Luxembourg (Academic Year 2025-2026)**.

The project focuses on approximating solutions to SDEs using the **Euler** and **Milstein** schemes.

$$X_{i+1} = X_{i} + a(X_{i})\Delta t + b(X_{i})\Delta W_{i}$$
$$X_{i+1} = X_{i} + a(X_{i})\Delta t + b(X_{i})\Delta W_{i} + \frac{1}{2}b(X_{i})b'(X_{i})[(\Delta W_{i})^{2} - \Delta t]$$

| Scheme | Theoretical Strong | Experimental Strong (Slope) | Theoretical Weak | Experimental Weak (Slope) |
| :--- | :---: | :---: | :---: | :---: |
| **Euler** | $0.5$ |$\approx 0.496$  | $1.0$ | $\approx 0.980$ |
| **Milstein** |$1.0$  |$\approx 0.992$  | $1.0$ | $\approx 0.995$ |
