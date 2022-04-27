# Numerical-methods
- # Task №1
My function is: z(x) = cos(2.8x + sqrt(1 + x)) * arctan(1.5x + 0.2), x = 0.1(0.01)0.2
Firstly, I calculated methodological errors for internal functions in file "20Б01_Проценко_Погрешности№1.pdf". 
Secondly, in "zad1.java" file, I implemented the root function through the Heron iterative formula and the cos and arctan functions through the expansion in a Taylor series.
Then I find the values of the internal functions and the external one with the previously calculated errors and compare them with the results of the library functions.
I output the results in Excel using the org.apache.poi library, they can be found in the file "20.Б01_Проценко_Погрешности№1.xls".
- # Task №2
In "zad2.java" file, I implemented the basic functions for working with matrices and four methods for finding a solution to a system of linear equations: Simple iteration method, Seidel's method, Householder method and Gauss method. 
And I output the results in Excel using the org.apache.poi library into the file "20.Б01_Проценко_СЛАУ№2.xls".
