from sympy import *


def f(x):
    return 6 * x ** 3 + 3 * x * x + x + 4


def g(x):
    return 18 * x * x + 6 * x + 1


x1 = Symbol('x1')
x2 = Symbol('x2')
m = Symbol('m')
a = Symbol('a')
b = Symbol('b')
c = Symbol('c')
sz = solve([-a + b - c - m,
            14 - a - b - c + m,
            f(x2) - a - b * x2 - c * x2 * x2 - m,
            f(x1) - a - b * x1 - c * x1 * x1 + m,
            g(x1) - b - 2 * c * x1,
            g(x2) - b - 2 * c * x2],
           (a, b, c, m, x1, x2))
for i in sz:
    if -1 < i[4] < i[5] < 1:
        print(i)