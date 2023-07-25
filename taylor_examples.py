import sympy as sp
import math
from poly import Polynomial, identity, zero
from taylor_expand import *
from interval import interval, imath


if __name__ == '__main__':
    # Define the function sin(x)
    x = sp.Symbol('x')
    f = sp.sin(x)

    # Calculate the first 5 Taylor expansions of sin(x) at x=0
    taylor_series = f.series(x, 0, 5).removeO()

    # output result
    print(taylor_series)

    x = sp.Symbol('x')
    f = sp.cos(x)

    taylor_series = f.series(x, 1, 2).removeO()
  
    #print(taylor_series)

def double_factorial(n: int):
    if n <= 0:
        return 1
    else:
        return n * double_factorial(n - 2)
        
def coefficient(function_name: str, n: int, x):
    func = imath
    if function_name == 'exp':
        return func.exp(x) / math.factorial(n)
    elif function_name == 'log':
        if n == 0:
            return func.log(x)
        else:
            return (-1) ** (n - 1) / ((x ** n) * n)
    elif function_name == 'sqrt':
        if n == 0:
            return func.sqrt(x)
        else:
            return func.sqrt(x) * ((-1) ** (n - 1)) * double_factorial(2 * n - 3) / (
                    math.factorial(n) * (2 ** n) * x ** n)
    elif function_name in ['sin', 'cos']:
        i = n % 4 if function_name == 'sin' else (n + 1) % 4
        if i == 0:
            return func.sin(x) / math.factorial(n)
        elif i == 1:
            return func.cos(x) / math.factorial(n)
        elif i == 2:
            return -func.sin(x) / math.factorial(n)
        else:
            return -func.cos(x) / math.factorial(n)
    elif function_name in ['sinh', 'cosh']:
        i = n % 2 if function_name == 'sinh' else (n + 1) % 2
        if i == 0:
            return func.sinh(x) / math.factorial(n)
        else:
            return func.cosh(x) / math.factorial(n)
    elif function_name == 'atan':
        if x == 0 or x == interval[0]:
            if n % 2 == 0:
                return 0
            else:
                return (-1) ** ((n - 1) / 2) / n
        else:
            if n == 0:
                return func.atan(x)
            elif n % 2 == 0:
                return (-1) ** (n - 1) / ((1 + x ** 2) ** (int(n / 2)) * n) * func.sin(n * (func.atan(1 / x)))
            else:
                return (-1) ** (n - 1) / (((1 + x ** 2) ** (int(n / 2)) * func.sqrt(1 + x ** 2)) * n) * func.sin(
                    n * (func.atan(1 / x)))



## R_n(x) = f^(n+1)(c) * (x - a)^(n+1) / (n+1)!
def element(func, n, center, domain):
    #coef_list = [0, 1, 0, -1/6]
    coef_list = []
    for i in range(n + 1):
        coef_list.append(coefficient(func, i, center))
    poly = Polynomial(max_order=n,params=coef_list)
    #coef1 = imath.sin(interval([-1, 1])) / math.factorial(4)
    coef1 = coefficient(func, n + 1, domain)
    func = eval('imath.' + func)
    if coef1[0].inf >= 0 or coef1[0].sup <= 0:
        a = interval(domain[0].inf)
        b = interval(domain[0].sup)
        error1 = func(a) - poly.bound_best(a - center, ['native', 'hornor'])
        error2 = func(b) - poly.bound_best(b - center, ['native', 'hornor'])
        error0 = func(interval(center)) - poly.bound(interval(0))
        error = (interval.hull((error1, error2, error0))) & (coef1 * ((domain - center) ** (n + 1)))
    else:
        error = coef1 * ((domain - center) ** (n + 1))
    return TaylorModel(
        poly=poly,
        error=error,
        center=center,
        domain=domain
    )