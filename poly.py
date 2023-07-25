from __future__ import annotations
import copy
import math
from interval import interval



class Polynomial:
    def __init__(self, max_order: int, params: list):
        self.max_order = max_order
        self.params = [0 for i in range(max_order + 1)]
        for index in range(min(len(params), max_order + 1)):
            self.params[index] = params[index]

    def __str__(self):
        return f'Polynomial(max order:{self.max_order},params:{self.params})'

    def __repr__(self):
        return f'Polynomial(max order:{self.max_order},params:{self.params})'

# First judge whether other is a polynomial, 
# if not, it must be converted to a poly

    def __add__(self, other):
        if not isinstance(other, Polynomial):
            other = Polynomial(0, [other])
        # max_poly_shape
        left, right, max_order = max_poly_shape(self, other)
        params = [left.params[i] + right.params[i] for i in range(max_order + 1)]
        return Polynomial(max_order, params)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if not isinstance(other, Polynomial):
            other = Polynomial(0, [other])
        # max_poly_shape
        left, right, max_order = max_poly_shape(self, other)
        params = [left.params[i] - right.params[i] for i in range(max_order + 1)]
        return Polynomial(max_order, params)

    def __rsub__(self, other):
        if isinstance(other, (int, float, interval)):
            other = Polynomial(0, [other])
        return other - self

    def __mul__(self, other):
        if not isinstance(other, Polynomial):
            if other != 0:
                other = Polynomial(0, [other])
            else:
                return 0
        # sum_poly_shape
        left, right, max_order = sum_poly_shape(self, other)
        params = []
        for i in range(max_order + 1):
            if i == 0:
                cur_param = left.params[0] * right.params[0]
            else:
                cur_param = 0
                for index in range(i + 1):
                    
                    cur_param += left.params[index] * right.params[i - index]
                    
            params.append(cur_param)
        return Polynomial(max_order, params)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        
        # polynomial division
        
        if not isinstance(other, (int, float, interval)):
            raise Exception('Can only divide non-polynomials')
        if other == 0:
            raise Exception('divisor cannot be 0')
        params = [param / other for param in self.params]
        return Polynomial(self.max_order, params)

    

    # Judging whether it is equal between polys
    def __eq__(self, other: Polynomial):
        if self.max_order != other.max_order:
            return False
        for index in range(self.max_order + 1):
            if self.params[index] != other.params[index]:
                return False
        return True
    
  # negtive
    def __neg__(self):
        return 0 - self
  
    def square(self):
        return self * self
    # power
    def __pow__(self, power, modulo=None):  
        max_order = self.max_order * power
        if power == 0:
            return identity(max_order)
        elif power == 1:
            return self
        elif power % 2 == 0:
            if power == 2:
                return self.square()
            else:
                half = int(power / 2)
                return Polynomial(max_order, self.__pow__(half).square().params)
        else:
            half = int(power / 2)
            return Polynomial(max_order, (self * self.__pow__(half).square()).params)



    def bound(self, x, method='native'):
        if method == 'root':
            return self.bound_root(x)
        elif method == 'hornor':
            return self(x)
        # doesn’t match is native
        return self.bound_native(x)


    def bound_best(self, x, methods):
        bound = self.bound(x, methods[0])
        for m in methods:
            bound = bound & self.bound(x, m)
        return bound

    def bound_root(self, x):
        derivate = self.derivative()
        bound = derivate.bound_native(x)
        # x[0].inf is the lower bound. x[0].sup is the upper bound
        if bound[0].inf >= 0 or bound[0].sup <= 0:
            return interval[interval(self(x[0].inf)), interval(self(x[0].sup))]
        else:
            derivate_2 = derivate.derivative()
            roots = x.newton(derivate, derivate_2)     
            bound = interval[interval(self(x[0].inf)), interval(self(x[0].sup))]
            # interval.hull method can take the union of multiple intervals
            for root in roots:  
                bound = interval.hull((bound, self(interval(root))))  
            return bound

    def bound_native(self, x):
        b = self.params[0]
        max_order = self.max_order
        for i in range(1, max_order + 1):
            b += self.params[i] * (x ** i)
        return b

    # find the first derivative
    def derivative(self):
        max_order = self.max_order
        params = [0 for _ in range(max_order)]
        params[0] = self.params[1]
        for index in range(2, max_order + 1):
            params[index - 1] = index * self.params[index]
        return Polynomial(max_order - 1, params)

    # Remove the coefficient after the nth order    
    def remove(self, n):
        return Polynomial(n, self.params[:n + 1])
    
    # Horner's rule
    def __call__(self, dx):
        
        max_order = self.max_order
        max_param = self.params[max_order]
        for index in range(max_order, 0, -1):
            max_param = max_param * dx + self.params[index - 1]
        return max_param


def sum_poly_shape(left: Polynomial, right: Polynomial):
    """
    Transform the maximum order of the two polynomials 
    left and right into the sum of the two
    """
    max_order = left.max_order + right.max_order
    left_params = copy.deepcopy(left.params)
    right_params = copy.deepcopy(right.params)
    left_params += [0 for _ in range(max_order - left.max_order)]
    right_params += [0 for _ in range(max_order - right.max_order)]
    return Polynomial(max_order, left_params), Polynomial(max_order, right_params), max_order


def max_poly_shape(left: Polynomial, right: Polynomial):
    """
    Transform the maximum order of the two polynomials 
    left and right into the maximum of the two
    """
    max_order = max(left.max_order, right.max_order)
    left_params = copy.deepcopy(left.params)
    right_params = copy.deepcopy(right.params)
    if left.max_order < max_order:
        left_params += [0 for _ in range(max_order - left.max_order)]
    if right.max_order < max_order:
        right_params += [0 for _ in range(max_order - right.max_order)]
    return Polynomial(max_order, left_params), Polynomial(max_order, right_params), max_order

def identity(max_order):
    params = [0 for i in range(max_order + 1)]
    params[0] = 1
    return Polynomial(max_order, params)


def zero(max_order):
    params = [0 for i in range(max_order + 1)]
    return Polynomial(max_order, params)


if __name__ == '__main__':
    p1 = Polynomial(3, [1.2593, -0.3962, -0.1539])
    p2 = Polynomial(1, [3, 2])
    #print(p1)
    # print(p1 + p2)
    # print(p1 - p2)
    #print(p1 * p2)
    # print(p1.square())
    # print(p1 * p1 * p1)
