from poly import Polynomial, identity, zero
from interval import interval, imath
from taylor_examples import *
import sys



class TaylorModel:
    def __init__(self, poly: Polynomial, center: interval, error: interval, domain: interval):
        # polynomial expression
        self.poly = poly
        # error bound
        self.error = interval(error)
        # center point
        self.center = interval(center)
        # domain
        self.domain = interval(domain)

        # Convert all coefficients to interval form
        self.poly.params = [interval[param] for param in self.poly.params]

        self.max_order = self.poly.max_order

    def __repr__(self):
        return f'Taylor Model:Poly : {self.poly},Center: {self.center},Error: {self.error},Domain: {self.domain}'

    def __str__(self):
        return f'Taylor Model:Poly : {self.poly},Center: {self.center},Error: {self.error},Domain: {self.domain}'

    def __neg__(self):
        return TaylorModel(
            poly=-self.poly,
            center=self.center,
            error=self.error,
            domain=self.domain
        )

    def __add__(self, other):
        """
        addition
        :param other:
        :return:
        """
        if isinstance(other, (int, float, interval)):
            other = TaylorModel(Polynomial(0, [other]), 0, 0, self.domain)
        return TaylorModel(
            poly=self.poly + other.poly,
            center=self.center,
            error=self.error + other.error,
            domain=self.domain & other.domain
        )

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        """
        subtraction
        """
        if isinstance(other, (int, float, interval)):
            other = TaylorModel(Polynomial(0, [other]), 0, 0, self.domain)
        return TaylorModel(
            poly=self.poly - other.poly,
            center=self.center,
            error=self.error - other.error,
            domain=self.domain & other.domain
        )

    def __rsub__(self, other):
        return self - other

    def __mul__(self, other):
        if not isinstance(other, (int, float, interval)):
            return self.muliply_helper(
                other=other,
                result_order=max(self.max_order, other.max_order)
            )
        if other == 0:
            return 0
        return TaylorModel(
            poly=self.poly * other,
            center=self.center,
            error=other * self.error,
            domain=self.domain
        )

    def __rmul__(self, other):
        return self * other

    def muliply_helper(self, other: "TaylorModel", result_order: int = -1, method='root'):
        """
        Multiplication helper function
        :param result_order: Integer, if it is -1 or greater 
          than the highest order, all parameters are reserved
        """
        poly = self.poly * other.poly
        poly_final = poly.remove(result_order)
        poly_trun = poly - poly.remove(result_order)
        
        bound0 = poly_trun.bound(self.domain, method)
        bound1 = self.poly.bound(self.domain, method) * other.error
        bound2 = other.poly.bound(other.domain, method) * self.error  
        error = bound0 + bound1 + bound2 + self.error * other.error
        taylor = TaylorModel(
            poly=poly_final,
            center=self.center,
            error=error,
            domain=self.domain & other.domain
        )
        return taylor

    

    def square(self):
        return self * self

    def identity(self):
        return TaylorModel(
            poly=identity(self.max_order),
            center=self.center,
            error=self.error,
            domain=self.domain
        )

    def __pow__(self, power: int, modulo=None):
        if power == 0:
            return self.identity()
        if power == 1:
            return self
        if power % 2 == 0:
            half = int(power / 2)
            return self.__pow__(half).square()
        half = int(power / 2)
        return self * self.__pow__(half).square()

    def zero_taylor(self):
        return TaylorModel(
            poly=zero(self.poly.max_order),
            center=self.center,
            error=self.error,
            domain=self.domain
        )

    def remove(self, n: int):
        poly = self.poly.remove(n)
        zero_poly = self.zero_taylor()
        zero_poly.poly.params[n + 1:] = self.poly.params[n + 1:]
        error = zero_poly.poly.bound(self.domain) + self.error
        return TaylorModel(
            poly=poly,
            center=self.center,
            error=error,
            domain=self.domain
        )

    def bound(self):
        """
        bound
        """
        return self.poly.bound(self.domain, method='root') + self.error

    def __call__(self):
        return self.bound()


sys.setrecursionlimit(100000) 


if __name__ == '__main__':
    
    operations = input("Which kind of operation?(addition, multi, power)")
    if operations == 'addition' or operations == 'multi':
        m_order = int(input("input max order of poly1: "))
        i = 0
        coe_list = []
        while i <= m_order:
            o = input("coefficients of poly1: ")
            coe_list.append(o)
            i += 1
        p = Polynomial(
            max_order = m_order,
            params = coe_list
        )
        error1 = float(input("error_lower_bound of tm1: "))
        error2 = float(input("error_upper_bound of tm1: "))
        center = float(input("center of tm1: "))
        domain1 = float(input("domain_lower_bound of tm1: "))
        domain2 = float(input("domain_upper_bound of tm1: "))
        tm1 = TaylorModel(
            poly=p,
            error=interval([error1, error2]),
            center=interval([center, center]),
            domain=interval([domain1, domain2])
        )
        m_order2 = int(input("input max order of poly2: "))
        i2 = 0
        coe_list2 = []
        while i2 <= m_order2:
            o2 = input("coefficients of poly2: ")
            coe_list2.append(o2)
            i2 += 1    
        p2 = Polynomial(
            max_order = m_order2,
            params = coe_list2
        )
        error3 = float(input("error_lower_bound of tm2: "))
        error4 = float(input("error_upper_bound of tm2: "))
        center2 = float(input("center of tm2: "))
        domain3 = float(input("domain_lower_bound of tm2: "))
        domain4 = float(input("domain_upper_bound of tm2: "))
        tm2 = TaylorModel(
            poly=p2,
            error=interval([error3, error4]),
            center=interval([center2, center2]),
            domain=interval([domain3, domain4])
        )
        if operations == 'addition':
            print('tm1: ', tm1)
            print('tm2: ', tm2)
            print('tm1 + tm2: ',tm1 + tm2)
        elif operations == 'multi':
            print('tm1: ',tm1)
            print('tm1: ',tm2)
            print('tm1 * tm2: ', tm1 * tm2)

    elif operations == 'power':
        m_order3 = int(input("input max order of poly: "))
        i3 = 0
        coe_list3 = []
        while i3 <= m_order3:
            o = input("coefficients of poly: ")
            coe_list3.append(o)
            i3 += 1   
        p = Polynomial(
            max_order = m_order3,
            params = coe_list3
        )
        error1 = float(input("error_lower_bound of tm: "))
        error2 = float(input("error_upper_bound of tm: "))
        center = float(input("center of tm: "))
        domain1 = float(input("domain_lower_bound of tm: "))
        domain2 = float(input("domain_upper_bound of tm: "))
        tm = TaylorModel(
            poly=p,
            error=interval([error1, error2]),
            center=interval([center, center]),
            domain=interval([domain1, domain2])
        )
        power = int(input("What is the power?"))
        print('tm: ',tm)
        print('power: ',power)
        print(tm ** power)
      

      

     
     

      

     
     
