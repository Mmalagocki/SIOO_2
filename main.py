from sympy import *
import numpy as np
import math

########################## GLOBALS #############################
a = Symbol('a')
b = Symbol ('b')
c = Symbol('c')
d = Symbol ('d')
e = Symbol('e')
f = Symbol ('f')
g = Symbol('h')
h = Symbol ('h')
i = Symbol('i')
j = Symbol ('j')
results = []
how_many_uknowns = 0
######################### CONVERSION #############################
def convert_function(fun):
    fun = sympify(fun, evaluate = False)
    return fun

def convert_values(values_string):
    values_array = []
    values_string.replace( " ", "" )
    separator = values_string.find("," , 0)
    values_array.append(float(values_string[0:separator]))
    values_array.append(float(values_string[separator+1:]))
    return values_array

########################## DERIVATIVES #########################
class gradient:
    def __init__(self, function, value, letter):
        self.function = function
        self.value = value
        self.letter = letter
        self.derivative = 0
        self.d_fun = 0
        self.direction = self.calculate_gradient()

    def calculate_derivative(self):
        function_str = str(self.function)
        if (function_str.find(self.letter , 0) != -1):
            xprime = self.function.diff(self.letter)
            #self.derivatives.append(xprime)
            self.derivative = xprime
            self.d_fun = xprime
            return 0

        else: 
            return 0

    def calculate_gradient(self):
        result = 0
        self.calculate_derivative()
        letter_value = self.value
        letter_value = str(letter_value)
        self.derivative = str(self.derivative)
        self.derivative = self.derivative.replace( " ", "" )
        self.derivative = self.derivative.replace(self.letter, letter_value)
        result = eval(self.derivative)
        return result
    
    def get_direction(self):
        return self.direction
    
    def get_derivative_fun(self):
        return self.d_fun

########################## FLETCHER-REEVES #########################  
class Fletcher_Reeves:
    def __init__(self, points, e, fun, values, n_order):
        points = convert_values(points)
        return_results = []
        self.eta = []
        self.points = values
        self.epsylon = e
        self.k = 1
        self.fun = fun
        self.d_fun = []
        self.grads = [None] * how_many_uknowns
        self.values = values
        #print("init values ", self.values)
        local_range = range(0 ,how_many_uknowns)
        for ik in local_range:
            self.values.append(None)
        self.ds = [None]*how_many_uknowns
        self.calculate()
        
    def calculate(self):
        letter = 'a'
        local_range = range(0 ,how_many_uknowns)

        for ik in local_range:
            grad = gradient(self.fun, self.values[ik], chr(ord(letter) + ik))
            self.grads[ik] = gradient.get_derivative_fun(grad)
            d = -gradient.get_direction(grad)
            self.ds[ik] = d
            self.d_fun.append(gradient.get_derivative_fun(grad))

        range_begin = 0
        range_end = 0
        norm = 10
        while ( norm > self.epsylon):
            for i in local_range:
                range_begin =  self.points[i]
                range_end =  self.points[i]
            str_fun = self.create_new_fun(self.fun, letter)
            dich = dichotomy(str_fun, range_begin, range_end)

            results = dichotomy.get_results(dich)
            results = np.average(results, axis = 0)
            point_value = self.calculate_value_in_point(results)

            norm = math.sqrt(pow(self.ds[0], 2) + pow(self.ds[1], 2))
            self.calculate_next_argument(results)
            self.count_eta()

            self.count_next_d(point_value)

        return 1

    def count_next_d(self, point_value):
        local_range = range(0 ,how_many_uknowns)
        letter = 'a'
        for i in local_range:
            letter_value = self.values[i]
            letter_value = str(letter_value)
            derivative = str(self.grads[i])
            derivative = derivative.replace( " ", "" )
            derivative = derivative.replace(chr(ord(letter) + i), letter_value)
            grad_value = eval(derivative)
            dk1 = -grad_value + self.eta[i]*self.ds[i]
            self.ds[i] = dk1

    def count_eta(self):
        local_range = range(0 ,how_many_uknowns)
        letter = 'a'
        for i in local_range:
            grad1 = gradient(self.fun, self.values[i], chr(ord(letter) + i))
            xk_value = gradient.get_direction(grad1)
            grad2 = gradient(self.d_fun[i], self.values[i], chr(ord(letter) + i))
            xk_plus_1_value = gradient.get_direction(grad2)
            single_eta = (xk_plus_1_value * xk_plus_1_value) / (xk_value * xk_value)
            self.eta.insert(i, single_eta)

    def calculate_next_argument(self, results):
        local_range = range(0 ,how_many_uknowns)
        for ik in local_range:
            prev_val = self.values[ik]
            self.values[ik + how_many_uknowns] = prev_val
            next_arg = self.values[ik] + results*self.ds[ik]
            self.values[ik] = round(next_arg, 2)

    def calculate_value_in_point(self, results):
        function_str = str(self.fun)
        function_str = function_str.replace( "alpha", str(results) )
        function_str = function_str.replace( " ", "" )
        value = compile(function_str, "<string>", "eval")
        return eval(value)

    def create_new_fun(self, fun, letter):
        function_str = str(fun)
        local_range = range(0 ,how_many_uknowns)
        for i in local_range:
            if (function_str.find(chr(ord(letter) + i), 0) != -1):
                letter_value = str(self.values[i])
                d_value = str(self.ds[i])
                new_value = letter_value + "+" + "alpha" + "*" + d_value
                fun = str(fun)
                fun = str(fun.replace(chr(ord(letter) + i), new_value))
            else: 
                print('Couldn"t create new function')
                return 0
        return fun

########################## DICHOTOMY #########################
class dichotomy:
    def __init__(self, fun, rb, re):
        self.function = fun
        self.rb = rb
        self.re = re
        self.minimum = 0
        self.result = self.calculate()

    def calculate(self):
        results = [None]*how_many_uknowns
        results_array_len = 0
        tolerance = 0.00001
        range_begin = self.rb
        range_end = self.re
        distance = np.abs(range_begin-range_end) 
        i = 0
        while (distance >= tolerance):
            delta = distance/ 4
            cl = 0.5 * (range_begin + range_end) - delta
            cr = 0.5 * (range_begin + range_end) + delta               
            frb = self.f(range_begin)
            fre = self.f(range_end)
            fcl = self.f(cl)
            fcr = self.f(cr)              
        
            if ((frb >= fcl) and (frb >= fre) and (fcl <= fcr)) :
                range_begin = cl
            elif ((frb >= fcl) and (frb >= fre) and (fcl >= fcr)) :
                range_begin =  cl
            elif ((frb <= fcl) and (frb <= fre) and (fcl >= fcr)):
                range_begin = cl    
            elif ((frb >= fcl) and (frb >= fre) and (fcl >= fcr)):
                range_begin = cl
            elif ((frb >= fcl) and (frb <= fre) and (fcl <= fcr)):
                range_end = cr
            elif ((frb <= fcl) and (frb <= fre) and (fcl <= fcr)):
                range_end = cr                    
            elif (fre < fcl) and (fcr> frb):
                self.not_unimodal(range_begin, cl)
                self.not_unimodal(cr, range_end)
                self.not_unimodal(cl, cr)
            elif ((frb >= fcl) and (frb < fre) and (fcl <= fcr)):
                range_begin = cl              
            elif ((frb >= fcl) and (frb < fre) and (fcl > fcr)):
                range_end = cr                    
            distance = np.abs(range_begin-range_end)
        
        results_array_len += 1
        results[results_array_len-1] = range_begin
        results_array_len += 1
        results[results_array_len-1] = range_end
        
        self.minimum = (abs(range_begin) + abs(range_end))/2
        return results

    def not_unimodal(self, rb, re):
        tolerance = 0.00001
        range_begin = rb
        range_end = re
        distance = np.abs(range_begin-range_end) 
    
        while (distance >= tolerance):
            delta = distance/ 4
            cl = 0.5 * (range_begin + range_end) - delta
            cr = 0.5 * (range_begin + range_end) + delta               
            frb = self.f(range_begin)
            fre = self.f(range_end)
            fcl = self.f(cl)
            fcr = self.f(cr)              
    
            if ((frb >= fcl) and (frb >= fre) and (fcl <= fcr)) :
                range_begin = cl
            elif ((frb >= fcl) and (frb >= fre) and (fcl >= fcr)) :
                range_begin =  cl
            elif ((frb <= fcl) and (frb <= fre) and (fcl >= fcr)):
                range_begin = cl
            elif ((frb >= fcl) and (frb >= fre) and (fcl >= fcr)):
                range_begin = cl
            elif ((frb >= fcl) and (frb <= fre) and (fcl <= fcr)):
                range_end = cr
            elif ((frb <= fcl) and (frb <= fre) and (fcl <= fcr)):
                range_end = cr
            elif (fre < fcl) and (fcr> frb):
                self.not_unimodal(range_begin, cl)
                self.not_unimodal(cr, range_end)
                self.not_unimodal(cl, cr)
            elif ((frb >= fcl) and (frb < fre) and (fcl <= fcr)):
                range_begin = cl              
            elif ((frb >= fcl) and (frb < fre) and (fcl > fcr)):
                range_end = cr                    
    
        distance = np.abs(range_begin-range_end)
        results_array_len += 1
        results.insert(results_array_len, round(range_begin, 10))
        results_array_len += 1
        results.insert(results_array_len, round(range_end, 10))

    def f(self, x):
        function_str = str(self.function)
        function_str = function_str.replace( "alpha", str(x) )
        function_str = function_str.replace( " ", "" )
        code = compile(function_str, "<string>", "eval")
        return eval(code)

    def get_results(self):
        #print(self.result)
        return self.result


########################## MAIN #########################
'''
print("Pass the function")
function_str = input()
print("Pass the point")
point = input()
print("Pass e")
e = input()
print("Pass the parameter or parameters")
parameters = input()
'''
function_str = "2*a**2 + 3*b**2"
points = "19, 20" 
e = "1"
parameters = '20,30'
how_many_uknowns = 2


fun = convert_function(function_str)
parameters = convert_values(parameters)
FR = Fletcher_Reeves(points, float(e), (fun), parameters, how_many_uknowns)


