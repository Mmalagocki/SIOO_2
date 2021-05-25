from sympy import *
import numpy as np
from numpy import linalg as LA
import math
import scipy.optimize as opt
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
NavigationToolbar2Tk)
from tkinter import *
import matlab.engine
from decimal import Decimal
import numdifftools as nd

########################## GLOBALS #############################
'''
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
'''
results = []
how_many_uknowns = 0
######################### CONVERSION #############################
def convert_function(fun):
    letter = 'a'
    local_range = range(0 ,how_many_uknowns)
    for ik in local_range:
        fun = fun.replace( chr(ord(letter) + ik), "x[" + str(ik) + "]" )
    fun =  lambda x:fun
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
        #print("self.derivative ", self.derivative)
        result = eval(self.derivative)
        return result
    
    def get_direction(self):
        return self.direction
    
    def get_derivative_fun(self):
        return self.d_fun

########################## FLETCHER-REEVES #########################  
class Fletcher_Reeves:
    def __init__(self, points, e, fun, values, n_order):
        print("sucess")
        points = convert_values(points)
        return_results = []
        self.eta = [None] * how_many_uknowns
        self.points = points
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
        fun =  lambda x:2*x[0]**2 + 3*x[1]**2
        grad2 = nd.Gradient(fun)([self.values[0], self.values[1]])
        for ik in local_range:
            grad = gradient(self.fun, self.values[ik],"x[" + str(ik) + "]" )
            self.grads[ik] = gradient.get_derivative_fun(grad)
            d = -gradient.get_direction(grad)
            self.ds[ik] = d
            self.d_fun.append(gradient.get_derivative_fun(grad))

        range_begin = 0
        range_end = 0
        norm = 10

        range_begin = self.points[0] + 50
        range_end = self.points[1] - 50
        combined_string = self.d_fun[0] + self.d_fun[1]
        x0 = 0.0
        while norm > self.epsylon:
            str_fun = self.create_new_fun(combined_string, letter)
            str_fminfun = self.create_new_fminmsearch(combined_string, letter)
            dich = dichotomy(str_fun, range_begin, range_end)
            min_alpha = dichotomy.get_results(dich)
            min_alpha = np.average(min_alpha, axis = 0)
            min_alpha = 0.00001
            self.calculate_value_in_point(min_alpha)
            norm = math.sqrt(pow(self.ds[0], 2) + pow(self.ds[1], 2)) 
            self.calculate_next_argument(min_alpha)
            self.count_eta()
            self.count_next_d()

        return 1

    def count_next_d(self):
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
            #fun =  lambda x:2*x[0]**2 + 3*x[1]**2
            xk_value = nd.Gradient(self.fun)([self.values[0], self.values[1]])
            xk_plus_1_value = nd.Gradient(self.fun)([self.values[2], self.values[3]])
            xk_value = np.array(xk_value)
            xk_plus_1_value = np.array(xk_plus_1_value)
            xk_product_plus_1 = np.dot(xk_plus_1_value, xk_plus_1_value)
            single_eta = np.divide(xk_product_plus_1, np.dot(xk_value, xk_value))
            self.eta[i] = single_eta

    def calculate_next_argument(self, min_alpha):
        local_range = range(0 ,how_many_uknowns)
        for i in local_range:
            prev_val = self.values[i]
            self.values[i + how_many_uknowns] = prev_val
            next_arg = self.values[i] + min_alpha*self.ds[i]
            self.values[i] = round(next_arg, 2)
    
    def calculate_value_in_point(self, min_alpha):
        local_range = range(0 ,how_many_uknowns)
        letter = 'a'
        for i in local_range:
            function_str = str(self.d_fun[i])
            new_value = self.values[i] + min_alpha*self.ds[i]
            function_str = function_str.replace( chr(ord(letter) + i), str(new_value) )
            function_str = function_str.replace( " ", "" )
            value = compile(function_str, "<string>", "eval")
            value = eval(value)
            self.values[i] = value

    def create_new_fun(self, fun, letter):
        function_str = str(fun)
        local_range = range(0 ,how_many_uknowns)
        for i in local_range:
            if (function_str.find(chr(ord(letter) + i), 0) != -1):
                #print(" create_new_fun.values[i] => ",self.values[i])
                letter_value = str(self.values[i])
                d_value = str(self.ds[i])
                new_value = letter_value + "+" + "x" + "*" + d_value
                fun = str(fun)
                fun = str(fun.replace(chr(ord(letter) + i), new_value))
            else: 
                print('Couldn"t create new function')
                return 0
        return fun

    def create_new_fminmsearch(self, fun, letter):
        function_str = str(fun)
        local_range = range(0 ,how_many_uknowns)
        for i in local_range:
            if (function_str.find(chr(ord(letter) + i), 0) != -1):
                letter_value = str(self.values[i])
                d_value = str(self.ds[i])
                new_value = letter_value + "+" + "x" + "*" + d_value
                fun = str(fun)
                fun = str(fun.replace(chr(ord(letter) + i), new_value))
                fun = str(fun.replace("**", "^"))
            else: 
                print('Couldn"t create new function')
                return 0
        return fun

########################## DICHOTOMY #########################
class dichotomy:
    def __init__(self, fun, rb, re):
        self.function = fun
        self.rb = float(rb)
        self.re = float(re)
        self.minimum = 0
        self.result = self.calculate()

    def calculate(self):
        solutions = [None]*how_many_uknowns
        results_array_len = 0
        tolerance = 0.00001
        range_begin = self.rb
        range_end = self.re
        
            
        distance = np.abs(range_begin-range_end) 
        i = 0
        while (distance >= tolerance):
            delta = distance/ 4
            frb = self.f(range_begin)
            fre = self.f(range_end)    
            #print("frb =>",frb)
            #print("fre =>",fre)
            if range_end < 0 :
                range_end = 0
            if range_begin < 0 :
                range_begin = 0
            if ((range_begin < 0)  and (range_end < 0)) :
                range_begin = range_begin + 50
                range_end = range_end - 50
                continue     
           # print("frb after change =>",frb)
           # print("fre  after change=>",fre)            
            cl = 0.5 * (range_begin + range_end) - delta
            cr = 0.5 * (range_begin + range_end) + delta
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
        #print("range_begin assignment =>",range_begin)
        #print("range_end  assignment=>",range_end)           
        solutions[0] = range_begin
        solutions[1] = range_end
        
        self.minimum = (abs(range_begin) + abs(range_end))/2
        #print("results =>",solutions)
        return solutions

    def f(self, x):
        function_str = str(self.function)
        function_str = function_str.replace( "x", str(x) )
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
points = "1, 2" 
e = "1"
parameters = '1,3'
how_many_uknowns = 2
eng = matlab.engine.start_matlab ()
eng.optimset('Display', 'off');

fun = convert_function(function_str)
parameters = convert_values(parameters)
FR = Fletcher_Reeves(points, float(e), fun, parameters, how_many_uknowns)