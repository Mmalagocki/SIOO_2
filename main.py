from sympy import *
import inspect
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
a = Symbol('a')
b = Symbol ('b')
results = []
how_many_uknowns = 0
######################### CONVERSION #############################
def convert_function(fun):
    letter = 'a'
    local_range = range(0 ,how_many_uknowns)
    for i in local_range:
        fun = fun.replace( chr(ord(letter) + i), "x[" + str(i) + "]" )
    return fun

def convert_values(values_string):
    values_array = []
    values_string.replace( " ", "" )
    separator = values_string.find("," , 0)
    values_array.append(float(values_string[0:separator]))
    values_array.append(float(values_string[separator+1:]))
    return values_array

def get_derivatives(func):
    arg_symbols = symbols(inspect.getargspec(func).args)
    sym_func = func(*arg_symbols)

    return [lambdify(arg_symbols, sym_func.diff(a)) for a in arg_symbols]

########################## DERIVATIVES #########################
class gradient:
    def __init__(self, function):
        function = sympify(function, evaluate = False)
        self.function = function
        self.d_fun = []
        self.calculate_derivative()

    def calculate_derivative(self):
        letter = 'a'
        function_str = str(self.function)
        local_range = range(0 ,how_many_uknowns)
        for i in local_range:
            current_letter = chr(ord(letter) + i)
            if (function_str.find(current_letter, 0) != -1):
                xprime = self.function.diff(current_letter)
                self.d_fun.append(str(xprime))
            else: 
                print("Something went wrong")
    
    def get_derivative_fun(self):
        return self.d_fun

########################## FLETCHER-REEVES #########################  
class Fletcher_Reeves:
    def __init__(self, points, e, fun, values, n_order, function_str):
        points = convert_values(points)
        return_results = []
        self.eta = 0
        self.points = points
        self.epsylon = e
        self.k = 1
        self.fun = fun
        self.derivative = []
        self.values = values
        self.function_str = function_str

        local_range = range(0 ,how_many_uknowns)
        for ik in local_range:
            self.values.append(None)
        self.d = [None]*how_many_uknowns
        self.calculate()
        
    def calculate(self):
        local_range = range(0 ,how_many_uknowns)
        fun =  lambda x: eval(self.fun)
        d = -nd.Gradient(fun)([self.values[0], self.values[1]])
        self.d = d
        grad = gradient(self.function_str)
        self.derivative = gradient.get_derivative_fun(grad)
        local_range = range(0 ,how_many_uknowns)
        for i in local_range:        
            self.derivative[i] = convert_function(self.derivative[i])
        range_begin = 0
        range_end = 0
        norm = 10

        range_begin = self.points[0] + 50
        range_end = self.points[1] - 50
        while norm > self.epsylon:
            str_vals = self.find_values_for_min(self.derivative)
            #str_fminfun = self.create_new_fminmsearch(self.derivative)
            str_fun = self.createn_new_fun(self.function_str ,str_vals)
            dich = dichotomy(str_fun, range_begin, range_end)
            min_alpha = dichotomy.get_results(dich)
            min_alpha = np.average(min_alpha, axis = 0)
            self.calculate_new_value(min_alpha)
            norm = math.sqrt(pow(self.d[0], 2) + pow(self.d[1], 2)) 
            print(norm)
            self.count_eta(fun)
            self.count_next_d(fun)

        return 1

    def createn_new_fun(self, fun, vals):
        letter = 'a'
        local_range = range(0 ,how_many_uknowns)
        for i in local_range:
            fun = fun.replace( chr(ord(letter) + i), vals[i] )
        return fun        
        
        
    def count_next_d(self, fun):
        local_range = range(0 ,how_many_uknowns)
        letter = 'a'
        xk_plus_1_value = nd.Gradient(fun)([self.values[2], self.values[3]])
        next_d = -xk_plus_1_value + self.eta * self.d
        self.d = next_d

    def count_eta(self, fun):
        local_range = range(0 ,how_many_uknowns)
        for i in local_range:
            xk_value = nd.Gradient(fun)([self.values[0], self.values[1]])
            xk_value = np.array(xk_value)
            xk_value_product = np.dot(xk_value, xk_value)
            
            xk_plus_1_value = nd.Gradient(fun)([self.values[2], self.values[3]])
            xk_plus_1_value = np.array(xk_plus_1_value)
            xk_product_plus_1 = np.dot(xk_plus_1_value, xk_plus_1_value)
            
            single_eta = np.divide(xk_product_plus_1, xk_value_product)
            self.eta = single_eta

    def calculate_next_argument(self, min_alpha):
        local_range = range(0 ,how_many_uknowns)
        for i in local_range:
            prev_val = self.values[i]
            self.values[i + how_many_uknowns] = prev_val
            next_arg = self.values[i] + min_alpha*self.d[i]
            self.values[i] = round(next_arg, 2)
    
    def calculate_new_value(self, min_alpha):
        local_range = range(0 ,how_many_uknowns)
        function_str = self.derivative
        for i in local_range:
            self.values[i + how_many_uknowns] = self.values[i]
            value = self.values[i] + min_alpha*self.d[i]
            self.values[i] = value

    def find_values_for_min(self, fun):
        function_str = str(fun)
        results = [0] * how_many_uknowns
        local_range = range(0 ,how_many_uknowns)
        for i in local_range:
            if (function_str.find("x[" + str(i) + "]", 0) != -1):
                letter_value = str(self.values[i])
                d_value = str(self.d[i])
                results[i] = letter_value + "+" + "alpha" + "*" + d_value
            else: 
                print('Couldn"t create new function')
                return 0
        return results

    def create_new_fminmsearch(self, fun):
        function_str = str(fun)
        local_range = range(0 ,how_many_uknowns)
        for i in local_range:
            if (function_str.find("x[" + str(i) + "]", 0) != -1):
                letter_value = str(self.values[i])
                d_value = str(self.d[i])
                new_value = letter_value + "+" + "x" + "*" + d_value
                fun = str(fun)
                fun = str(fun.replace("x[" + str(i) + "]", new_value))
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
points = "1, 2" 
e = "1"
values = '1,3'
how_many_uknowns = 2
eng = matlab.engine.start_matlab ()
eng.optimset('Display', 'off');

fun = convert_function(function_str)
values = convert_values(values)
FR = Fletcher_Reeves(points, float(e), fun, values, how_many_uknowns, function_str)