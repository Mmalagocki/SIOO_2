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
#print(type(i))
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
        self.direction = self.calculate_gradient()

    def calculate_derivative(self):
        function_str = str(self.function)
        if (function_str.find(self.letter , 0) != -1):
            xprime = self.function.diff(self.letter)
            #print(xprime)
            #self.derivatives.append(xprime)
            self.derivative = xprime
            return 0
    
        else: 
            print('Your function must contain atleast one x parameter')
            return 0

    def calculate_gradient(self):
        result = 0
        self.calculate_derivative()
        #for derivative in self.derivatives:
        letter_value = self.value
        letter_value = str(letter_value)
        self.derivative = str(self.derivative)
        self.derivative = self.derivative.replace( " ", "" )
        self.derivative = self.derivative.replace(self.letter, letter_value)
        result = eval(self.derivative)
        return result
    
    def get_direction(self):
        #print(self.direction)
        return self.direction

########################## FLETCHER-REEVES #########################  
class Fletcher_Reeves:
    def __init__(self, points, e, fun, values, n_order):
        points = convert_values(points)
        #print(points)
        self.points = points
        self.epsylon = e
        self.k = 1    
        self.fun = fun
        self.values = values
        self.ds = []
        self.calculate()
        
    def calculate(self):
        letter = 'a'
        local_range = range(0 ,how_many_uknowns)
        for i in local_range:
            grad = gradient(self.fun, self.values[i], chr(ord(letter) + i))
            d = gradient.get_direction(grad)
            self.ds.append(d)

        #print(self.ds)
        range_begin = 0
        range_end = 0

        for point in self.points:
            range_begin = point + 50
            range_end = point - 50
        #print(range_begin[0], "   ", range_end[0])
        #print('passed fun' , self.fun)
        self.create_new_fun(letter)
        dich = dichotomy(self.fun, range_begin, range_end)
        results = dichotomy.get_results(dich)
        print(results)
        '''
        point_value = self.calculate_value_in_point()
        min_alpha = dich.get_minimum()
        norm = math.sqrt(pow(self.ds[0], 2) + pow(self.ds[1], 2))
        if (norm< self.epsylon):
            self.return_results()
        dk = -self.ds[0]

        grad2 = gradient(self.fun, self.values)
        #self.point = point_value + min_alpha*dk
        eta = (self.ds[0])
        dk_1
        '''
        '''
        print(point_value)
        print(dk)
        print(min_alpha)
        '''
        
        return 1
        
    def calculate_value_in_point(self):
        function_str = str(self.fun)
        function_str = function_str.replace( "x", str(self.point) )
        function_str = function_str.replace( " ", "" )
        value = compile(function_str, "<string>", "eval")
        return eval(value)

    def create_new_fun(self, letter):
        function_str = str(self.fun)
        local_range = range(0 ,how_many_uknowns)
        for i in local_range:
            if (function_str.find(chr(ord(letter) + i), 0) != -1):
                letter_value = str(self.values[i])
                d_value = str(self.ds[i])
                new_value = letter_value + "+" + "alpha" + "*" + d_value
                self.fun = str(self.fun)
                self.fun = str(self.fun.replace(chr(ord(letter) + i), new_value))
                #print(self.fun)
            else: 
                print('Couldn"t create new function')
                return 0
        
    def return_results(self, return_values):
        #unused
        return e

########################## DICHOTOMY #########################
class dichotomy:
    def __init__(self, fun, rb, re):
        self.function = fun
        self.rb = rb
        self.re = re
        self.minimum = 0
        self.result = self.calculate()



    def calculate(self):
        results = []
        results_array_len = len(results)
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
                #print ("1")
                range_begin = cl
            elif ((frb >= fcl) and (frb >= fre) and (fcl >= fcr)) :
                #print ("2")
                range_begin =  cl
            elif ((frb <= fcl) and (frb <= fre) and (fcl >= fcr)):
                #print ("3")
                range_begin = cl    
            elif ((frb >= fcl) and (frb >= fre) and (fcl >= fcr)):
                #print ("4")
                range_begin = cl
            elif ((frb >= fcl) and (frb <= fre) and (fcl <= fcr)):
                #print ("5")
                range_end = cr
            elif ((frb <= fcl) and (frb <= fre) and (fcl <= fcr)):
                #print ("6")
                range_end = cr                    
            elif (fre < fcl) and (fcr> frb):
                #print ("7")
                self.not_unimodal(range_begin, cl)
                self.not_unimodal(cr, range_end)
                self.not_unimodal(cl, cr)
            elif ((frb >= fcl) and (frb < fre) and (fcl <= fcr)):
                #print ("8")
                range_begin = cl              
            elif ((frb >= fcl) and (frb < fre) and (fcl > fcr)):
                #print ("9")
                range_end = cr                    
            distance = np.abs(range_begin-range_end)
        
        results_array_len += 1
        results.insert(results_array_len, range_begin)
        results_array_len += 1
        results.insert(results_array_len, range_end)
        
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
                #print ("1")
                range_begin = cl
            elif ((frb >= fcl) and (frb >= fre) and (fcl >= fcr)) :
                #print ("2")
                range_begin =  cl
            elif ((frb <= fcl) and (frb <= fre) and (fcl >= fcr)):
                #print ("3")
                range_begin = cl    
            elif ((frb >= fcl) and (frb >= fre) and (fcl >= fcr)):
                #print ("4")
                range_begin = cl
            elif ((frb >= fcl) and (frb <= fre) and (fcl <= fcr)):
                #print ("5")
                range_end = cr
            elif ((frb <= fcl) and (frb <= fre) and (fcl <= fcr)):
                #print ("6")
                range_end = cr                    
            elif (fre < fcl) and (fcr> frb):
                #print ("7")
                self.not_unimodal(range_begin, cl)
                self.not_unimodal(cr, range_end)
                self.not_unimodal(cl, cr)
            elif ((frb >= fcl) and (frb < fre) and (fcl <= fcr)):
                #print ("8")
                range_begin = cl              
            elif ((frb >= fcl) and (frb < fre) and (fcl > fcr)):
                #print ("9")
                range_end = cr                    
    
        distance = np.abs(range_begin-range_end)
        results_array_len += 1
        results.insert(results_array_len, range_begin)
        results_array_len += 1
        results.insert(results_array_len, range_end)        

    def f(self, x):
        function_str = str(self.function)
        function_str = function_str.replace( "alpha", str(x) )
        function_str = function_str.replace( " ", "" )
        code = compile(function_str, "<string>", "eval")
        return eval(code)

    def get_minimum(self):
        return self.minimum

    def get_results(self):
        print(self.result)
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
function_str = "2*a**3 + 3*b"
points = "19, 20" 
e = "10"
parameters = '20,30'
how_many_uknowns = 2


fun = convert_function(function_str)
#print(fun.diff(a))
parameters = convert_values(parameters)
FR = Fletcher_Reeves(points, float(e), (fun), parameters, how_many_uknowns)


