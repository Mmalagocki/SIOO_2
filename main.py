from sympy import *
import numpy as np
import math

########################## GLOBALS #############################
x = Symbol('x')
y = Symbol ('y')    
results = []
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
    def __init__(self, function, values):
        self.function = function
        self.derivatives = []
        self.values = values
        self.direction = self.calculate_gradient()

    def calculate_derivative(self):
        function_str = str(self.function)

        if (function_str.find('x' , 0) != -1):
            xprime = self.function.diff(x)
            self.derivatives.append(xprime)
            
            if (function_str.find('y', 0) != -1):
                yprime = self.function.diff(y)
                self.derivatives.append(yprime)
            return 0
    
        else: 
            print('Your function must contain atleast one x parameter')
            return 0

    def calculate_gradient(self):
        vector = []
        self.calculate_derivative()
        for derivative in self.derivatives:
            x_value = self.values[0]
            x_value = str(x_value)
            y_value = self.values[1]
            y_value = str(y_value)            
            derivative = str(derivative)
            derivative = derivative.replace( " ", "" )
            derivative = derivative.replace("x", x_value)
            derivative = derivative.replace("y", y_value)
            vector.append(eval(derivative))
        return vector
    
    def get_direction(self):
        return self.direction

########################## FLETCHER-REEVES #########################  
class Fletcher_Reeves:
    def __init__(self, points, e, fun, values):
        points = convert_values(points)
        print(points)
        self.points = points
        self.epsylon = e
        self.k = 1    
        self.fun = fun
        self.values = values
        self.calculate()
        
    def calculate(self):
        grad = gradient(self.fun, self.values)
        d = gradient.get_direction(grad)
        range_begin = []
        range_end = []

        for point in self.points:
            range_begin.append(point + 50)
            range_end.append(point - 50)
        print(range_begin[0], "   ", range_end[0])
        dich = dichotomy(self.fun, range_begin, range_end)
        point_value = self.calculate_value_in_point()
        min_alpha = dich.get_minimum()
        
        norm = math.sqrt(pow(d[o], 2) + pow(d[1], 2))
        if (norm< self.epsylon):
            self.return_results()
        dk = -d[0]  
        grad2 = gradient(self.fun, self.values)
        self.point = point_value + min_alpha*dk
        eta = (d[0])
        dk_1    
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
        self.calculate()
        

    def calculate(self):
        results = []
        results_array_len = len(results)
        tolerance = 0.00001
        ranges_begins = self.rb
        ranges_ends = self.re
        distance = np.abs(range_begin[0]-range_end[0]) 
        i = 0

        for range_begin in ranges_begining:
            range_end = range_ends[i]
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
    
            return  range_begin, range_end    
    
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

    def f(self, x, y = False):
        function_str = str(self.function)
        function_str = function_str.replace( "x", str(x) )
        function_str = function_str.replace( " ", "" )
        code = compile(function_str, "<string>", "eval")
        return eval(code)
    
    def get_minimum(self):
        return self.minimum


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
function_str = "2*x**3"
points = "19, 20" 
e = "10"
parameters = '20'


fun = convert_function(function_str)
parameters = convert_values(parameters)
FR = Fletcher_Reeves(points, float(e), (fun), parameters)


