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
eng = matlab.engine.start_matlab ()
eng.optimset('Display', 'off');

########################## GUI #############################
def iterations_chosen():
    Label(frame, text="You have chosen accuracy").grid(row = 0)
    
    Label(frame, text="Function:").grid(row = 1)
    function_input = Entry(frame, width = 20, cursor = 'hand2')
    function_input.insert(0,'')
    function_input.grid(row=1 , column=1, pady = 10)
    
    Label(frame, text="Iterations:").grid(row = 2)
    iterations_input = Entry(frame, width = 20, cursor = 'hand2')
    iterations_input.insert(0,'')
    iterations_input.grid(row = 2 , column=1, pady = 10)
    
    Label(frame, text="Pass the values (2,2):").grid(row = 3)
    values_input = Entry(frame, width = 20, cursor = 'hand2')
    values_input.insert(0,'')
    values_input.grid(row = 3 , column = 1, pady = 10)
    
    Label(frame, text="How many unkowns:").grid(row = 4)
    how_many_uknowns_input = Entry(frame, width = 20, cursor = 'hand2')
    how_many_uknowns_input.insert(0,'')
    how_many_uknowns_input.grid(row = 4 , column = 1, pady = 10)    
     
    Label(frame, text="Pass the points:").grid(row = 5)
    points_input = Entry(frame, width = 20, cursor = 'hand2')
    points_input.insert(0,'')
    points_input.grid(row = 5 , column = 1, pady = 10)    
         
    
    Button_submit = Button(frame, text = "Submit", command = lambda: set_iterations_values(function_input.get(), 
                                                                               iterations_input.get(), 
                                                                               values_input.get(),
                                                                               how_many_uknowns_input.get(),
                                                                               points_input.get()
                                                                               ))
    Button_submit.grid(row=6 , column=1)

def accuracy_chosen():
    Label(frame, text="You have chosen accuracy").grid(row = 0)
    
    Label(frame, text="Function:").grid(row = 1)
    function_input = Entry(frame, width = 20, cursor = 'hand2')
    function_input.insert(0,'')
    function_input.grid(row=1 , column=1, pady = 10)
    
    Label(frame, text="Accuracy:").grid(row = 2)
    accuracy_input = Entry(frame, width = 20, cursor = 'hand2')
    accuracy_input.insert(0,'')
    accuracy_input.grid(row = 2 , column=1, pady = 10)
    
    Label(frame, text="Pass the values (2,2):").grid(row = 3)
    values_input = Entry(frame, width = 20, cursor = 'hand2')
    values_input.insert(0,'')
    values_input.grid(row = 3 , column = 1, pady = 10)
    
    Label(frame, text="How many unkowns:").grid(row = 4)
    how_many_uknowns_input = Entry(frame, width = 20, cursor = 'hand2')
    how_many_uknowns_input.insert(0,'')
    how_many_uknowns_input.grid(row = 4 , column = 1, pady = 10)    
     
    Label(frame, text="Pass the points:").grid(row = 5)
    points_input = Entry(frame, width = 20, cursor = 'hand2')
    points_input.insert(0,'')
    points_input.grid(row = 5 , column = 1, pady = 10)    
         
    
    Button_submit = Button(frame, text = "Submit", command = lambda: set_accuracy_values(function_input.get(), 
                                                                               accuracy_input.get(), 
                                                                               values_input.get(),
                                                                               how_many_uknowns_input.get(),
                                                                               points_input.get()
                                                                               ))
    Button_submit.grid(row=6 , column=1)

def choose_main_condition(chosen_condition):
    if (chosen_condition == "Get the result based on  accuracy") :
        accuracy_chosen()
    elif(chosen_condition == "Get the result based on  number of iterations"):
        iterations_chosen()
    else:
        something_went_wront(chosen_condition)

def something_went_wront(chosen_condition):
    Label(frame, text="Sorry! Something went wrong. Here is the codition: " + chosen_condition).grid(row=0)

def set_accuracy_values(function, accuracy, values, hmu, points):
    set_accuracy(accuracy)
    set_function(function)
    set_iterations(1000000)
    set_values(values)
    set_hmu(hmu)
    set_points(points)
    fun = convert_function(function_str)
    values = convert_values(values)
    FR = Fletcher_Reeves(points, float(e), fun, values, how_many_uknowns, function_str)


def set_iterations_values(function, interations, values, hmu, points):
    set_accuracy('0.00001')
    set_iterations(interations)
    set_function(function)
    set_values(values)
    set_hmu(hmu)
    set_points(points)
    fun = convert_function(function_str)
    values = convert_values(values)
    FR = Fletcher_Reeves(points, float(e), fun, values, how_many_uknowns, function_str)

def set_points(input_string):
    global points
    points = input_string

def set_hmu(input_string):
    global how_many_uknowns
    how_many_uknowns = int(input_string)
    
def set_accuracy(input_string):
    global e
    e = input_string

def set_values(input_string):
    global values
    values = input_string
    
def set_function(input_string):
    global function_str
    function_str = input_string

def set_iterations(input_string):
    global iterations_limit
    iterations_limit = float(input_string)

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
        #### Initialization ####
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
        #### Step 2 ####

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

        range_begin = self.points[0] - 50
        range_end = self.points[1] + 50
        iterations = 0
        function_values = []
        xopt = [None]*how_many_uknowns
        lowest_norm = 0
        iterations_since_finding_min = 0

        #### Step 3 ####

        #while norm > self.epsylon and iterations < iterations_limit:
        while norm > self.epsylon:
            str_vals = self.find_values_for_min(self.derivative)
            str_fun = self.createn_new_fun(self.function_str ,str_vals)
            print("str_fun", str_fun)
            dich = dichotomy(str_fun, range_begin, range_end)
            min_alpha = dichotomy.get_results(dich)
            min_alpha = np.average(min_alpha, axis = 0)
            function = lambda alpha: eval(str_fun)
            xopt = opt.fminbound(function, self.points[0],self.points[1])
            print("   min_alpha result =>", min_alpha)
            print("   xopt result =>", xopt)
            our_result = self.function_value(min_alpha, str_fun)

            #### Step 4 ####

            self.calculate_new_value(min_alpha)
            norm = LA.norm(self.d)
            print("Current nom ==>", norm)

            #### Step 5 ####

            self.count_eta(fun)
            self.count_next_d(fun)
            iterations = iterations + 1
            if iterations == 1:
                lowest_norm = norm + 1
            function_values.append(norm)

            if lowest_norm > norm:
                print("****************************** I FOUND MINIMUM ******************************")
                lowest_norm = norm
                iterations_since_finding_min = 0
            print("     lowest_norm ==>", lowest_norm)

            print(iterations)

        #### Drawing Plot ####
        local_range = range(0, math.floor(iterations))
        x_array = []
        for i in local_range:
            x_array.append(i)
        plt.plot(x_array,function_values, 'ro')
        plt.grid()  
        plt.show()
        return lowest_norm

    '''
    Function used to calculate value of function in point

    float: min_alpha - value local minimum found with optimizing function
    string: function - string of fuction with all values set except alpha
    '''
    def function_value (self, min_alpha, function):
        function_str = str(function)
        function_str = function_str.replace( "alpha", str(min_alpha))
        result = eval(function_str)
        return result


    '''
    Function used to swap unknowns with xk + alhpa * d form

    array of strings: vals - value local minimum found with optimizing function
    string: fun - string of fuction
    '''
    def createn_new_fun(self, fun, vals):
        letter = 'a'
        local_range = range(0 ,how_many_uknowns)
        for i in local_range:
            fun = fun.replace( chr(ord(letter) + i), vals[i] )
        return fun

    '''
    Function used to find new direction of minimalization

    string: fun - string of fuction
    '''
    def count_next_d(self, fun):
        local_range = range(0 ,how_many_uknowns)
        letter = 'a'
        xk_plus_1_value = nd.Gradient(fun)([self.values[0], self.values[1]])
        next_d = -xk_plus_1_value + self.eta * self.d
        self.d = next_d

    '''
    Function used to calculate eta

    string: fun - string of fuction
    '''
    def count_eta(self, fun):
        if how_many_uknowns != 1:
            local_range = range(0 ,how_many_uknowns)
            for i in local_range:
                xk_value = nd.Gradient(fun)([self.values[2], self.values[3]])
                xk_value = np.array(xk_value)
                xk_value_product = np.dot(xk_value, xk_value)
                
                xk_plus_1_value = nd.Gradient(fun)([self.values[0], self.values[1]])         
                xk_plus_1_value = np.array(xk_plus_1_value)
                xk_product_plus_1 = np.dot(xk_plus_1_value, xk_plus_1_value)
                
                single_eta = np.divide(xk_product_plus_1, xk_value_product)
                self.eta = single_eta
        else:
            xk_value = nd.Gradient(fun)(self.values[1])
            xk_value = np.array(xk_value)
            xk_value_product = np.dot(xk_value, xk_value)
            
            xk_plus_1_value = nd.Gradient(fun)(self.values[0])         
            xk_plus_1_value = np.array(xk_plus_1_value)
            xk_product_plus_1 = np.dot(xk_plus_1_value, xk_plus_1_value)
            
            single_eta = np.divide(xk_product_plus_1, xk_value_product)
            self.eta = single_eta            

    '''
    Function used to calculate value of function in point

    float: min_alpha - value local minimum found with optimizing function
    '''
    def calculate_new_value(self, min_alpha):
        local_range = range(0 ,how_many_uknowns)
        function_str = self.derivative
        for i in local_range:
            self.values[i + how_many_uknowns] = self.values[i]
            value = self.values[i] + min_alpha*self.d[i]
            self.values[i] = value

    '''
    Function used to swap unknowns with xk + alhpa * d form

    float: min_alpha - value local minimum found with optimizing function
    '''
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
        tolerance = 0.000000000001
        range_begin = self.rb
        range_end = self.re
          
        distance = np.abs(range_begin-range_end) 
        i = 0
        while (distance >= tolerance):
            delta = distance/ 4
            frb = self.f(range_begin)
            fre = self.f(range_end)
            if range_end < 0 :
                range_end = 0
            if range_begin < 0 :
                range_begin = 0
            if ((range_begin < 0)  and (range_end < 0)) :
                range_begin = range_begin + 50
                range_end = range_end - 50
                continue
            if frb < 0:
                frb = 0
                fre = abs(fre)
            if fre < 0:
                frb = abs(frb)
                fre = 0
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
        solutions[0] = range_begin
        solutions[1] = range_end
        
        self.minimum = (abs(range_begin) + abs(range_end))/2
        return solutions

    def f(self, x):
        function_str = str(self.function)
        function_str = function_str.replace( "alpha", str(x) )
        function_str = function_str.replace( " ", "" )
        code = compile(function_str, "<string>", "eval")
        return eval(code)

    def get_results(self):
        return self.result


########################## MAIN #########################
root = Tk()
root.geometry("1200x900")
root.title("MM&MJ")

 
frame = Frame(root)
B = Button(root, text = "Choose condition", command = lambda: choose_main_condition(tkvarq.get()) )

options = ["Get the result based on  accuracy",
           "Get the result based on  number of iterations"
           ]

## SELECT MENU
tkvarq = StringVar(root)
tkvarq.set(options[1])
question_menu = OptionMenu(root, tkvarq, *options)
question_menu.pack()
B.pack()
frame.pack()
### DISPLAYS CHOSEN VERSION
root.mainloop()
'''
function_str = "4*a**2 + 4*b**2"
points = "-2, 2" 
e = "0.01"
parameters = '3,5'
how_many_uknowns = 2
eng = matlab.engine.start_matlab ()
eng.optimset('Display', 'off');

fun = convert_function(function_str)
parameters = convert_values(parameters)
FR = Fletcher_Reeves(points, float(e), fun, parameters, how_many_uknowns, function_str)
'''