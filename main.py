from sympy import *
import numpy as np

########################## GLOBALS #############################
x = Symbol('x')
y = Symbol ('y')    
########################## DERIVATIVES #########################
class derivative:
    def __init__(self, function):
        self.function = function
        self.calculate_derivative()
        
    def calculate_derivative(self):
        function_str = str(self.function)
        #print(function_str)
        if (function_str.find('x' , 0) != -1):
            xprime = self.function.diff(x)
            print(xprime)
            if (function_str.find('y', 0) != -1):
                yprime = self.function.diff(y)
                print(yprime)      
        else: 
            print('Your function must contain atleast one x parameter')
                
                #function = 3*x**2 + 1 - 1*x + 20*y**4            
                #f = lambdify(x, yprime, 'numpy')
                #print(f(np.ones(5)))

def convert_function(fun):
    '''
    parameters_places = []
    parameters_places_array_length = 0
    last_found = 0
    
    while (fun.find("x" , last_found + 1) != -1):
        print('Im looking for x')
        found = (fun.find("x" , last_found))
        parameters_places.insert(parameters_places_array_length, found)
        last_found = found + 1
        parameters_places_array_length += 1
        #fun = fun.replace( "x", "1" )
        print(fun.replace( "x", "1" ))
     
    last_found = 0
    
    while (fun.find("y" , last_found) != -1):
        print('Im looking for y')
        found = (fun.find("y" , last_found))
        parameters_places.insert(parameters_places_array_length, found)
        last_found = found + 1
        parameters_places_array_length += 1
        fun.replace( "y", "1" )
    '''   

    #print('Before sympify => ',fun)   
    fun = sympify(fun, evaluate = False)
    #print(fun)
    return fun


########################## FLETCHER-REEVES #########################
                
class fletcher_Reeves:
    def __init__(self, point, e):
        self.derivative = []
        self.point = point
        self.epsylon = e
        self.calculate()
        
    def calculate():
        if(abs(grad) < e):
            self.results()
        xi = xi_1 +labmdai*si
        
            
    def results(e):
        #unused
        return e
    
    
    

###################### MAIN #################
#function = 3*y**2 + 1 - 1*y
print("Pass the function")
function_str = input()
print("Pass the point")
point = input()
fun = convert_function(function_str)
#function_str = function_str.replace( "x", x )
#function_str = function_str.replace( "y", y )
#function = float(function_str)


derivatives = derivative(fun)
