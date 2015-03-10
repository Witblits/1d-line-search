#!/usr/bin/python
#This script is based on the Particle Swarming algorithm which is used to approximate the discrete global optimum.

from __future__ import division				#ensures division always returns the correct answer (eg. 1/2 = 0.5 instead of 1/2 = 0)
import argparse
import pylab as pl
from sympy import *

###GLOBAL VARIABLES
##########################################################################################
#total_runs = 0
#pop = 0                  		#population size
#n_dim = 0                 		#number of dimensions

a = -1			#line search start position
b = 1			#line search end position
epsilon1 = 0.01		#prescribed accuracy
#i = 0			#iteration number
L = b - a		#Length of line search is conducted on
r = 0.618034		#golden ratio
algorithm_string = ""
example_string = ""
graph = []
graph1 = []
graph2 = []
x_initial = 0
##########################################################################################


###PARSE INFORMATION
##########################################################################################
parser = argparse.ArgumentParser(description="Line Search Algorithms",formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("algorithm", help="Select algorithm type:\n"+
		    "1 - 1-D line search based on the Newton-Raphson Method\n"+
		    "2 - Modification of the Newton-Raphson we discussed in class\n"+
		    "3 - Golden Section Method\n\n", type=int, choices=[1,2,3])					#add a new POSITIONAL argument with a help message included

parser.add_argument("choice", help="Select example function to test the Line Search Algorithms on:\n"+
		    "1 - Test 1\n"+
		    "2 - Test 2\n"+
		    "3 - Test 3\n"+
		    "4 - Test 4\n\n", type=int, choices=[1,2,3,4])					#add a new POSITIONAL argument with a help message included

#parser.add_argument("i", help="Enter total number of iterations the algorithm should complete", type=int)	#add a new POSITIONAL argument with a help message included
#parser.add_argument("p", help="Enter desired population the algorithm should use", type=int)			#add a new POSITIONAL argument with a help message included
#parser.add_argument("n", help="Enter number of dimensions the example problem should have", type=int)		#add a new POSITIONAL argument with a help message included
#parser.add.argument("filename", help="Enter the name of the file the output figure should be saved to", type=argparse.FileType('w'))

#parser.add_argument("-a","--average", default=1, help="Enter total number of independant runs to calculate average solution", type=int)	#add a new OPTIONAL argument with a help message included

args = parser.parse_args()													#store the argument parsed into the variable "args"

#if args.i != None:
  #total_runs = args.i
  
#if args.p != None:
  #pop = args.p
  
#if args.n != None:
  #n_dim = args.n
##########################################################################################


###OBTAIN TEST_PARAMETERS
##########################################################################################

#1 Test function 1
if args.choice == 1:
  a = 0
  b = 2
  epsilon1 = 0.01
  example_string = "Test function 1"

#2 Test function 2
elif args.choice == 2:
  a = 0
  b = (pl.np.pi)/2
  epsilon1 = 0.001
  example_string = "Test function 2"

#3 Test function 3
elif args.choice == 3:
  a = -1.9
  b = 0.9
  function_evaluations = 10
  example_string = "Test function 3"

#4 Test function 4
elif args.choice == 4:
  a = 0
  b = 20
  epsilon1 = 1e-5
  example_string = "Test function 4"
  
  
L = b - a		#Length of line search is conducted on
##########################################################################################


###OPTIMIZATION METHODS
##########################################################################################
### 1-D line search based on the Newton-Raphson Method
def lineSearchNewtonRapshon(a,b):
  solution = pl.np.zeros((3,1))
  solution_found = False
  i = 0
  x_old = (b+a)/2 		#assume a solution nearby (works good for all test functions except Test Function 4 which does not converge)
  #x_old = 0.9*b			#works for all test functions (just requires more iterations for some)
  x = x_old
  global x_initial
  x_initial = x_old
  #print "x_old: %.3f" %x_old
  #print "x    : %.3f" %x
  
  while (pl.np.abs(x-x_old) > epsilon1) or (i == 0):
    x_old = x
    x = x_old - (testFunction(1,x_old)/testFunction(2,x_old))
    i = i + 1
    #print "pl.np.abs(x-x_old) : %.3f" %pl.np.abs(x-x_old)
    #print "i	: %d" %i
    #print "x_old: %.3f" %x_old
    #print "x    : %.3f" %x
    #print "f(x) : %.3f" %testFunction(0,x_old)
  
  solution[0] = x
  solution[1] = testFunction(0,x)
  solution[2] = i
  solution_found = True 
  return solution

### Modification of the Newton-Raphson we discussed in class
def modifiedNewtonRapshon(a,b,L,r):
  solution = pl.np.zeros((3,1))
  solution_found = False
  epsilon1 = 1e-5
  i = 0
  #x_old = (b+a)/2 		#assume a solution nearby (works good for all test functions except Test Function 4 which does not converge)
  x_old = 0.9*b			#works for all test functions (just requires more iterations for some)
  x = x_old
  global x_initial
  x_initial = x_old
  #print "x_old: %.3f" %x_old
  #print "x    : %.3f" %x
  
  while (pl.np.abs(x-x_old) > epsilon1) or (i == 0):
    x_old = x
    numerator = testFunction(1,x_old)
    denominator = maximize_f2ndDerivative(x_old,epsilon1,r,L)
    x = x_old - numerator/denominator		#numerator/denominator
    i = i + 1
    #print "pl.np.abs(x-x_old) : %.3f" %pl.np.abs(x-x_old)
    #print "i	: %d" %i
    #print "x_old: %.3f" %x_old
    #print "x    : %.3f" %x
    #print "f(x) : %.3f" %testFunction(0,x_old)
  
  solution[0] = x
  solution[1] = testFunction(0,x)
  solution[2] = i
  solution_found = True 
  return solution  
  
  
  return solution

### Golden Section Method
def goldenSection(a,b,L,r):
  solution = pl.np.zeros((3,1))
  i = 0
  lambda1 = a + pl.np.power(r,2)*L
  lambda2 = a + r*L
  solution_found = False
  order = 0
  
  while(solution_found == False):
    f1 = testFunction(order,lambda1)
    f2 = testFunction(order,lambda2)
    i = i + 1
    
    if(f1 > f2):
      a = lambda1
      lambda1 = lambda2
      L = b - a
      lambda2 = a + r*L
    elif (f1 < f2):
      b = lambda2
      lambda2 = lambda1
      L = b - a
      lambda1 = a + pl.np.power(r,2)*L
    
    if L < epsilon1 or (i == 10 and args.choice == 3):
      solution[0] = (b+a)/2
      solution[1] = testFunction(order,solution[0])
      solution[2] = i
      solution_found = True
      return solution
##########################################################################################

###ADDITIONAL HELPER FUNCTION
def maximize_f2ndDerivative(x_old,epsilon,r,L):
  #use goldenSection method to determine maximum point on 2nd derivative function within range
  #goldenSection is selected instead of Newton-Raphsody method because we assume that we only have
  #a twice differentiable (smooth) f(x) available. Thus we cant apply Newton-Raphsody to a function
  #which is already the 2nd derivative.
  a = x_old - epsilon
  b = x_old + epsilon
  L = b - a
  solution = pl.np.zeros((3,1))
  i = 0
  lambda1 = a + pl.np.power(r,2)*L
  lambda2 = a + r*L
  solution_found = False
  order = 2			#this is changed to 2 because we are already working with the 2nd order function
  
  while(solution_found == False):
    f1 = testFunction(order,lambda1)
    f2 = testFunction(order,lambda2)
    i = i + 1
    
    if(f1 > f2):
      a = lambda1
      lambda1 = lambda2
      L = b - a
      lambda2 = a + r*L
    elif (f1 < f2):
      b = lambda2
      lambda2 = lambda1
      L = b - a
      lambda1 = a + pl.np.power(r,2)*L
    
    if L < epsilon or (i == 10 and args.choice == 3):
      return testFunction(order,(b+a)/2)


###SYMBOLIC TEST_FUNCTIONS f(x), f'(x), f''(x)
##########################################################################################
def test1(order,lambdax):
  x = Symbol('x')	#declare symbolic variable
  f = Function('f')(x)
  
  f = x**2 + 2*exp(-x)
  
  if order == 0:
    return f.subs(x,lambdax)
  else:
    f = f.diff(x,order)
    return f.subs(x,lambdax)
    
def test2(order,lambdax):
  x = Symbol('x')	#declare symbolic variable
  f = Function('f')(x)
  
  f = x*cos(x)
  
  if order == 0:
    return f.subs(x,lambdax)
  else:
    f = f.diff(x,order)
    return f.subs(x,lambdax)

def test3(order,lambdax):
  x = Symbol('x')	#declare symbolic variable
  f = Function('f')(x)
  
  f = 4*(x-7)/(x**2 + x - 2)
  
  if order == 0:
    return f.subs(x,lambdax)
  else:
    f = f.diff(x,order)
    return f.subs(x,lambdax)
  
def test4(order,lambdax):
  x = Symbol('x')	#declare symbolic variable
  f = Function('f')(x)
  
  f = x**4 - 20*x**3 + 0.1*x
  
  if order == 0:
    return f.subs(x,lambdax)
  else:
    f = f.diff(x,order)
    return f.subs(x,lambdax) 

'''
###NUMERICAL TEST_FUNCTIONS f(x)
##########################################################################################
def test1(lambdax):
  return pl.np.power(lambdax,2) + 2*pl.np.exp(-lambdax)

def test2(lambdax):
  return -lambdax*pl.np.cos(lambdax)	#maximize

def test3(lambdax):
  return 4*(lambdax-7)/(pl.np.power(lambdax,2) + lambdax - 2)

def test4(lambdax):
  return pl.np.power(lambdax,4) - 20*pl.np.power(lambdax,3) + 0.1*lambdax
'''

def testFunction(order,lambdax):
  if args.choice == 1:
    return test1(order,lambdax)
  elif args.choice == 2:
    return -test2(order,lambdax)	
  elif args.choice == 3:
    return test3(order,lambdax)
  elif args.choice == 4:
    return test4(order,lambdax)
##########################################################################################  


###RUN SELECTED ALGORITHM
##########################################################################################  
if args.algorithm == 1:
  sol = lineSearchNewtonRapshon(a,b)
  algorithm_string = "1-D line search based on the Newton-Raphson Method"
elif args.algorithm == 2:
  sol = modifiedNewtonRapshon(a,b,L,r)
  algorithm_string = "Modification of the Newton-Raphson we discussed in class"
elif args.algorithm == 3:
  sol = goldenSection(a,b,L,r)
  algorithm_string = "Golden Section Method"
  order = 0

print "The following algorithm was used : " + algorithm_string
print "The following test function was tested : " + example_string
print "Total number of iterations:	%d" %sol[2]
print "lambda solution = %.3f" %sol[0]
print "f solution = %.3f" %sol[1]

fig1 = pl.figure(1)
pl.title(algorithm_string + " evaluated by " + example_string)
pl.xlabel("lambda")
pl.ylabel("f(lambda)")
#pl.figtext(0.40,0.882, "Total iterations: %d\nLambda solution: %d\nf(lambda solution): %d\n" %(sol[2],sol[0],sol[1]),fontsize=10, bbox=dict(facecolor = 'white', alpha=1, linewidth = 1, edgecolor = 'black'), horizontalalignment = 'left',verticalalignment='top')

#is hierdie a en b weer global waardes?
print "a = %.2f" %a
print "b = %.2f" %b

j = a
k = 0
while k < 50:
  if args.choice == 2:
    graph.append(-testFunction(0,j))
  else:
    graph.append(testFunction(0,j))
    
  #graph1.append(testFunction(1,j))
  #graph2.append(testFunction(2,j))  
    
  j = j + L/50.0
  k = k + 1
  

pl.plot(pl.np.arange(a,b,L/50.0),graph)
#pl.plot(pl.np.arange(a,b,L/50.0),graph1)
#pl.plot(pl.np.arange(a,b,L/50.0),graph2)
if args.choice == 2:
  pl.figtext(0.50,0.35, "Total iterations: %d\nLambda solution: %.3f\nf(lambda solution): %.3f\na: %.2f\nb: %.2f\nepsilon: %.6f\nNewton-Raphson initial lambda value: %.3f" %(sol[2],sol[0],sol[1],a,b,epsilon1,x_initial),fontsize=10, bbox=dict(facecolor = 'white', alpha=1, linewidth = 1, edgecolor = 'black'), horizontalalignment = 'center',verticalalignment='top')
  pl.plot(sol[0], -sol[1], 'g.', markersize=20.0)			#adds a dot where the solution was found
else:
  pl.figtext(0.50,0.882, "Total iterations: %d\nLambda solution: %.3f\nf(lambda solution): %.3f\na: %.2f\nb: %.2f\nepsilon: %.6f\nNewton-Raphson initial lambda value: %.3f" %(sol[2],sol[0],sol[1],a,b,epsilon1,x_initial),fontsize=10, bbox=dict(facecolor = 'white', alpha=1, linewidth = 1, edgecolor = 'black'), horizontalalignment = 'center',verticalalignment='top')
  pl.plot(sol[0], sol[1], 'g.', markersize=20.0)			#adds a dot where the solution was found


F = pl.gcf()
#F.savefig("./graphs/algorithm2/testFunction4.pdf",bbox_inches='tight',format='PDF',pad_inches=0.1)
  
pl.show()