import numpy as np
import parser
from math import sin, cos, tan, asin, acos, atan
from math import exp, pow, log, log10, sqrt
from math import pi, e
from decimal import *
import matplotlib.pyplot as plot 
import time
import csv


def write_to_csv(name,  data_set):
	with open(name, 'w') as file:
		w = csv.writer(file, delimiter=',')
		w.writerows(data)


### This function class allows for some you to easily input
### Mathematical equations for processing

class Function:
	def __init__(self, f = "0"):
		self.f = f
		self.f_code = parser.expr(self.f).compile()

	def redefine(self, f = "0"):
		self.f = f
		self.f_code = parser.expr(self.f).compile()

	def eval(self,x,t=0,y = 0, z = 0):
		return eval(self.f_code)

def graph(x_data, y_data, x_ax = "x", y_ax ="y",title="", name = "", show=True, save=True, precision=0.0001):
	for i in range(len(y_data)):
		plot.plot(x_data[i],y_data[i])
	plot.title(title)
	plot.xlabel(x_ax)
	plot.ylabel(y_ax)
	if show: plot.show()
	if save: plot.savefig("diffeqs/"+name+".jpg")
	plot.gcf().clear()

def exact_graph(x,y,x_ax = "X", y_ax = "Y", title="", name="", show=True, save=True, precision=0.0001):
	plot.plot(x,y)
	plot.title(title)
	plot.xlabel(x_ax)
	plot.ylabel(y_ax)
	if show: plot.show()
	if save: plot.savefig("diffeqs/exact_"+name+".jpg")
	plot.gcf().clear()

def Matrix_A(n):
	l1 = [1 for i in range(n-1)]
	l2 = [4 for i in range(n)]
	l3 = (l1)

	A = np.diag(l2)
	B = np.diag(l1,k=1)
	C = np.diag(l3,k=-1)

	M = A+B+C
	return M

def Matrix_B(n):
	l1 = [-1 for i in range(n-1)]
	l2 = [2 for i in range(n)]

	A = np.diag(l2)
	B = np.diag(l1,k=1)
	C = np.diag(l1,k=-1)

	return A+B+C

def Matrix_C(n):
	l1 = [1 for i in range(n-1)]
	l2 = [-1 * i for i in l1]

	A = np.diag(l1,k=1)
	B = np.diag(l2, k=-1)

	return A+B

def Matrix_D(n):
	l1 = [1 for i in range(n)]
	l2 = [-1 for i in range(n-1)]

	A = np.diag(l1)
	B = np.diag(l2,k=-1)
	return A + B

def eq1(n,eps,f):
	K = Matrix_B(n)
	A = np.eye(n)
	f = Function(f)

	M = (K*((eps**2)*(n**2)))+(A)

	F = [f.eval(i*(1/float(n))) for i in range(n)]
	return np.linalg.solve(M,F)

def eq2(n,eps,f):
	K = Matrix_B(n)
	D = Matrix_C(n)
	f = Function(f)

	M = (K*(eps*(n**2)))+(D*(n/2.0))

	F = [f.eval(i*(1/float(n))) for i in range(n)]
	return np.linalg.solve(M,F)
def eq2_fixed(n,eps,f):
	K = Matrix_B(n)
	D = Matrix_D(n)
	f = Function(f)

	M = (K*(eps*(n**2)))+(D*(n))

	F = [f.eval(i*(1/float(n))) for i in range(n)]
	return np.linalg.solve(M,F)

def main_loop():
	eqs = ['1','x','x**2']
	actual2 = [
		'x + (e**(-1/t))/(1-e**(-1/t)) - (e**((x-1)/t))/(1-e**(-1/t))',
		'((x**2)/2.0) + (t*x) + ((0.5+t)/(1-e**(-1/t)))*(e**(-1/t)-e**((x-1)/t))',
		'(x**3/3.0)+(t*x**2)+(2*t*x) + (((1/3.0)+t+(2*t**2))/(1-e**(-1/t)))*(e**(-1/t)-e**((x-1)/t))'
	]
	actual1 = [
		'1 + ((e**(-1/t)-1)/(1-e**(-2/t)))*(e**((x-1)/t) + e**(-x/t))',
		'x + (e**((x-1)/t)/(e**(-2/t)-1)) + (e**((-x-1)/t)/(1-e**(-2/t)))',
		'x**2+(2*t**2) + ((1+2*t**2-2*t**2*e**(-1/t))/(e**(-2/t)-1))*e**((x-1)/t)+((e**(-1/t)+2*t**2*e**(-1/t)-2*t**2)/(1-e**(-2/t))) * e**(-x/t)'
	]	
	errors = []
	h = [log(i) for i in range(3,12)]
	for k in range(4):
		sub_error = []
		y = []
		x = []
		FX =[]
		eps = 10**(-1*k)
		for j in range(3,12):
			U = eq2_fixed(2**j,eps,'1')
			X = [i*(1/float(2**j)) for i in range(2**j)]
			y.append(U)
			x.append(X)

			fx=Function(actual2[0])
			FX = [fx.eval(i*(1/float(2**j)),eps) for i in range(2**j)]
			error = [(abs(FX[i]-U[i])) for i in range(2**j)]
			sub_error.append(max(error))

		fx = Function(actual2[0])
		FX = [fx.eval(i*(1/float(10000)),eps) for i in range(10000)]
		X = [i*(1/float(10000)) for i in range(10000)]
		errors.append(sub_error)
		exact_graph(h,sub_error,"Log(h)","Log(error)", "Error Graph", "fx=1_error_e="+str(eps))
		graph(x,y,"x","y","Approximate Differential Equation","fx=1_e="+str(eps))
		exact_graph(X,FX,"x","y","Actual Differential Equation", "fx=1_actuale="+str(eps))
main_loop()