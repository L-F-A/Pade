import numpy as np
#################################################################
#   This implement pade fitting for analytical continuation	#
#								#
# wn :	  Positive Matsubara frequency wn. Numpy 1d array	#
# fwn :   Function f(iwn) to be continued. Numpy 1d array	#
# N :	  Use N wn for the fitting				#
#################################################################

class Pade():

	def __init__(self,wn,fwn,N):
		self.wn=wn[:N]
		self.fwn=fwn[:N]
		self.N=N

	#Findthe coefficients of the continued fraction expansion from 
	#effective matrix approach
	def coeffs_mod(self):
	#Contructing the g's
    	#This matrix is artificial is the sense that the lower triangular part is
    	#NAN always but that does not matter as we are only interested, at the end,
    	#at the diagonal
		coef=np.zeros(self.N)
		g=np.zeros(self.N,self.N)
		g[0,:]=self.fwn
		g[1,:]=(self.fwn[0]-self.fwn)/(( 1j*self.wn-1j*self.wn[0] )*self.fwn)
		for k in range(2,self.N):
			g[2,:]=( g[k-1,k-1]-g[k-1,:] )/((1j*self.wn-1j*self.wn[k-1])*g[k-1,:])
		self.coef=np.diag(g)
	#Find the coefficients of the continued fraction expansion
	def coeffs(self):
		coef=np.zeros(self.N)
		for j in range(self.N):
			coef(j)=__pade_recursion(self,j,j,coef)
		self.coef=coef

	def __pade_recursion(self,idFunc,idx,coef):		
		if idFunc > 0:
			f0=coef(idFunc-1)
			f1=self.__pade_recursion(self,idFunc-1,idx,coef)
			v=(f0-f1)/( (1j*self.wn[idx]-1j*self.wn[idFunc-1])*f1  )
		else:
			v=self.fwn[idx]
		return v

	def query(self,w,eta):
	# w :	Real frequency at which we want to know the function
	# eta :	Small imaginary part. Note that it can often safely be put equal 
	#	to zero since it appears only in expression of the form w+i*eta-iwn
	#	and thus unless wn is very small or eta very big, it will not change anything
		z=w+1j*eta
		A=1.
		B=1. + self.coef[1]*(z-1j*self.wn[0])
		P=A/B
		for c in range(2,self.N):
			A=1.+self.coef[c]*(z-1j*self.wn[c-1])/A
			B=1.+self.coef[c]*(z-1j*self.wn[c-1])/B
			P=P*A/B
		return P*a[0] 
