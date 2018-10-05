import numpy as np
'''
#################################################################
#   This implement pade fitting for analytical continuation	#
#								#
#   written by Louis-François Arsenault				#
#								#
# wn  :	  Positive Matsubara frequency wn. Numpy 1d array	#
# fwn :   Function f(iwn) to be continued. Numpy 1d array	#
# N   :	  Use N wn for the fitting				#
#								#
# wn and fwn must be 1d array, very important as I did not put a#
# testing test throwing warning etc for the moment		#
#								#
# example:							#
#	#Say we have a function fwn in Matsubara frequency 	#
#	#wn and we want to do Pade using the first 100 wn	#
#	N=100							#
#	pade=Pade(wn,fwn,N)					#
#	pade.coeffs()						#
#								#
#	#If for example we want to know the function on real	#
#	axis for 500 freqencies from -1 to 1 with eta=0		#
#	w=np.linspace(-1,1,500)					#
#	fw=np.zeros(500)					#
#	for r in range(500):					#
#		fw[r]=pade.query(w[r],0)			#
#################################################################
'''
class Pade():

	def __init__(self,wn,fwn,N):
		self.wn=wn[:N]
		self.fwn=fwn[:N]
		self.N=N

	#Find the coefficients of the continued fraction expansion
	def coeffs(self):
		coef=np.zeros(self.N)*(1.+1j)
		for j in range(self.N):
			coef[j]=self.__pade_recursion(j,j,coef)
		self.coef=coef

	def __pade_recursion(self,idFunc,idx,coef):		
		if idFunc > 0:
			f0=coef[idFunc-1]
			f1=self.__pade_recursion(idFunc-1,idx,coef)
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
		return P*self.coef[0] 
