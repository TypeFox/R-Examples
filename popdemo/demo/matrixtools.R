#=================================================
#Useful tools for working with matrices
#References:
#Caswell 2001 Matrix population models
#(primitivity)
#Stott et al. 2010 Methods Ecol. Evol., 1, 242-252
#(reducibility and ergodicity)
#=================================================
#
#
#
#-------------------------------------------------
#USING MATLAB STYLE MATRIX NOTATION
#-------------------------------------------------
#
#Matlab style notation for matrices
#popdemo contains a function enabling Matlab
#style notation of matrices. As well as being
#easier to use, it crucially enables import of 
#many matrices simultaneously, using comma-
#seperated csv files imported as dataframes.
#
#This matrix is for the desert tortoise
#(Gopherus Agassizii) with medium fecundity
#(Doak et al. 1994 Ecol. Appl., 4, 446-480).
Tort<-Matlab2R("[0 0 0 0 0 1.3 1.98 2.57;0.716 0.567 0 0 0 0 0 0;0 0.149 0.567 0 0 0 0 0;0 0 0.149 0.604 0 0 0 0;0 0 0 0.235 0.56 0 0 0;0 0 0 0 0.225 0.678 0 0;0 0 0 0 0 0.249 0.851 0;0 0 0 0 0 0 0.016 0.86]")
Tort
#
readline("Hit <Enter> for next topic, or <Esc> to cancel")
#
#
#
#-------------------------------------------------
#PRIMITIVITY, REDUCIBILITY & ERGODICITY
#-------------------------------------------------
#
#popdemo provides a few tools for diagnosing the
#ergodic properties of the matrix.
#
#Is the matrix primitive?
is.matrix_primitive(Tort)
#
#Is the matrix reducible?
is.matrix_irreducible(Tort)
#
#Is the matrix ergodic?
is.matrix_ergodic(Tort)
#
readline("Hit <Enter> for next topic, or <Esc> to cancel")
#
#
#
#-------------------------------------------------
#DEALING WITH REDUCIBLE MATRICES
#-------------------------------------------------
#
#Create a reducible matrix
TortR<-Tort
TortR[7,6]<-0
is.matrix_irreducible(TortR)
#
#Block-permute the reducible matrix
blockmatrix(TortR)
TortRblock<-blockmatrix(TortR)$blockmatrix
eigen(TortR)$values
eigen(TortRblock[1:6,1:6])$values
eigen(TortRblock[7:8,7:8])$values
#
#The first diagonal block has the dominant
#eigenvalue and the second diagonal block has
#the first subdominant eigenvalue. Therefore
#the matrix should still be ergodic.
is.matrix_ergodic(TortR)
#
#
#
#END
