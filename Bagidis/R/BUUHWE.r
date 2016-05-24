BUUHWE= function(x){

#===============================================================================
# FUNCTION : BUUHWE
#-------------------------------------------------------------------------------
# BY: Catherine Timmermans, Institute of Statistics, Biostatistics and Actuarial Sciences, UCLouvain. 
# Last update: 8/15/2010 10:59:17 AM
# Contact: catherine.timmermans@uclouvain.be
#-------------------------------------------------------------------------------
# DESCRIPTION: 
# function of a time series x returning its unbalanced Haar wavelet expansion, as obtained with a bottom-up algortithm.
# BUUHWE is an acronym for Bottom-Up Unbalanced Haar Wavelets Expansion.
#-------------------------------------------------------------------------------
# USAGE: 
# BUUHWE(x)
#-------------------------------------------------------------------------------
# ARGUMENTS: 
# x: a numeric vector of length N
#-------------------------------------------------------------------------------
# DETAILS:
#
# Unbalanced Haar wavelet bases are orthonormal bases that are made up of one constant vector and a set of Haar-like (i.e. \textit{up-and-down}-shaped ) orthonormal wavelets whose discontinuity point between positive and negative parts is not necessarily located at the middle of its support. Using the notation of \citet{Piotr}, the general mathematical expression of one of those Haar-like wavelets is given by
# \begin{equation}\label{HaarVec}
# \boldsymbol{\phi}_{e,b,s}(t) = (\frac{1}{b-s+1}- \frac{1}{e-s+1})^{1/2} \, \Ind_{s \le t \le b} 
#                + (\frac{1}{e-b}- \frac{1}{e-s+1})^{1/2} \, \Ind_{b+1 \le t \le e},  \nonumber
# \end{equation}
# where $t =1 \ldots N$ is a discrete index along the abscissa axis, and where $s < b < e$ are values along this axis that determine the particular shape of one wavelet ($s$, $b$ and $e$ stand for \textit{start}, \textit{breakpoint} and \textit{end} respectively). Each wavelet $\boldsymbol{\phi}_{e,b,s}(t)$ is thus associated with a level change from one observation (or group of observations) to the consecutive one, and the projection of the series $\mathbf{x}(t)$ on the wavelet $\boldsymbol{\phi}_{e,b,s}(t)$ encodes the importance of the related  level change in the series. Different choices of $N-1$ sets of values $s$, $b$ and $e$ leading to orthonormal wavelets define a whole family of bases $\{ \boldsymbol{\phi}_k\}_{k=0 \ldots N-1}$.
# 
#\paragraph{The Basis Pursuit Algorithm and the Property of Hierarchy.} In 2007, \citet{Piotr} proposed an algorithm for building the unbalanced Haar wavelet basis $\{ \boldsymbol{\phi}_k\}_{k=0 \ldots N-1}$ that is best suited to a given series, according to the principle of hierarchy - namely, the vectors of this basis and their associated coefficients are ordered following the importance of the level change they encode  for describing the global shape of the series. He called it the \textit{bottom-up unbalanced Haar wavelet transform} (here-after BUUHWT). The resulting expansion is organized in a hierarchical way and avoids the dyadic restriction that is typical for classical wavelets. The family of unbalanced Haar wavelets is thus really adaptive to the shape of the series.
#
#-------------------------------------------------------------------------------
# VALUE:
# a list whose elements are
# detail: detail coefficients starting from rank k=0 up to k=N-1 
# basis:  unbalanced Haar basis vectors ordered by colums. First column is the constant vector of rank k=0. This is a matrix of dim N X N.
# split.abs: localization index. For consistency with basis matrix and detail vector dimensions, this is a vector of length N but first coefficient is NA. Index are then ordered by increasing rank k.
# series: the initial series x.
#-------------------------------------------------------------------------------
# WARNING/NOTE :
# This function can be compared to function uh.bu in package unbalhaar of Piotr Fryzlewicz. 
#-------------------------------------------------------------------------------
# REFERENCES: 
# Timmermans C., von Sachs R., 2010, BAGIDIS, a New Method for Statistical Analysis of Differences Between Curves with Sharp Discontinuities, discussion paper 1030 of Institute of Statistics, Biostatistics and Actuarial Sciences. 
# Fryzlewicz P., 2007,  Unbalanced Haar Technique for Non Parametric Function Estimation, Journal of the American Statistical Association, 102, 1318-1327.
#-------------------------------------------------------------------------------
# SEE ALSO:  
# BAGIDIS.dist, semimetric.BAGIDIS, BUUHWE.plot,...
#-------------------------------------------------------------------------------
# EXAMPLES:
#> x= c(1,7,3,0,-2,6,4,0,2)
#> BUUHWE(x)
#$detail
#[1]  7.000000  2.828427 -4.618802  4.000000
#[5] -3.265986  2.828427 -1.414214  1.414214
#[9]  1.414214
#
#$basis
#                  new.vect   new.vect
# [1,] 0.3333333  0.4714045  0.0000000
# [2,] 0.3333333  0.4714045  0.0000000
# [3,] 0.3333333  0.4714045  0.0000000
# [4,] 0.3333333 -0.2357023  0.5773503
# [5,] 0.3333333 -0.2357023  0.5773503
# [6,] 0.3333333 -0.2357023 -0.2886751
# [7,] 0.3333333 -0.2357023 -0.2886751
# [8,] 0.3333333 -0.2357023 -0.2886751
# [9,] 0.3333333 -0.2357023 -0.2886751
#      new.vect   new.vect   new.vect
# [1,]      0.0  0.8164966  0.0000000
# [2,]      0.0 -0.4082483  0.7071068
# [3,]      0.0 -0.4082483 -0.7071068
# [4,]      0.0  0.0000000  0.0000000
# [5,]      0.0  0.0000000  0.0000000
# [6,]      0.5  0.0000000  0.0000000
# [7,]      0.5  0.0000000  0.0000000
# [8,]     -0.5  0.0000000  0.0000000
# [9,]     -0.5  0.0000000  0.0000000
#        new.vect   new.vect   new.vect
# [1,]  0.0000000  0.0000000  0.0000000
# [2,]  0.0000000  0.0000000  0.0000000
# [3,]  0.0000000  0.0000000  0.0000000
# [4,]  0.0000000  0.0000000  0.7071068
# [5,]  0.0000000  0.0000000 -0.7071068
# [6,]  0.0000000  0.7071068  0.0000000
# [7,]  0.0000000 -0.7071068  0.0000000
# [8,]  0.7071068  0.0000000  0.0000000
# [9,] -0.7071068  0.0000000  0.0000000
#
#$split.abs
#[1] NA  3  5  7  1  2  8  6  4
#
#$series
#[1]  1  7  3  0 -2  6  4  0  2
#
#===============================================================================
#===============================================================================



x.init =x 

N =length(x)

axis.weights =rep(1,N)
group =rep(1,N)


split.rel = NULL
sin.phi.vect = NULL
detail.vect = NULL

new.basis = NULL
split.abs =NULL

for (step in (1:(N-1))){ ## START LOOP

#-------------------------------------------------------------------------------
# 1/ Look for next basis vector
#------------------------------------------------------------------------------- 

M = N-step +1

x.plus = x[-1]
x.minus=x[-M]                                         
axis.weights.plus =axis.weights[-1]
axis.weights.minus =axis.weights[-M]

norm = sqrt(axis.weights.plus^2 + axis.weights.minus^2)
cos.phi = axis.weights.plus/norm
sin.phi = axis.weights.minus/norm

detail = cos.phi*x.minus - sin.phi*x.plus
somme  = sin.phi*x.minus + cos.phi*x.plus

i = min(which(abs(detail) == min(abs(detail))))

#-------------------------------------------------------------------------------       
# 2/ Store information about the new basis.
#-------------------------------------------------------------------------------       

# > detail coefficient, relative split index, sin.phi.vect
split.rel = c(split.rel,i)
sin.phi.vect = c(sin.phi.vect, sin.phi[i])
detail.vect =  c(detail.vect,detail[i])

# > basis vector.

new.vect =rep(0,N)
middle.index = sum(group[1:i]) +1
start.index = middle.index - group[i]
end.index = middle.index + group[i+1]

na=group[i]
nb=group[i+1]
ntot =na+nb

new.vect[start.index:(middle.index-1)]=  sqrt(nb/na*(1/ntot)*(1/cos.phi[i]^2) ) * cos.phi[i]

new.vect[middle.index:(end.index-1)] = - sqrt(na/nb*(1/ntot)*(1/sin.phi[i]^2)  ) *sin.phi[i]

new.basis = cbind(new.vect, new.basis)

split.abs = c(split.abs,(middle.index-1))

#------------------------------------------------------------------------------- 
# 3/ Evolution
#------------------------------------------------------------------------------- 


# x evolution 
x=x[-i]
x[i]=somme[i]

# weights evolution
axis.weights = axis.weights[-i]
axis.weights[i] = norm[i]

# group evolution

group[i] = group[i]  + group[i+1]
group = group[-(i+1)]   

}## END LOOP.

#-------------------------------------------------------------------------------       
# 4/last basis vector and coefficient
#-------------------------------------------------------------------------------
       
new.basis =cbind(rep(sqrt(N)/N,N),new.basis)
detail.vect = c(detail.vect,x)
split.rel = c(split.rel,NA)
sin.phi.vect = c(sin.phi.vect, NA)

split.abs =c(split.abs,NA)

#------------------------------------------------------------------------------- 
# 5/ Output
#------------------------------------------------------------------------------- 

#Piotr =rbind(split.rel ,sin.phi.vect, detail.vect)    #could be calculated for comparisons with Piotr's function uh.bu() outputs.

return(list("detail"=rev(detail.vect),"basis" = new.basis, "split.abs" = rev(split.abs), "series" = x.init))

}

