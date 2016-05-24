### hplugin.R  (2006-01)
###
###                 Estimation of bandwith for step B
###             (plugin method)
###
### Copyright 2006-01 Sophie Lambert-Lacroix and Julie Peyre
###
### the method was described by Fan, J. and Gijbels, I. in
### "Local Polynomial Modelling and its Applications", Chapman and Hall,
### 1996, London.
### Implementation is largely inspired of the one proposed by Carroll et al in
### "Generalized Partially Linear Single-Index Models", Journal of the American
### Statistical Association, 92 (438), pages 477+, 1997.
###
### This file is part of the `plsgenomics' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

hplugin <- function(newXtrain,Ytrain, tau, intgsize,NbrIterMax=15)
{
##  INPUT variables
####################
##
##  newXtrain   : vector of length ntrain
##      train data matrix projected on the estimated projection direction
##  Ytrain   : vector of length ntrain
##      response variable {0,1}-valued vector
##  Xtest   : vector of length ntest
##      test data matrix projected on the estimated projection direction
##  tau, intgsize :
##
##  NbrIterMax : positive integer
##      max number of iteration in the WIRRLS part
##
##
##  OUTPUT variables
#####################
##  Structure with fields
##      h : value of hB estimated by plugin for step B
##      Cvg : logical
##      1 if convergence in WIRRLS 0 otherwise
#################################
## Choice of hB by plugin
#################################
ntrain <- length(newXtrain)
Umat <- cbind(rep(1,length=ntrain),newXtrain,newXtrain**2)

pfit <- (wirrls(Y=Ytrain,Z=Umat,Lambda=0,NbrIterMax=NbrIterMax,WKernel=diag(rep(1,ntrain))))
cvg <- pfit$Cvg
rm(Umat)
pcoefs <- pfit$Coefficients

q0 = pcoefs[1]
q1 = pcoefs[2]
q2 = pcoefs[3]

# Work out the summation involving the weight function
# and the denominator
wlow <-  min(newXtrain) + tau*(max(newXtrain) - min(newXtrain))
wupp <- max(newXtrain) - tau*(max(newXtrain) - min(newXtrain))

wsum = length(newXtrain[newXtrain<=wupp]) -length(newXtrain[newXtrain<wlow])
denom = 4*q2^2*wsum/ntrain

if (denom==0)
{
   h <- Inf  
   warning("Message from hplugin.R: Data well separated, parametric fit in second step")
}
else
{
 #  Now work out the integral required for the numerator 
   intgrid <- seq(wlow,wupp,length=intgsize)
   provmEta <- -q0 - q1*intgrid - q2*intgrid**2
   heights <- (exp(provmEta/2)+exp(-provmEta/2))**2

   gap <- (wupp - wlow)/(intgsize - 1)
   numint <- gap*sum(heights)

   # for Gaussian Kernels
   Ck <-(4*pi)**(-1/10);

   h = Ck*(numint/(denom*ntrain))**(1/5);
}
#############
# Conclusion
#############
return(list(h=h,Cvg=cvg))
}
