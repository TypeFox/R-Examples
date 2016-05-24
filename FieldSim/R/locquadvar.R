### quadvar.R  (2006-16)
###
###    Estimation of the multifractional function of the multi-fractal Brownian field
###             by the localized quadratic variations method
###
### Copyright 2006-16 Alexandre Brouste and Sophie Lambert-Lacroix
###

locquadvar <- function(Z,t,h){

##    INPUT VARIABLES
#########################
##  Z   : matrix of size NxN 
##        matrix associated with the sample path of one  
##        multi-fractal Brownian field
##        (N=2^(nblevel)+1 and Z[i,j] gives the value of 
##        the process at point (2^(-i*nblevel),2^(-j*nblevel))) 
##  t   : point where H(t) must be estimated
##  h   : the neighborhood used to estimate H(t) is of the following form
##        (2^(-i*nblevel),2^(-j*nblevel)), i,j=1,...,N such that
##        |2^(-i*nblevel)-t_1|<=h and |2^(-j*nblevel)-t_2|<=h


##    OUTPUT VARIABLES
##########################
## H: real in ]0,1[
##      Estimation of the multifractional function of the multi-fractal Brownian field
##      at the point t


##  TEST ON INPUT VARIABLES
##############################

##On Z
if ((is.numeric(Z)==FALSE)||(is.matrix(Z)==FALSE)){
 stop("Message from locquadvar.R: Z is not of valid type")}

if (dim(Z)[1]!=dim(Z)[2]){
 stop("Message from locquadvar.R: Z is not of valid type")}

##On h 
if (is.vector(h)==FALSE)
   stop("Message from locquadvar.R: h is not of valid type")
if (length(h)!=1)
   stop("Message from locquadvar.R: only one value can be specified for h")
if ((is.numeric(h)==FALSE)||(h<=0)){
 stop("Message from locquadvar.R: h is not of valid type")}

##On t 
if (is.vector(t)==FALSE)
   stop("Message from locquadvar.R: t is not of valid type")
if (length(t)!=2)
   stop("Message from locquadvar.R: t must be a vector of length 2")
if ((is.numeric(t)==FALSE)||(t[1]<0)||(t[2]<0)||(t[1]>1)||(t[2]>1)){
 stop("Message from locquadvar.R: t is not of valid type (2)")}


#Find the points in the neighborhood of t
############################################
x1 <- max(0,t[1]-h)
x2 <- min(1,t[1]+h)
y1 <- max(0,t[2]-h)
y2 <- min(1,t[2]+h)
N <- dim(Z)[1]
M <- N-1
i1 <- floor(x1*M)+1
i2 <- floor(x2*M)+1
j1 <- floor(y1*M)+1
j2 <- floor(y2*M)+1
Vh <- Z[i1:i2,j1:j2]

#Compute the estimate
#########################
H <- quadvaraux(Vh)

return(H)}
