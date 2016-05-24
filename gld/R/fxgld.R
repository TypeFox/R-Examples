dgl <- function(x,lambda1=0,lambda2=NULL,lambda3=NULL,lambda4=NULL,param="fkml",
  lambda5=NULL,inverse.eps=.Machine$double.eps,max.iterations=500)
{
# Tidy the parameters so gl.check.lambda will work
lambdas <- .gl.parameter.tidy(lambda1,lambda2,lambda3,lambda4,param,lambda5)
# Check the parameters
if(!gl.check.lambda(lambdas,param=param,vect=TRUE)) {
        stop(paste("The parameter values", lambdas,
"\ndo not produce a proper distribution with the",param,
"parameterisation - see \ndocumentation for gl.check.lambda"))
        }
# calculate u=F(x) numerically, then use dqgl
# Unless x is outside the range, then density should be zero
extreme<-qgl(c(0,1),lambda1=lambdas,param=param)
# It may be better to change this to simply  
# (x <= extreme[2])*(x >= extreme[1])
outside.range <- !as.logical((x<=extreme[2])*(x>=extreme[1]))
# possibly calculate the end points here, rather than in the C code?
u <- pgl(x,lambdas,param=param,inverse.eps=inverse.eps,max.iterations=max.iterations)
dens <- dqgl(u,lambda1=lambdas,param=param)
dens[outside.range] <- 0
dens
}

pgl <- function(q,lambda1=0,lambda2=NULL,lambda3=NULL,lambda4=NULL,param="fkml",
    lambda5=NULL,inverse.eps=.Machine$double.eps,max.iterations=500)
{
# Thanks to Steve Su, see GLDEX package for contact details, for improvements to this code
# Tidy the parameters so gl.check.lambda will work
lambdas <- .gl.parameter.tidy(lambda1,lambda2,lambda3,lambda4,param,lambda5)
# Check the parameters
if(!gl.check.lambda(lambdas,param=param,vect=TRUE)) {
	stop(paste("The parameter values", lambda1, lambda2, lambda3,lambda4,
    	"\ndo not produce a proper distribution with the",param,
    	"parameterisation - see \ndocumentation for gl.check.lambda"))
    	} 
jr <- q; jr[sort.list(q)] <- seq(along=q) 
order.x<-order(q) 
xx<-sort(q) 
# xx has the sorted data, and jr & order.x the information to get back to the
# original order.
extreme<-qgl(c(inverse.eps,1-inverse.eps),lambda1=lambdas,param=param)
max.e<-extreme[2]
min.e<-extreme[1]
ind.min<-xx<=min.e
ind.max<-xx>=max.e 
# This simpler comparison works here because we are using inverse.eps as our
# tolerance
q<-xx[as.logical((xx<max.e)*(xx>min.e))] 
# We only want to calculate the probabilities for q values inside the support
length.of.vector <- length(q) 
# Need a blank u to send to C
u <- 0*q 
result <- switch(param, 
	freimer=, # allows for alternate expressions 
	frm=, # allows for alternate expressions 
	FMKL=, # Notes on .C call - the "numerics", lambdas and inverse.eps don't need the as.??? call as they are implicitly double
	FKML=, 
	fkml=, 
	fmkl=.C("gl_fmkl_distfunc",lambdas[1],lambdas[2],lambdas[3],lambdas[4], 
		as.double(0),as.double(1),inverse.eps,
		as.integer(max.iterations),as.double(q),as.double(u),
		as.integer(length.of.vector),PACKAGE="gld"),
    	ramberg=, # Ramberg & Schmeiser 
    	ram=, 
    	RS=, 
    	rs=.C("gl_rs_distfunc",lambdas[1],lambdas[2],lambdas[3],lambdas[4], 
    		as.double(0),as.double(1),inverse.eps,
		as.integer(max.iterations),as.double(q),as.double(u),
		as.integer(length.of.vector),
		PACKAGE="gld"), 
	fm5=.C("gl_fm5_distfunc",lambdas[1],lambdas[2],lambdas[3],
		lambdas[4],lambdas[5], 
    		as.double(0),as.double(1),inverse.eps,
		as.integer(max.iterations),as.double(q),as.double(u),
		as.integer(length.of.vector),PACKAGE="gld"),
    VSK=,
    gpd=,
    GPD=,
    vsk=.C("gl_vsk_distfunc",lambdas[1],lambdas[2],lambdas[3],lambdas[4],
        as.double(0),as.double(1),inverse.eps,
        as.integer(max.iterations),as.double(q),as.double(u),
        as.integer(length.of.vector),PACKAGE="gld"),
    stop("Error: Parameterisation must be one of fmkl, rs, fm5 or vsk")
    	) # closes "switch"
if (!(is.numeric(result[[1]]))){ 
	stop("Values for quantiles outside range. This shouldn't happen") 
} 
if (param=="fm5") {
	u <- result[[11]]
	}
else 	{
	u <- result[[10]]
	}
xx[as.logical((xx<max.e)*(xx>min.e))]<-u 
xx[ind.min]<-0 
xx[ind.max]<-1 
# Revert to the original order of the dataset: 
xx[jr] 
} 
