# R package for air multivariate time series data imputation based o EM algorithm
# IMS - UERJ - Brasil
# written by Washington Junger <wjunger@ims.uerj.br>
# restarted on 01/09/2004
# original work: EM Library for R
# Last change: 10/07/2002; 30/05/2003 (R version)
#
#---------------------------------
# last changes list: (date first)
# 01/09/2004 - first beta release - v. 0.1.0
# 21/10/2004 - added smooth spline method for filtering - v. 0.1.1
# 04/05/2005 - new interface, new features, C code for fastening, code rewriting etc - v. 0.2.0
# 01/12/2005 - variance window, now it returns some measure of uncertainty, added plot function - v. 0.2.1
# 30/06/2006 - fixed bug where data were not time series
# 25/07/2006 - documentation fixed - v. 0.2.2
# 25/10/2006 - documentation and some minor bugs fixed - v. 0.2.3
# 04/01/2007 - some improvements and bug fixes - v. 0.2.4
# 14/06/2007 - some improvements and bug fixes - v. 0.2.5
# 20/06/2007 - some improvements and bug fixes - v. 0.2.6
# 08/11/2012 - fix namespace, startup functions, and warnings in em.spline - v. 0.3.3
#

#---------------------------------
# To do list:
# 
#
#



mstats <- function(dataset)
# Carries some missing data statistics out
# dataset - data matrix
{
n <- dim(dataset)[1]
p <- dim(dataset)[2]

pollnames <- names(dataset)
colnames <- c("missing", "missing(%)")
rownames <- c(1:n)

# getting status of columns
colmissing <- matrix(0,p,2)
dimnames(colmissing) <- list(pollnames,colnames)

for (j in 1:p)
    {
	colmissing[j,1] <- sum(is.na(dataset[,j]))
	colmissing[j,2] <- 100*sum(is.na(dataset[,j]))/n
	}

# getting status of rows
rowmissing <- matrix(NA,n,2)
dimnames(rowmissing) <- list(rownames,colnames)

rowcounts <- em.extractcoord(dataset)[[2]][,2]
relrowcounts <- 100*rowcounts/p

rowmissing[,1:2] <- cbind(rowcounts,relrowcounts)
# a little bit of information
maxmissing <- max(rowmissing[,1],na.rm=TRUE)
counts <- matrix(0,p+1,2)
dimnames(counts) <- list(as.character(c(0:p)),colnames)

for (k in 0:p)
	{
	counts[(k+1),1] <- sum(ifelse(rowcounts[1:n]==k,1,0))
	counts[(k+1),2] <- 100*counts[(k+1),1]/n
	}

retval <- list(rows=rowmissing,columns=colmissing,pattern=counts)
return(retval)
}


edaprep <- function(dataset)
# Prepares a dataset for exploratory data analysis
# ds - data set
{
mean <- em.mean(dataset)	# estimates Xi mean
coord <- em.extractcoord(dataset)	# extracts NA coordinates
retval <- em.replacewmean(dataset,mean,coord,TRUE)	# replaces NA with mean
return(retval)	 # returning value
}


mkjnw <- function()
# Creates Johnson & Wichern example matrix
{
retval <- as.data.frame(matrix(c(NA,0,3,7,2,6,5,1,2,NA,NA,5),ncol=3,byrow=TRUE))	# returning value
colnames(retval) <- c("X1","X2","X3")
return(retval)
}
	
	
elapsedtime <- function(st,et)
# Computes the elapsed time between t1 and t2
# t1 - start time
# t2- end time
{
time <- et[3]-st[3]		# gets time from the third position of the vectors
h <- trunc(time/3600)
if (h<10)
	hs <- paste("0",h,sep="",collapse=" ")
else
	hs <- h
time <- time-h*3600
min <- trunc(time/60)
if (min<10)
	mins <- paste("0",min,sep="",collapse=" ")
else
	mins <- min
time <- time-min*60
sec <- trunc(time)
if (sec<10)
	secs <- paste("0",sec,sep="",collapse=" ")
else
	secs <- sec

retval <- paste(hs,":",mins,":",secs,sep="",collapse=" ")
return(retval)
}


em.mean <- function(m)
# Calculates mean vector regardless NA
# m - data matrix 
{
retval <- t(apply(m,2,mean,na.rm=TRUE))	# returns the mean function applied to all series
return(retval)
}


em.det <- function(x)
# Calculates the determinant of a matrix
# x - matrix to calculate the determinant
{
retval <- prod(eigen(x)$values)		# returns the product of the eigenvalues
return(retval)
}	


em.trace <- function(x)
# Calculates the trace of a matrix
# x - matrix to calculate the trace
{
retval <- sum(diag(x))	# returns the summation of main diagonal elements
return(retval)
}	


em.existna <- function(v)
# Checks existence and counts missing values in a given vector
# v - vector to be checked
{
na.count <- sum(is.na(v))
		
retval <- list(naexist=na.count>0,nacount=na.count)		# returning values
return(retval)
}


em.extractcoord <- function(m)
# Extracts NA coordinates
# m - data matrix	
{
r <- dim(m)[1]		# internal function variables initialization
c <- dim(m)[2]      
n <- 0

for (i in 1:r)		# checks if the ith element is NA and counts
	{
	for (j in 1:c)
		if (is.na(m[i,j]))
			n <- n+1
	}
mm.1 <- matrix(NA,nrow=n,ncol=2,byrow=TRUE)	 # internal function variables initialization
mm.2 <- matrix(0,nrow=r,ncol=2,byrow=TRUE)
n <- 0  # internal function variables initialization
		
for (i in 1:r)	# extracts NA coordinates
	{
	k <- 0
	for (j in 1:c)
		if (is.na(m[i,j]))
			{
			n <- n+1
			k <- k+1
			mm.1[n,1] <- i
			mm.1[n,2] <- j
			}
	mm.2[i,1] <- i
	mm.2[i,2] <- k
	}
retval <- list(nacoord=mm.1,nainrow=mm.2)  # returning values
return(retval)			
}
	
	
em.countmvec <- function(x)
# Counts missing vector in the given data matrix
# x - data matrix
{
mr <- apply(x,1,em.existna)	 # count missing values per row
	
retval <- sum(mr$nacount>=dim(x)[2])	 # if all the elements in the vector are missing then counts it as missing vector
return(retval)
}


em.countmnear <- function(x)
# counts consecutive missing values in columns
{
rows <- nrow(x)
cols <- ncol(x)
mc <- matrix(0,nrow=rows,ncol=cols)
for(j in 1:cols)
	{
	for(i in 1:rows)
		{                                                                                
		count <- 0                                                                                    
		k <- i                                                                                        
		while(is.na(x[k,j])) 
			{                                                         
			count <- count+1                                                                            
			k <- k-1
			if(k==0)
				break
			}                                                                               
		mc[i,j] <- count
		}                                                                       
	} 
return(mc)
}


em.replacewmean <- function(m,mu,c,diag=FALSE)
# Replaces NA with mean
# m - data matrix
# mu - mean vector
# c - coordinates list
# diag - diagnostics (used for arima diagnostics)
{
n <- dim(c$nacoord)[1]	# internal function variable initialization
for (i in 1:n)
	# if ((c$nainrow[c$nacoord[i,1],2] != length(mu)) || (diag == TRUE))	 # if not all the elements are 	missing or diagnostics call (deprecated)
		# then fill in the missing elements with series mean
		m[c$nacoord[i,1],c$nacoord[i,2]] <- mu[c$nacoord[i,2]]
	
return(m)
}


em.dispersion <- function(m,v=TRUE)
# Estimates variance or correlation
# m - data matrix
# v - TRUE:variance FALSE:correlation
{
if (v == TRUE)
	retval <- cov(m,use="complete.obs")		# returns variance matrix
else
	retval <- cor(m,use="complete.obs")		# 	returns correlation matrix
return(retval)
}


em.arrangevec <- function(x,sv)
# Arranges a vector based in the given order
# x - vector to be arranged
# sv - split vector
{
lna <- dim(sv$part1)[1]		# internal function variables initialization
lknown <- dim(sv$part2)[1]
newx <- numeric(length(x))

for (i in 1:lna)	# stores missing elements at the beginning of the vector
	newx[i] <- x[sv$part1[i,1]]	
for (i in 1:lknown)		# stores known elements at the end of the vector
	newx[i+lna] <- x[sv$part2[i,1]]

return(newx)
}
	
	
em.rearrangevec <- function(x,sv)
# Rearranges a vector to its original order
# x - vector to be rearranged
# sv - split vector
{
lna <- dim(sv$part1)[1]		# internal function variables initialization
lknown <- dim(sv$part2)[1]
originalx <- numeric(length(x))
	
for (i in 1:lna)	# restores estimated elements to the original missing location
	originalx[sv$part1[i,1]] <- x[i]	
for (i in 1:lknown)		# restores known elements to their original location
	originalx[sv$part2[(i),1]] <- x[i+lna]

return(originalx)
}


em.arrangemat <- function(x,sv)
# Arranges a matrix based on a given order
# x - matrix to be arranged
# sv - split vector
{
lna <- dim(sv$part1)[1]		# internal function variables initialization
lknown <- dim(sv$part2)[1]
newx.t <- matrix(NA,nrow=dim(x)[1],ncol=dim(x)[2],byrow=TRUE)
newx <- matrix(NA,nrow=dim(x)[1],ncol=dim(x)[2],byrow=TRUE)
	
for (i in 1:lna)	# moves rows related to the missing elements to the beginning of the matrix
	newx.t[i,] <- x[sv$part1[i,1],]	
for (i in 1:lknown)		# moves rows related to the known elements to the end of the matrix
	newx.t[(i+lna),] <- x[sv$part2[i,1],]
for (j in 1:lna)	# moves columns related to the missing elements to the beginning of the matrix
	newx[,j] <- newx.t[,sv$part1[j,1]]	
for (j in 1:lknown)		# moves columns related to the missing elements to the end of the matrix
	newx[,(j+lna)] <- newx.t[,sv$part2[j,1]]
		
return(newx)
}
	
	
em.rearrangemat <- function(x,sv)
# Rearranges a matrix to its original order
# x - matrix to be rearranged
# sv - split vector
{
lna <- dim(sv$part1)[1]		# internal function variables initialization
lknown <- dim(sv$part2)[1]
newx.t <- matrix(NA,nrow=dim(x)[1],ncol=dim(x)[2],byrow=TRUE)
newx <- matrix(NA,nrow=dim(x)[1],ncol=dim(x)[2],byrow=TRUE)

for (i in 1:lna)	# restores rows related to estimated values to the original missing location
	newx.t[sv$part1[i,1],] <- x[i,]	
for (i in 1:lknown)		# restores rows related to known values to the original location
	newx.t[sv$part2[i,1],] <- x[(i+lna),]
for (j in 1:lna)	# restores columns related to estimated values to the original missing location
	newx[,sv$part1[j,1]] <- newx.t[,j]	
for (j in 1:lknown)		# restores columns related to known values to the original location
	newx[,sv$part2[j,1]] <- newx.t[,(j+lna)]

return(newx)
}


em.partmu <- function(mu,x)
# Divides mean vector in 2 partitions
# mu - mean vector
# x - partitioning position
{
mu.1 <- mu[1:x]	 # stores first partition of the vector
mu.2 <- mu[(x+1):length(mu)]	 # stores second partition of the vector

retval <- list(mu1=mu.1,mu2=mu.2)
return(retval)
}
	

em.partsigma <- function(sigma,x)
# Divides covariance matrix in 4 partitions
# sigma - covariance matrix
# x - partitioning position
{
c <- dim(sigma)[2]		# sub-matrices related to the partitions initialization

sigma.11 <- sigma[1:x,1:x]	 # stores values to the partition S11
sigma.22 <- sigma[(x+1):c,(x+1):c]	 # stores values to the partition S22
sigma.12 <- sigma[1:x,(x+1):c]	 # stores values to the partition S12
sigma.21 <- sigma[(x+1):c,1:x]	 # stores values to the partition S21

retval <- list(sigma11=sigma.11,sigma12=sigma.12,sigma21=sigma.21,sigma22=sigma.22)
return(retval)
}
	

em.splitvec <- function(v)
# Splits a vector in 2 parts: NA part and non-NA one
# v - vector to be split
{
l <- length(v)		# internal function variables initialization
n <- 0

for (j in 1:l)		# counts missing elements
	if (is.na(v[j]))
		n <- n+1
	
v.1 <- matrix(NA,nrow=n,ncol=2,byrow=TRUE)		# internal function variables initialization
v.2 <- matrix(NA,nrow=l-n,ncol=2,byrow=TRUE)
jv1 <- 1	# internal function variables initialization
jv2 <- 1

for (j in 1:l)
	if (is.na(v[j]))	# stores missing 	elements in the first part
		{
		v.1[jv1,1] <- j
		v.1[jv1,2] <- v[j]
		jv1 <- jv1+1
		}
	else	# stores known 	elements in the second part
		{
		v.2[jv2,1] <- j
		v.2[jv2,2] <- v[j]
		jv2 <- jv2+1
		}

retval <- list(part1=v.1,part2=v.2)
return(retval)		
}
	

em.correctmat <- function(m,l)
# Substitutes original values in the matrix by those estimated by EM algorithm
# m - matrix to be corrected
# l - list with the estimated values
{
r.l.1 <- dim(l$cont2na)[1]		# internal function variables initialization
c.l.1 <- dim(l$cont2na)[2]
r.l.2 <- dim(l$cont2obs)[1]
c.l.2 <- dim(l$cont2obs)[2]
for (i in 1:r.l.1)		# overwrite with the missing X missing contribution to T2
	{
	for (j in 1:c.l.1)
	m[i,j] <- l$cont2na[i,j]
	}
for (i in 1:r.l.2)		# overwrite with the missing X known contribution to T2
	{
	for (j in 1:c.l.2)
		{
		m[i,(c.l.1+j)] <- l$cont2obs[i,j]
		m[(c.l.1+j),i] <- l$cont2obs[i,j]
		}
	}

return(m)
}
	

em.contribt1 <- function(mp,sp,sv)
# Estimates the NA vector (contribution to T1)
# mp - mean vector
# sp - partitioned covariance matrix
# sv - NA/non <- NA split vector
{
retval <- mp$mu1+sp$sigma12%*%solve(sp$sigma22)%*%(sv$part2[,2]-mp$mu2)		 # returns the estimated vector
return(retval)
}
	

em.putvalue <- function(pv,sv)
# Replaces missing values with predicted values
# pv - predicted vector
# sv - split vector
{
sv$part1[,2] <- pv		# internal function variable initialization

retval <- c(sv$part1[,2],sv$part2[,2])		# returns the complete vector
return(retval)
}
	
	
em.contribt2 <- function(sp,pv,sv)
# Estimates the contribution of NA vector to T2
# sp - partitioned covariance matrix
# pv - predicted vector
# sv - split vector
{
c1 <- sp$sigma11-sp$sigma12%*%solve(sp$sigma22)%*%sp$sigma21+pv%*%t(pv)		 # estimates missing X missing contribution to T2
c2 <- pv%*%t(sv$part2[,2])		# estimates missing X known contribution to T2

retval <- list(cont2na=c1,cont2obs=c2)		# returning values
return(retval)
}
	

em.arima <- function(xn,o,s,eps,maxit)
# Estimates one-step ahead ARIMA predictions for each time series of a given multivariate time series
#  09/07/2002 - Modified to support seasonal series 
# xn - data matrix
# o - order matrix
# s - seasonal period
# mit - optimizer max iterations
# eps - optimizer relative tolerance
{
rows <- dim(xn)[1]
cols <- dim(xn)[2]
models <- as.list(rep(NA,cols))

ar.pred <- matrix(NA,nrow=rows,ncol=cols)		 # internal function variable initialization
for (j in 1:cols)		# applies models to all series
	{
	if (is.null(s))
		{
		order <- o[1:3,j]	# assigns jth order vector to jth series - without seasonality
		seasonal <- list(order=c(0,0,0),period=NA)
		}
	else
		{
		order <- o[1:3,j]
		seasonal <- list(order=o[4:6,j],period=s)		# assigns jth order vector to jth series - with 		seasonality
		}
			
	models[[j]] <- arima(xn[,j],order=order,seasonal=seasonal,xreg=NULL,optim.control=list(maxit=maxit,reltol=eps))	# fits 	arima model to jth series
	#ar.filter <- predict(models[[j]],xreg=NULL)		# applies filter to jth series
	#ar.pred[,j] <- ar.filter$pred		# stores filtered series to internal variable
	ar.pred[,j] <- xn[,j] - residuals(models[[j]]) # compute the one-step ahead prediction using the inovation
	}

retval <- list(ar.pred=ar.pred,models=models)	 # returning value
return(retval)
}


em.spline <- function(xn,df,w)
# Estimates smooth spline for each time series of a given multivariate time series
# xn - data matrix
# df - degrees of freedom for the splines. If NULL, cross-validation is used 
# w - weights matrix
{
rows <- dim(xn)[1]
cols <- dim(xn)[2]

t <- seq(1:rows)	# time index for smoothing
if (length(df) == 1)
df <- rep(df,cols)	# extend df vector if length=1
	
sp.pred <- matrix(NA,nrow=rows,ncol=cols)	 # internal function variable initialization
if (is.null(w))
	w <- matrix(1,nrow=rows,ncol=cols)	# creates null weights matrix
		
models <- as.list(rep(NA,cols))
if (is.null(df))
	for (j in 1:cols)		# applies models to all series
		{
		models[[j]] <- smooth.spline(x=t,y=xn[,j],w=w[,j],cv=TRUE)
		sp.pred[,j] <- predict(models[[j]],t)$y
		}
else
	for (j in 1:cols)		# applies models to all series
		{
		models[[j]] <- smooth.spline(x=t,y=xn[,j],w=w[,j],df=df[j])
		sp.pred[,j] <- predict(models[[j]],t)$y
		}
			
retval <- list(sp.pred=sp.pred,models=models)	 # returning value
return(retval)
}

	
em.gam <- function(formula,xn,names,dataset,w,eps,maxit,bf.eps,bf.maxit)
# Estimates a gam model for each time series of a given multivariate time series
# formula - gam formula
# xn - dependent variables matrix
# names - variable names
# dataset - dataset for covariates
# w - weights matrix
{
rows <- dim(xn)[1]
cols <- dim(xn)[2]

t <- seq(1:rows)	# time index for smoothing
ga.pred <- matrix(NA,nrow=rows,ncol=cols)	 # internal function variable initialization
if (is.null(w))
	w <- rep(1,rows)	# creates null weights 	vector
	
models <- as.list(rep(NA,cols))
XN <- as.data.frame(xn)
dimnames(XN) <- names
dataset <- cbind.data.frame(xn,dataset)
formula <- paste("XN$",formula,sep="")
for (j in 1:cols)		# applies models to all series
	{
	models[[j]] <- 	gam(as.formula(formula[j]),family=gaussian(),data=dataset,weights=w,na.action=na.exclude,epsilon=eps,maxit=maxit,bf.epsilon=bf.eps,bf.maxit=bf.maxit)
	ga.pred[,j] <- fitted(models[[j]])
	}

retval <- list(ga.pred=ga.pred,models=models)	# returning value
return(retval)
}


em.filter <- function(XN,names,dataset,method,ar.control,sp.control,ga.control,f.eps,f.maxit,ga.bf.eps,ga.bf.maxit)
# Performs filtering on the time series
{
if (method == "arima")
	{
	m.obj <- em.arima(xn=XN,o=ar.control$order,s=ar.control$period,eps=f.eps,maxit=f.maxit)		# then predicts one-step ahead arima estimates
	MA <- m.obj$ar.pred
	}
else if (method == "spline")
	{
	m.obj <- em.spline(xn=XN,df=sp.control$df,w=sp.control$weights)		# then predicts smooth spline estimates
	MA <- m.obj$sp.pred
	}
else if (method == "gam")
	{
	m.obj <- em.gam(formula=ga.control$formula,xn=XN,names=names,dataset=dataset,w=ga.control$weights,eps=f.eps,maxit=f.maxit,bf.eps=ga.bf.eps,bf.maxit=ga.bf.maxit)	# gam estimates
	MA <- m.obj$ga.pred
	}
else
	stop("Filtering method not implemented.")

retval <- list(MA=MA,models=m.obj$models)
return(retval)
}


em.nofilter <- function(M,by,rows,cols)
# Create mean matrix when it is not time series
{
nbylevels <- nlevels(by)
MA <- matrix(NA,nrow=rows,ncol=cols)
for (i in 1:rows)
	MA[i,] <- M[,,by[i]]
	
retval <- list(MA=MA)
return(retval)
}


em.recursion <- function(X,XN,MA,S,C,rows,cols,by,nbylevels,n)
# Performs algorithm recursion over the observations
{
weights <- matrix(1.0,rows,cols)
SCT1 <- array(0.0,dim=c(1,cols,nbylevels))	# internal function variables initialization
SCT2 <- array(0.0,dim=c(cols,cols,nbylevels))
		
for (i in 1:rows)		# applies to all vectors
	{
	if((C$nainrow[i,2] != 0) &&  (C$nainrow[i,2] != cols))	# if there is at least one missing 	value but not intire vector is missing (mincol is deprecated)
		{     # then runs E step
		SV <- em.splitvec(X[C$nainrow[i,1],])	# splits i-th vector at position nainrow
		MP <- em.partmu(em.arrangevec(MA[i,],SV),C$nainrow[i,2])	# assigns predicted vector
		SP <- em.partsigma(em.arrangemat(S[,,by[i]],SV),C$nainrow[i,2])		# assigns partitioned covariance matrix
		CT1 <- em.contribt1(MP,SP,SV)	# estimates contribution to T1
		TV <- em.putvalue(CT1,SV)	# puts estimated values in a temporary vector
		XN[C$nainrow[i,1],] <- em.rearrangevec(TV,SV)	# replaces existing vector in XN with the original order vector
		CT2 <- em.contribt2(SP,CT1,SV)		# estimates contributions to T2
		SCT1[,,by[i]] <- SCT1[,,by[i]]+XN[i,]	# increments summation vector of T1 in order to update mean vector
		SCT2INC <- em.rearrangemat(em.correctmat(TV%*%t(TV),CT2),SV)  # increment in SCT2, but matrix is ok
		SCT2[,,by[i]] <- SCT2[,,by[i]]+SCT2INC	# increment summation matrix of T2 in order to update covariance matrix
        }
	else if((C$nainrow[i,2] != 0) &&  (C$nainrow[i,2] == cols))	# if whole vector is missing (mincol is deprecated)
		{
		XN[C$nainrow[i,1],] <- MA[i,]
		SCT1[,,by[i]] <- SCT1[,,by[i]]+XN[i,] 	# \tilde{x}=\tilde{\mu}
		SCT2INC <- S[,,by[i]]+XN[i,]%*%t(XN[i,])	# \tilde{x}_{j}\tilde{x}_{j}^{'}=\tilde{\Sigma}+\tilde{\mu}_{j}\tilde{\mu}_{j}^{'}
		SCT2[,,by[i]] <- SCT2[,,by[i]]+SCT2INC
		}
	else	# else (no missing values)
		{
		SCT1[,,by[i]] <- SCT1[,,by[i]]+XN[i,]	# increments summation vector of T1 in order to update mean vector
		SCT2INC <- XN[i,]%*%t(XN[i,])
		SCT2[,,by[i]] <- SCT2[,,by[i]]+SCT2INC	# increments summation matrix of T2 in order to update covariance matrix
		}
    weights[i,] <- (cols-0.5*C$nainrow[i,2])/cols   # crude measure of uncertanty
	}
retval <- list(XN=XN,SCT1=SCT1,SCT2=SCT2,weights=weights)
return(retval)
}


mnimput <- function(formula,dataset,by=NULL,log=FALSE,log.offset=1,eps=1e-3,maxit=1e2,ts=TRUE,method="spline",sp.control=list(df=NULL,weights=NULL),ar.control=list(order=NULL,period=NULL),ga.control=list(formula,weights=NULL),f.eps=1e-6,f.maxit=1e3,ga.bf.eps=1e-6,ga.bf.maxit=1e3,verbose=FALSE,digits=getOption("digits"))
# EM algorithm with time series option
#  09/07/2002 - Modified to support multiplicative seasonal series
# dataset - input data matrix
# eps - stop criterion
# maxit - maximum iterations
# ts - is time series
# method - filtering method for time series imputation. it can be "smooth", "gam" or "arima"
# sp.control()
# df - degrees of freedom for splines
# weights - matrix holding weights for smooth spline
# ar.control()
# order - order matrix (required if time series is true)
# period - seasonal period (required if time series is true)
# ar.maxit - ARIMA optimizer max iterations (required if time series is true)
# ar.eps - ARIMA optimizer relative tolerance (required if time series is true)
# gam.control()
#
{
if (missing(formula))
	stop("Model formula is missing")
if (missing(dataset))
	stop("Dataset is missing")

t1 <- proc.time()	# start time
call <- match.call()	# stores call sintaxe
if (verbose)
	cat("Expectation-Maximization Algorithm\n")
mdframe <- model.frame(formula,dataset,na.action=na.pass)  # get missing data frame
names <- dimnames(mdframe)		# keeps dataset variables names

rows <- dim(mdframe)[1]		# internal function variables initialization
cols <- dim(mdframe)[2]
if (is.null(by))
	{
	by <- as.factor(rep(1,rows))
	nbylevels <- nlevels(by)
	}
else if (sum(is.na(by))>0)
		stop("BY indicator must not have missing values")
else
	{
	by <- as.factor(by)
	nbylevels <- nlevels(by)
	}

X <- matrix(unlist(mdframe),nrow=rows,ncol=cols)	 #copies input data matrix into X matrix to avoid conflicts with variables names

if (log)
   if (!is.null(log.offset))
      X <- log(X+log.offset)
   else
      stop("Log offset is set to NULL")

M <- array(sapply(by(X,by,em.mean),as.matrix),dim=c(cols,1,nbylevels))	 # assigns mean vector
dimnames(M) <- list(names[[2]],"mean",levels(as.factor(by)))
C <- em.extractcoord(X)		# stores missing values coordinates
XN <- em.replacewmean(X,M,C,FALSE)		# replaces jth column missing data in X matrix with jth column mean
#S <- em.dispersion(XN,TRUE)	# assigns the covariance matrix of XN
S <- array(sapply(by(XN,by,em.dispersion),as.matrix),dim=c(cols,cols,nbylevels))  # creates an array of covariance matrices
dimnames(S) <- list(names[[2]],names[[2]],levels(as.factor(by)))
cc <- 1e35		# initializes convergence criterion variable
#det <- em.det(S)		# assigns determinant of S
det <- 0
for (i in 1:nbylevels) 
	det <- det+em.det(S[,,i])  # sum determinants of all covariance matrices
k <- 0	
n <- c(by(by,by,length))
#m <- c(by(X,by,em.countmvec))		# assigns missing vector counts
models <- NULL

if ((method=="spline") && (length(sp.control$df)==1))
	sp.control$df <- rep(sp.control$df,cols)
	
if (verbose)
	cat("Iteration   Covergence Criterion   Elapsed Time (hh:mm:ss)\n")
	
while ((cc > eps) && (k < maxit))	# repeat until the convergence criterion or maximum iterations is reached
	{
	olddet <- det	# reassigns determinant of S
	if (ts)
		{		# if is time series
		filtered <- em.filter(XN,names,dataset,method,ar.control,sp.control,ga.control,f.eps,f.maxit,ga.bf.eps,ga.bf.maxit)
		MA <- filtered$MA
		models <- filtered$models
		}
	else
		{
		MA <- em.nofilter(M,by,rows,cols)$MA		# if it is not a time series problem, use the actual estimate of column means
		models <- NULL
		}
	recursed <- em.recursion(X,XN,MA,S,C,rows,cols,by,nbylevels,n)
	XN <- recursed$XN
	for (j in 1:nbylevels)
		M[,,j] <- recursed$SCT1[,,j]/n[j]	# estimates sufficient statistic T1 and computes the mean vector
	for (j in 1:nbylevels)
		S[,,j] <- recursed$SCT2[,,j]/n[j]-M[,,j]%*%t(M[,,j])	# estimates sufficient statistic T2 and computes the covariance matrix
	weights <- recursed$weights
	
	det <- 0
    for (i in 1:nbylevels) 
        det <- det+em.det(S[,,i])  # sum determinants of all covariance matrices
	cc <- abs((det-olddet)/olddet)		# updates covergence criterion
	k <- k+1	# increments iterations counter
	t2 <- proc.time()	# end time
	
	if (verbose)
		cat(paste(k,paste("         ",format(round(cc,digits),nsmall=digits),sep=""),elapsedtime(t1,t2),"\n", 		sep="          "))	# updates log file
	}
	
converged <- ifelse(cc<=eps,TRUE,FALSE)	
if (verbose)
	{
	if (converged)
		cat("\nConverged\n")	# updates log file
	else
		cat("\nDid not converged\n")	# updates log file
	}	
dimnames(MA) <- names	# gives variables names back
dimnames(XN) <- names	# gives variables names back
dimnames(weights) <- names
for (j in 1:nbylevels)
	names(M[,,j]) <- names[[2]]
for (j in 1:nbylevels)
	dimnames(S[,,j]) <- list(names[[2]],names[[2]])
	 
XN <- as.data.frame(XN)		# converts imputed dataset in dataframe
etime <- elapsedtime(t1,proc.time())	# getting elapsed time
missings <- C$nainrow[,2]
	
retval <- list(call=call,filled.dataset=XN,dataset=mdframe,muhat=M,sigmahat=S,level=MA,weights=weights,missings=missings,iterations=k,convergence=cc,converged=converged,time=etime,models=models,log=log,log.offset=log.offset)	 # assigns return values to temporary object
class(retval) <- "mtsdi"
return(retval)		# returning value
}


getmean <- function(object,weighted=TRUE,mincol=1,maxconsec=3)
# gets mean of row when mincol sastifies
{
if (!inherits(object,"mtsdi"))
    stop("Object class is not mtsdi")
     
if (object$log)
    filled.dataset <- exp(as.matrix(object$filled.dataset))-object$log.offset
else
    filled.dataset <- as.matrix(object$filled.dataset)
weights <- as.matrix(object$weights)
nm <- as.vector(object$missings)
n <- nrow(filled.dataset)
p <- ncol(filled.dataset)

misscount <- em.countmnear(object$dataset)
for(j in 1:p)
	for(i in 1:n)
		filled.dataset[i,j] <- ifelse(misscount[i,j]<=maxconsec,filled.dataset[i,j],NA)

cmean <- double(n)

for (t in 1:n)
    {
	if ((p-nm[t]) >= mincol)
		{
     	if (weighted)
        	cmean[t] <- weighted.mean(filled.dataset[t,],weights[t,],na.rm=TRUE)
     	else
        	cmean[t] <- mean(filled.dataset[t,],na.rm=TRUE)
     	}
	else
		cmean[t] <- NA
	}
return(cmean)
}


predict.mtsdi <- function(object,...)
# puts values on its original scale if log is used
{
if (!inherits(object,"mtsdi"))
    stop("Object class is not mtsdi")
     
if (object$log)
    filled.dataset <- exp(as.matrix(object$filled.dataset))-object$log.offset
else
    filled.dataset <- as.matrix(object$filled.dataset)
return(filled.dataset)
}


plot.mtsdi <- function(x,vars="all",overlay=TRUE,level=TRUE,points=FALSE,leg.loc="topright",horiz=FALSE,at.once=FALSE,...)
# plot method for imputed data plotting
{
if (!inherits(x,"mtsdi"))
    stop("The class of this object is not mtsdi.")
if(x$log)
	{
	ylab <- "Axis on the log scale"
	dataset <- log(x$dataset+x$log.offset)
	}
else
	{
	ylab <- "Axis on the original scale"
	dataset <- x$dataset
	}
filled.dataset <- as.matrix(x$filled.dataset)
llevel <- x$level

rows <- dim(dataset)[[2]]

if (vars[1]=="all")
	vars <- seq(1:rows)
else if ((max(vars)>rows) || (length(vars)>rows))
	stop("Argument vars must be \"all\" or a valid vector")

#if(typeof(leg.loc)!="list")	
#	if (leg.loc=="auto") 
#		{ 
#		leg.loc <- list(x=18.5,y=16.1)
#		horiz <- FALSE
#		}

for (i in vars)
	{
	if (at.once)
		get(getOption("device"))()
	else
		if (interactive() && (i!=vars[1]) && (length(vars)!=1))
			{
			prompt <- readline("Press ENTER for next page or X to exit and keep this page... ")
			if((prompt == "x") || (prompt == "X"))
				break
			}
	plot(filled.dataset[,i],ylab=ylab,type="l",col="red",...)
	if (overlay)
		lines(dataset[,i],col="black",...)
	if (level)
		lines(llevel[,i],col="blue",...)
	if (points)
		points(filled.dataset[,i],col="red",...)
	legend(x=leg.loc,legend=c("observed","imputed","level"),col=c("black","red","blue"),lty=1,horiz=horiz) # alternate location c(18.5,16.1)
	title(main=paste("Imputed data for series",colnames(filled.dataset[i])))
	}
}


summary.mtsdi <- function(object,...)
# summary
{
call <- object$call
muhat <- object$muhat
sigmahat <- object$sigmahat
iterations <- object$iterations
convergence <- object$convergence
converged <- object$converged
time <- object$time
models <- object$models
log <- object$log
log.offset <- object$log.offset

retval <- list(call=call,muhat=muhat,sigmahat=sigmahat,iterations=iterations,convergence=convergence,converged=converged,time=time,models=models,log=log,log.offset=log.offset)	 # assigns return values to temporary object
class(retval) <- "summary.mtsdi"
return(retval)
}


print.summary.mtsdi <- function(x,digits=getOption("digits"),print.models=TRUE,...)
# summary print method
{
cat("\nCall:\n")
print(x$call)
cat("\nEstimated mean vector:\n")
for (i in 1:dim(x$muhat)[3])
	{
	if(dim(x$muhat)[3]!=1)
		cat("\nBY factor level: ",dimnames(x$muhat)[[3]][i],"\n",sep="")
	print(x$muhat[,,i])
	}
cat("\nEstimated covariance matrix:\n")
for (i in 1:dim(x$muhat)[3])
	{
	if(dim(x$muhat)[3]!=1)
		cat("\nBY factor level: ",dimnames(x$muhat)[[3]][i],"\n",sep="")
	print(x$sigmahat[,,i])
	}

if (x$log)
	cat("\nData are on the log scale with an offset of ",x$log.offset,".\n",sep="")
else
	cat("\nData are on the original scale.\n")

conv <- ifelse(x$converged==TRUE,"converged","did not converge")
cat("\nThe algorithm ",conv," after ",x$iterations," iterations with relative diference in covariance matrix equal to ",round(x$convergence,digits),".\n",sep="")
cat("The process took ",x$time,".\n",sep="")

if ((print.models)&&(!is.null(x$models)))
	{
	cat("\nTime filtering models:")
	for (i in 1:length(x$models))
		{
		cat("\n\nFilter model of variable: ",dimnames(x$muhat)[[1]][i],"\n",sep="")
		print(x$models[[i]])
		}
	}
}


print.mtsdi <- function(x,digits=getOption("digits"),...)
# print method for mtsdi class
{
cat("\nCall:\n")
print(x$call)
if (x$log)
	cat("\nData are on the log scale with an offset of ",x$log.offset,".\n",sep="")
else
	cat("\nData are on the original scale.\n")

conv <- ifelse(x$converged==TRUE,"converged","not converged")
cat("\nThe algorithm ",conv," after ",x$iterations," iterations with relative diference in covariance matrix equal to ",round(x$convergence,digits),".\n",sep="")
cat("The process took ",x$time,".\n",sep="")
}




