#########################################################################
# AUXILIAR FUNCTIONS
#########################################################################

#########################################################################

# FUNCTION:  EPANECHNIKOV KERNEL 

quad <- function(u)  
{  
    if ( abs(u)  > 1 ) {v <-0} else { v <- 0.75*(1 - (u)^2) }
    return(v)
}

#########################################################################

# FUNCTION: THE LOCAL CONSTANT (NADARAYA-WATSON)

# From a bivariate sample (data) and a sequence of bandwidths (h.seq), this 
# procedure estimates the regression fuction in the values of a explanatory 
# variable (let's say x).

# data ==> sample : data[,1] contains the scalar response yi,
#                   data[,2] contains the scalar explanatory variable ti
# x: values of the explanatory variable where the estimates will be computed
# h.seq: sequence of bandwidths 
# kind.of.kernel ==> allows us to choose between quadratic ("quad") or 
#                    triangular ("triangular") kernel

nadaraya.watson <- function(data, x, h.seq, kind.of.kernel)
    
{
    kernel <- get(kind.of.kernel)
    n <- nrow(data)
    t <- data[,1]
    y <- data[,2] 
    num.band <- length(h.seq)
    Yhat<-matrix(0, length(x), num.band)
    
    if ((sum(y)==Inf) | (sum(y)=="NaN")) 
    {
        Yhat[,]<-0/0
        return(Yhat)
    } else {          
        for(i in 1:length(x)) 
        {
            # We obtain, for each value in x and each bandwidth in h.seq, 
            # the estimate of E(y|x) using the Nadaraya-Watson estimator
            diff <- t-x[i]                              
            Zmat <- matrix(rep(diff, num.band), nrow = num.band, byrow = T)
            Umat <- Zmat/h.seq
            Kmat <- kernel(Umat)
            Kmat[(Kmat<0)] <- 0
            S0 <- apply(Kmat, 1, sum)
            Kmat <- Kmat/S0
            Ymat <- matrix(rep(y, num.band), nrow = num.band, byrow = T)
            Yhat[i,] <- apply(Ymat * Kmat, 1, sum)
        }
    }
    return(Yhat)
}

#########################################################################

# FUNCTION:  CROSS-VALIDATION BANDWIDTH SELECTOR FOR NADARAYA-WATSON ESTIMATOR 

# data (data matrix), r (2r+1 observations around each point will be left out), 
# h.seq (sequence of bandwidths where cross-validation function will be evaluated)

h_cv <- function(data,r,h.seq)
    
{ 
    n <- nrow(data)
    t <- data[,1]
    y <- data[,2] 
    
    mat <- array(1,dim=c(n,n))
    for (i in 1:n)   mat[i,max(i-r,1):min(i+r,n)] <-0
    cv.t <- array(t,dim=c(n,n))*mat
    cv.y <- array(y,dim=c(n,n))*mat
    
    h.cv.mat <- matrix(0, n, length(h.seq))
    h.cv.vec <- rep(0,length(h.seq))
    for (i in 1:n) 
    {
        dat <- cbind(cv.t[,i],cv.y[,i])
        h.cv.mat[i,] <- ( y[i]-nadaraya.watson(dat,t[i],h.seq,kind.of.kernel = "dnorm") )^2
    }
    h.cv.vec <- apply(h.cv.mat,2,sum)/n
    h.cv.min <- h.seq[which.min(h.cv.vec)]
    return(h.cv.min)
}

#########################################################################

# FUNCTION: THE MODIFIED LOCAL CONSTANT (NADARAYA-WATSON)

# As nadaraya.watson function but incorporating a truncation argument     

nadaraya.watson.mod <- function(data, x, h.seq, kind.of.kernel, CM)
    
{
    
    #num.band <- length(h.seq)
    #Yhat<-matrix(0, length(x), num.band)
    #for (i in 1:length(h.seq)) {
    #	nadarw =ksmooth( data[,1], data[,2], "normal", h.seq[i], x.points=x)$y
    #nadarw[nadarw > CM] <- CM
    #nadarw[nadarw < -CM] <- -CM
    #	Yhat[,i] = nadarw
    #}
    #return(Yhat)
    
    kernel <- get(kind.of.kernel)
    n <- nrow(data)
    t <- data[,1]
    y <- data[,2] 
    num.band <- length(h.seq)
    Yhat<-matrix(0, length(x), num.band)
    estado<-matrix(0, length(x), num.band)
    
    if ((sum(y)==Inf) | (sum(y)=="NaN")) 
    {
        Yhat[,]<-0/0
        return(Yhat)
    } else {          
        for(i in 1:length(x)) 
        {
            # We obtain, for each value in x and each bandwidth in h.seq, 
            # the estimate of E(y|x) using the Nadaraya-Watson estimator
            diff <- t-x[i]                              
            Zmat <- matrix(rep(diff, num.band), nrow = num.band, byrow = T)
            Umat <- Zmat/h.seq
            Kmat <- kernel(Umat)
            Kmat[(Kmat<0)] <- 0
            S0 <- apply(Kmat, 1, sum)
            Kmat <- Kmat/S0
            Ymat <- matrix(rep(y, num.band), nrow = num.band, byrow = T)
            SD <- apply(Kmat, 1, sum)
            SN <- apply(Ymat * Kmat, 1, sum)
            Yhat[i,] <- SN/SD
            for (j in 1:num.band) 
            {
                if ((Yhat[i,j]==Inf) | (Yhat[i,j]=="NaN")) 
                {
                    if (is.nan(SN[j]) | is.na(SN[j])) SN[j] <- 0
                    if (SN[j]>0) Yhat[i,j] <- CM
                    if (SN[j]==0) Yhat[i,j] <- 0
                    if (SN[j]<0) Yhat[i,j] <- -1*CM
                }
                if ( Yhat[i,j] > CM ) Yhat[i,j] <- CM
                if ( Yhat[i,j] < -1*CM ) Yhat[i,j] <- -1*CM
            }
        }
    }
    return(Yhat)
}

#########################################################################

# FUNCTION: BOOTSTRAP RESAMPLES OF THE RESIDUALS

bootstrap.res <- function(centred.res, band, kind.of.kernel="dnorm", resamples.number, tamano)
    
    # This function provides "resamples.number" bootstrap resamples of i.i.d. 
    # observations from the estimated density of the centred residuals. 
    
    # "centred.res" is the sequence of centred residuals.
    # The bandwidth and the kernel function used to estimate the density are "band" 
    # and "kind.of.kernel", respectively.
    # "resamples.number" is the number of bootstrap resamples required.
    # "tamano" is the size of each resample. 
    
{
    kernel <- get(kind.of.kernel)
    l <- tamano * resamples.number
    boot.res <- numeric(length= l)
    boot.res.mat <- matrix(0, tamano, resamples.number)
    boot.res <- sample(centred.res, size=l, replace=TRUE) + ( band * rnorm(l) )
    boot.res.mat <- matrix(boot.res, tamano, resamples.number)
    return(boot.res.mat)
}

#########################################################################

# FUNCTION: KERNEL DENSITY

estimator.density <- function( data, x, h.seq, kind.of.kernel="dnorm")
    
    # From a sample "data" and a sequence of bandwidths "h.seq", this procedure 
    # computes the kernel density function evaluated on s set of points "x".
    
{
    
    result <- approx(data[,1], data[,2], x)$y
    result[ is.na(result)] <- 0
    result
    
    #kernel <- get(kind.of.kernel)
    #n <- length(data)
    #num.band <- length(h.seq)
    #f.hat <- matrix(0, length(x), num.band)
    #for(i in 1:length(x)) 
    #	{
    #	diff <- data - x[i]                              
    #	Zmat <- matrix(rep(diff, num.band), nrow = num.band, byrow = T)
    #	Umat <- Zmat/h.seq
    #	Kmat <- kernel(Umat)
    #	Kmat[(Kmat<0)] <- 0
    #	f.hat[i,] <- (apply(Kmat, 1, sum))/(n*h.seq)
    #	}
    #return(f.hat)
}

#########################################################################


