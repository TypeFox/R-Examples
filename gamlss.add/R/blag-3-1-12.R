#----------------------------------------------------------
# a function to create lags 
# author MS
# last change 14-8-11 
#----------------------------------------------------------
#----------------------------------------------------------
# this functions creates a matrix of lag basis for the vector x
# lag : how many lag values (that is the number of colums in the matrix if x is not included)
# omit.na : if true the first "lag" rows of the resulting matrix are omited
# no.x : If TRUE only lag values are included if FALSE  x is included in the matrix b as the first colunm
# value : what values should be set at the begining of the lags colunms, by default is set to NA
# but this has to changed to something else like zero if the lag basis is to be used withing GAMLSS.
# In this case the wlag and llag functions are usefull
# llag produdes a list with a matrix and weights while wlag takes the basis for lag matrix and produdes 
# the prior weights to be used in the regression
#------------------------------------------------------------ 
blag <- function(x, lags=1, from.lag=0, omit.na = FALSE,  value=NA, ...)
{
## local function 	
 lag1 <- function(x)
   {
    l <- c(value, x[1:(length(x)-1)])
    l
   }
  ## main function starts here
  lenx <- length(x)
 start <- if (from.lag==0) 1 else from.lag
 ncolX <- lags+start
     d <- matrix(0, nrow=length(x), ncol=ncolX)
 d[,1] <- x
 xname <- deparse(substitute(x))
 cname <- xname
 #if (lag==1) 
 # {names(d[,1]) <- c("lag1")
 # if (omit.na) d <-na.omit(d)
 # return(d)
 # }
for(i in 2:ncolX) 
 {
 d[,i] <- lag1(d[,i - 1])
 cname <- c(cname,paste(xname, paste("lag",i-1, sep=""), sep=""))
 }
colnames(d ) <- cname
#d <- as.matrix(d, nrows=length(x))
if (omit.na) d <-na.omit(d)
 start1 <- if (from.lag==0) 1 else from.lag+1
 d <- d[,(start1:ncolX)]
class(d) <- c("mlags", "matrix")
attr(d, "no.lags") <- lags
d
 }
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
llag<-function(x, ...) 
{
L <- blag(x, ...)
W <- rep(1, length(x))
W[1:attr(L,"no.lags")] <- 0	
list(matrix=L, weights=W)
}
#-------------------------------------------------------------------
#-------------------------------------------------------------------
wlag <- function(x, lags=NULL)
{
if (is(x, "mlags")) 
  { #stop("the object need to be an mlags object")
     W <- if (is.matrix(x)) rep(1, dim(x)[1]) else rep(1, length(x))
     W[1:attr(x,"no.lags")] <- 0	
  }
else 
 {
  W <- if (is.null(lags)) stop("the lag argument is needed here")
        else rep(1, length(x))
   W[1:lags] <- 0
 }
 W
}
#-------------------------------------------------------------------
#-------------------------------------------------------------------
# Example 
#y <- arima.sim(500, model=list(ar=c(.4,.3,.1))) 
#X <- blag(y, 5, value=0)
#head(X)
#w<-wlag(X)
#m1<-gamlss(y~X, weights=w )
#summary(m1)
