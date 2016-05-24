#<<BEGIN>>
dempiricalD <- function(x,values,prob=NULL,log=FALSE)
#TITLE The Discrete Empirical Distribution
#NAME empiricalD
#KEYWORDS distribution
#DESCRIPTION
#Density, distribution function and random generation
#for a discrete empirical distribution. This function is vectorized to accept
#different sets of \samp{values} or \samp{prob}.
#INPUTS
#{x, q}<<Vector of quantiles.>>
#{p}<<Vector of probabilities.>>
#{n}<<Number of random values. If length(n) > 1, the length is taken to be the number required.>>
#{values}<<Vector or matrix of numerical values. See details.>>
#[INPUTS]
#{prob}<<Optionnal vector or matrix of count or probabilities. See details.>>
#{log, log.p}<<logical; if \samp{TRUE}, probabilities \samp{p} are given as \samp{log(p)}.>>
#{lower.tail}<<logical; if \samp{TRUE} (default), probabilities are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.>>
#DETAILS
#If \samp{prob} is missing, the discrete distribution is obtained directly from the vector of \samp{values},
#otherwise \samp{prob} is used to weight the values. \samp{prob} is normalized before use. Thus, \samp{prob}
#may be the count of each \samp{values}. \samp{prob} values should be non negative and their sum should not be 0.</>
#\samp{values} and/or \samp{prob} may vary: in that case, \samp{values} and/or \samp{prob} should be sent
#as matrixes, the first row being used for the first element of \samp{x}, \samp{q}, \samp{p} or the first random value, the
#second row for the second element of \samp{x}, \samp{q}, \samp{p} or random value, ...
#Recycling is permitted if the number of rows of \samp{prob} and \samp{values} are equal or
#if the number of rows of \samp{prob} and/or \samp{values} are one.</>
#\samp{rempiricalD(n, values, prob)} with \samp{values} and \samp{prob} as vectors
#is equivalent to \samp{sample(x=values, size=n, replace=TRUE, prob=prob)}.
#NOTE
#In the future, the fonctions should be written for non numerical values.
#SEE ALSO
#\code{\link{sample}}.
#\code{\link{empiricalC}}.
#VALUE
#\samp{dempiricalD} gives the density, \samp{pempiricalD} gives the distribution function,
#\samp{qempiricalD} gives the quantile function and \samp{rempiricalD} generates random deviates.
#EXAMPLE
#dempiricalD(1:6,2:6,prob=c(10,10,70,0,10))
#pempiricalD(1:6,2:6,prob=c(10,10,70,0,10))
#qempiricalD(seq(0,1,0.1),2:6,prob=c(10,10,70,0,10))
#table(rempiricalD(10000,2:6,prob=c(10,10,70,0,10)))
#
### Varying values
#(values <- matrix(1:10,ncol=5))
### the first x apply to the first row : p = 0.2
### the second x to the second one: p = 0
#dempiricalD(c(1,1),values)
#
#
###Use with mc2d
###Non Parameteric Bootstrap
#val <- c(100, 150, 170, 200)
#pr <- c(6, 12, 6, 6)
#out <- c("min", "mean", "max")
###First Bootstrap in the uncertainty dimension
#(x <- mcstoc(rempiricalD, type = "U", outm = out, nvariates = 30, values = val, prob = pr))
###Second one in the variability dimension
#mcstoc(rempiricalD, type = "VU", values = x)
#CREATED 08-06-15
#--------------------------------------------
{
  if(length(x) == 0) return(numeric(0))
  if(is.vector(values)) values <- matrix(values,nrow=1)
  if(is.null(prob)) prob <- array(1,dim=dim(values))
  if(is.vector(prob)) prob <- matrix(prob,nrow=1)

  if(!is.numeric(values) || !is.numeric(prob) || !is.numeric(x)) stop("Non numeric arguments in a mathematical function")

  if(any(is.na(prob)) || any(apply(prob,1,"<",0)) || any(rowSums(prob) == 0))
    stop("Nas in prob, or sum(prob)=0 or negative values of prob")

  ncv <- ncol(values)
  ncp <- ncol(prob)
  if(ncv != ncp) stop("prob and values should be of same length or have the same number of columns.")

  l <- length(x)
  nrv <- nrow(values)
  nrp <- nrow(prob)
  mnr <- min(max(nrv,nrp),l)
  if(nrv != 1 && nrp != 1 && nrv != nrp)
    stop("values/prob should be vector(s), matrix(es) of 1 row or matrix(es) of the same number of rows")

  prob   <- lapply(1:min(nrp,mnr), function(x) prob[x,])
  values <- lapply(1:min(nrv,mnr), function(x) values[x,])

  prob <- mapply(function(pr,val) tapply(pr,val,sum), prob, values, SIMPLIFY=FALSE)
  prob <- lapply(prob, function(y) y /sum(y))
  values <- lapply(values,function(y) sort(unique(y)))

  find <- Vectorize(function(x,val,pr) {
    if(is.na(x)) return(NA)
    quel <- x==val
    if(any(quel)) return(pr[quel])
    return(0)}
    ,vectorize.args = "x")

  # to gain time and memory
  n1 <- l %/% mnr
  n2 <- l %% mnr
  l2 <- l-n2
  x1 <- x[1:l2]

  res <- vector(mode = "numeric", length = l)
  res[1:l2] <- mapply(find,x1,values,prob)
  if(n2!=0){
    x2 <- x[(l2+1):l]
    res[(l2+1):l] <- mapply(find,x2,values[1:n2],prob[1:n2])
    }

  if(any(is.na(res))) warning("NaN in dempiricalD")
  if(log) res <- log(res)
  return(res)}

#<<BEGIN>>
pempiricalD <- function(q,values,prob=NULL,lower.tail = TRUE, log.p = FALSE)
#ISALIAS dempiricalD
#--------------------------------------------
{
  if(length(q) == 0) return(numeric(0))
  if(is.vector(values)) values <- matrix(values,nrow=1)
  if(is.null(prob)) prob <- array(1,dim=dim(values))
  if(is.vector(prob)) prob <- matrix(prob,nrow=1)

  if(!is.numeric(values) || !is.numeric(prob) || !is.numeric(q)) stop("Non numeric arguments in a mathematical function")

  if(any(apply(prob,1,is.na)) || any(apply(prob,1,function(x) x < 0)) || any(rowSums(prob)==0))
    stop("Nas in prob, or sum(prob)=0 or negative values of prob")

  ncv <- ncol(values)
  ncp <- ncol(prob)
  if(ncv != ncp) stop("prob and values should be of same length or have the same number of columns.")

  l <- length(q)
  nrv <- nrow(values)
  nrp <- nrow(prob)
  mnr <- min(max(nrv,nrp),l)
  if(nrv != 1 && nrp != 1 && nrv != nrp)
    stop("values/prob should be vector(s), matrix(es) of 1 row or matrix(es) of the same number of rows")

  values <- lapply(1:min(nrv,mnr),function(y) values[y,])
  prob <- lapply(1:min(nrp,mnr),function(y) prob[y,])

  prob <- mapply(function(pr,val) tapply(pr,val,sum),prob,values,SIMPLIFY=FALSE)
  prob <- lapply(prob, function(x) cumsum(x)/sum(x))
  values <- lapply(values,function(y) sort(unique(y)))

  # the "-1" combined with the "<" allows duplicated values
  find <- Vectorize(function(x,val,pr){
        if(is.na(x)) return(NA)
        if(x < val[1]) return(0)
        return(max(pr[x >= val]))},vectorize.args = "x")

  # to gain time and memory
  n1 <- l %/% mnr
  n2 <- l %% mnr
  l2 <- l-n2
  q1 <- q[1:l2]

  res <- vector("numeric",length=l)
  res[1:l2] <- mapply(find,q1,values,prob)
  if(n2!=0){
    q2 <- q[(l2+1):l]
    res[(l2+1):l] <- mapply(find,q2,values[1:n2],prob[1:n2])
    }

  if(any(is.na(res))) warning("NaN in pempiricalD")
  if(!lower.tail) res <- 1-res
  if(log.p) res <- log(res)
  return(res)}


#<<BEGIN>>
qempiricalD <- function(p,values,prob=NULL,lower.tail = TRUE, log.p = FALSE)
#ISALIAS dempiricalD
#--------------------------------------------
{
  if(length(p) == 0) return(numeric(0))
  if(is.vector(values)) values <- matrix(values,nrow=1)
  if(is.null(prob)) prob <- array(1,dim=dim(values))
  if(is.vector(prob)) prob <- matrix(prob,nrow=1)

  if(!is.numeric(values) || !is.numeric(prob) || !is.numeric(p)) stop("Non numeric arguments in a mathematical function")

  if(any(apply(prob,1,is.na)) || any(apply(prob,1,function(x) x < 0)) || any(rowSums(prob)==0))
    stop("Nas in prob, or sum(prob)=0 or negative values of prob")

  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p

  ncv <- ncol(values)
  ncp <- ncol(prob)
  if(ncv != ncp) stop("prob and values should be of same length or have the same number of columns.")

  l <- length(p)
  nrv <- nrow(values)
  nrp <- nrow(prob)
  mnr <- min(max(nrv,nrp),l)
  if(nrv != 1 && nrp != 1 && nrv != nrp)
    stop("values/prob should be vector(s), matrix(es) of 1 row or matrix(es) of the same number of rows")

  values <- lapply(1:min(nrv,mnr),function(y) values[y,])
  prob <- lapply(1:min(nrp,mnr),function(y) prob[y,])

  prob <- mapply(function(pr,val) tapply(pr,val,sum),prob,values,SIMPLIFY=FALSE)
  prob <- lapply(prob, function(x) cumsum(x)/sum(x))
  values <- lapply(values,function(y) sort(unique(y)))

  find <- Vectorize(function(x,val,pr){
    if(x > 1 || x < 0) return(NaN)
    return(min(val[x <= pr]))},vectorize.args = "x")

  # to gain time and memory
  n1 <- l %/% mnr
  n2 <- l %% mnr
  l2 <- l-n2
  p1 <- p[1:l2]

  res <- vector("numeric",length=l)
  res[1:l2] <- mapply(find,p1,values,prob)
  if(n2!=0){
    p2 <- p[(l2+1):l]
    res[(l2+1):l] <- mapply(find,p2,values[1:n2],prob[1:n2])
    }

  if(any(is.na(res))) warning("NaN in qempiricalD")
  return(res)
}

#<<BEGIN>>
rempiricalD <- function(n,values,prob=NULL)
#ISALIAS dempiricalD
#--------------------------------------------
{
  if(length(n) > 1) n <- length(n)
  if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  n <- as.integer(n)
  if(n < 0) stop("integer(n) cannot be negative in rempiricalD")
 
  if(is.vector(values)) values <- matrix(values,nrow=1)
  if(is.null(prob))     prob <- rep(1,length(values[1,]))
  if(is.vector(prob))   prob <- matrix(prob,nrow=1)

  if(!is.numeric(values) || !is.numeric(prob) || !is.numeric(n)) stop("Non numeric arguments in a mathematical function")

  ncv <- ncol(values)
  ncp <- ncol(prob)
  if(ncv!=ncp) stop("Prob and values should be of same length or have the same number of columns.")

  nrv <- nrow(values)
  nrp <- nrow(prob)
  mnr <- min(max(nrv,nrp),n)
  if(nrv!= 1 && nrp!=1 && nrv!=nrp) stop("values/prob should be vector(s), matrix(es) of 1 row or matrix(es) of the same number of rows")

  values <- lapply(1:min(nrv,mnr), function(y) values[y,])
  prob <- lapply(1:min(nrp,mnr), function(y) prob[y,])

  n1 <- n %/% mnr
  n2 <- n%% mnr
  # To gain time, send sample on n1
  res <- mapply(sample, x=values, prob=prob, MoreArgs=list(size=n1,replace=TRUE))
  dim(res) <- c(n1,mnr)
  res <- t(res)
  dim(res) <- NULL
  # then the remaining n2, if any
  if(n2!=0)
    res <- c(res,mapply(sample,x=values[1:n2],prob=prob[1:n2],MoreArgs=list(size=1,replace=TRUE)))
  return(as.vector(res))
}

