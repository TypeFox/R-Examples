setClass("board",
         representation = representation(x="matrix"),
         prototype      = list(x=matrix())
)

setValidity("board",
            function(object){
              x <- object@x
              non.nas <- x[!is.na(x)]
              if (!is.matrix(x)){
                return("not a matrix")
              } else if (any(non.nas != round(non.nas))){
                return("matrix includes non-integers")
              } else if (any(non.nas < 0)){
                return("matrix includes negative numbers")
              } else if (any(is.nan(x))){
                return("matrix includes one or more NaN elements")
              } else if (dof(x)<1){
                return("less than one (1) degree of freedom")
              } else {
                return(TRUE)
              }
            }
            )

"is.board" <- function(x){is(x,"board")}

"as.board" <- function(x){
  if(is.board(x)){
    return(x)
  } else {
    return(new("board",x=x))
  }
}

"marginals" <- function(x){
  x <- as.board(x)@x
  mode(x) <- "integer"
  list(
       rs = apply(x,1,sum,na.rm=TRUE),
       cs = apply(x,2,sum,na.rm=TRUE),
       na = which(is.na(x),arr.ind=TRUE),
       x  = x
       )
}

".Cargs" <- function(x){
  jj <- marginals(x)
  list(
       rs   = as.integer(       jj$rs ),
       nrow = as.integer(length(jj$rs)),
       cs   = as.integer(       jj$cs ),
       ncol = as.integer(length(jj$cs)),
       na   = as.integer(       jj$na ),
       nna  = as.integer(nrow(  jj$na))
       )
} 

"aylmer.test" <- 
  function (x,  alternative = "two.sided", simulate.p.value = FALSE, n = 1e5, B = 2000, burnin=100, use.brob=FALSE) 
{
  DNAME <- deparse(substitute(x))
  if(is.function(alternative)){
    return(aylmer.function(x=x,  func=alternative, simulate.p.value = simulate.p.value, n = n, B = B, burnin=burnin, use.brob=use.brob, DNAME=DNAME))
  }
  METHOD <- "Aylmer test for count data"
  if(!any(is.na(x))){
    warning("supplied matrix has no NAs. Consider using 'stats:::fisher.test()'")
  }
  x.dof <- dof(x)
  stopifnot(x.dof>0)
  almost.1 <- 1 + 64 * .Machine$double.eps
  if(simulate.p.value){
    stopifnot(identical(length(B),1L))
    STATISTIC <- prob(x, use.brob=use.brob)  
    METHOD <-
      paste(METHOD, "with simulated p-value\n\t (based on", B, "replicates)")
    random_probs <- randomprobs(x, B, burnin=burnin, use.brob=use.brob)
    PVAL <-
      as.numeric((1+sum(random_probs <=  STATISTIC/almost.1))/(B+1))
  } else {
    STATISTIC <- prob(x, use.brob=use.brob, give.log=FALSE)
    a <- allprobs(x, n=n, normalize=FALSE)
    PVAL <- sum(a[a <= STATISTIC*almost.1])/sum(a)
    if(x.dof == 1){
      alternative <-
        char.expand(alternative, c("two.sided", "less", "greater"))
      PVAL <- 
        switch(alternative,
               two.sided = PVAL,
               greater   = .pval.1dof(x, greater=TRUE),
               less      = .pval.1dof(x, greater=FALSE)
               )
    }
  }
  RVAL <- list(p.value = PVAL, alternative = alternative, method = METHOD, 
               data.name = DNAME)
  attr(RVAL, "class") <- "htest"
  return(RVAL)
}

"aylmer.function" <-
  function (x, func, simulate.p.value = FALSE, n = 1e5, B = 2000, burnin=100, use.brob=FALSE, DNAME=NULL) 
{
  if(is.null(DNAME)){
    DNAME <- deparse(substitute(x))
  }
  METHOD <- "Aylmer functional test for count data"
  stopifnot(dof(x)>0)
  if(simulate.p.value){ # Monte Carlo ...
    stopifnot(identical(length(B),1L))
    STATISTIC <- func(x)  
    METHOD <-
      paste(METHOD, "with simulated p-value\n\t (based on", B, "replicates)")
    random_probs <- randomprobs(x, B, burnin=burnin, use.brob=use.brob, func=func)
    almost.1 <- 1 + 64 * .Machine$double.eps
    PVAL <-
      ##  as.numeric((1+sum(random_probs <=  STATISTIC/almost.1))/(B+1))
      as.numeric((1+sum(random_probs >=  STATISTIC*almost.1))/(B+1))
  } else {  # ... enumeration
    STATISTIC <- func(x)
    a <- allprobs(x, n=n, normalize=FALSE)
    allfuncs <- apply(allboards(x,n=n),3,func)
    ## PVAL <- sum(a[allfuncs <= STATISTIC])/sum(a)
    PVAL <- sum(a[allfuncs >= STATISTIC])/sum(a)
  }

  RVAL <- list(p.value = PVAL, alternative = "test function exceeds observed",
               method = METHOD, data.name = DNAME)
  attr(RVAL, "class") <- "htest"
  return(RVAL)
}

".pval.1dof" <- function(x,greater){
  almost.1 <- 1 + 64 * .Machine$double.eps
  
  jj <- allboards(x)
  or <- apply(jj,3,odds.ratio)
  p <- allprobs(x)
  x.or <- odds.ratio(x) / almost.1

  if(greater){
    return(sum(p[or > x.or]))
  } else {
    return(sum(p[or < x.or]))
  }
}

"dof" <- function(x){(nrow(x)-1)*(ncol(x)-1)-sum(is.na(x))}

"odds.ratio" <- function(x){
  stopifnot(is.1dof(x))
  n <- nrow(x)
  ind <- cbind(1:n,c(2:n,1))
  return(prod(diag(x))/prod(x[ind]))
}

"maxlike" <- function(x){
  warning("not coded up in C")
  allboards(x)[,,which.max(allprobs(x))]
}

"is.1dof" <- function(x){
  n <- nrow(x)
  if(!is.matrix(x) | n != ncol(x)){
    return(FALSE)
  }
  ind <- cbind(1:n,c(2:n,1))
  if(all(!is.na(diag(x))) & all(!is.na(x[ind])) & sum(is.na(x))==n*(n-2)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

"as.pairwise" <- function(x){
  stopifnot(nrow(x)==ncol(x))
  n <- nrow(x)
  k <- n * (n - 1) / 2
  out <- matrix(NA, k, n)
  upper.indexes <- which( lower.tri( x ), arr.ind=TRUE )
  from.mat <- rbind( upper.indexes, upper.indexes[ , 2:1 ] )
  to.mat <-
    cbind(rep(1:nrow(upper.indexes),2), as.vector(upper.indexes[, 2:1]))
  out[ to.mat ] <- x[ from.mat ]
  colnames(out) <- colnames(x)
  return(out)
}
  
"randomprobs" <- function(x, B=2000, n=100, burnin=0, use.brob=FALSE, func=NULL){
  x <- as.board(x)@x
  out <- rep(0,B)
  if(use.brob){
    out <- as.brob(out)
  }

  default <- FALSE
  if(is.null(func)){
    func <- function(x){prob(x, give.log=TRUE, use.brob=use.brob)}
    default <- TRUE
  }
  
  old <- x
  out[1] <- func(x)
  
  if(out[1] == -Inf){
    if(use.brob){
      stop("This cannot happen unless the board has astronomically large entries")
    } else {
      stop("Board has probability of zero (well, less than .Machine$double.xmin).  Consider setting use.brob to TRUE")
    }
  }
  
  for(i in seq_len(B+burnin)[-1]){
    proposed <- candidate(old, n=n)
    num <- prob(proposed, give.log=TRUE, use.brob=use.brob)
    den <- prob(old     , give.log=TRUE, use.brob=use.brob)
    
    if ((num == -Inf) & (den == -Inf)) {  #zero probability
      stop("this cannot happen.")
    } 
    alpha <- min(as.numeric(exp(num-den)),1)  #num, den are logs
    if (runif(1) < alpha){
      if(default){
        out[i] <- num
      } else {
        out[i] <- func(proposed)
      }
      old <- proposed
    } else {
      if(default){
        out[i] <- den
      } else {
        out[i] <- func(old)
      }
    }
  }
  if(burnin>0){
    out <- out[-seq_len(burnin)]
  }
  return(out)
}

"randomboards" <- function(x, B=2000, n=100, burnin=0){
  x <- as.board(x)@x
  out <- array(0L,c(nrow(x),ncol(x),B+burnin))

  old <- x
  out[,,1] <- x
  
  for(i in seq_len(B+burnin)[-1]){
    proposed <- candidate(old, n=n)
    num <- prob(proposed, give.log=TRUE)
    den <- prob(old     , give.log=TRUE)
    
    if ((num == -Inf) & (den == -Inf)) {  #zero probability
      stop("this cannot happen.")
    } 
    alpha <- min(as.numeric(exp(num-den)),1)  #num, den are logs
    if (runif(1) < alpha){
      out[,,i] <- proposed
      old <- proposed
    } else {
      out[,,i] <- old
    }
  }
  
  if(burnin>0){
    out <- out[,,-seq_len(burnin)]
  }
  dimnames(out) <- dimnames(x)
  return(out)
}

"best" <- function(x, func=NULL, n=100, ...){

  if(is.null(func)){
    func <- function(x){-prob(x)}
  }
  
  dims <- dim(x)
  ind <- which(is.na(x) , arr.ind=TRUE)

  tovec <- function(x){
    x[ind] <- -1
    as.vector(x)
  }
  
  tomat <- function(x){
    dim(x) <- dims
    x[ind] <- NA
    x
  }
  
  out <- optim(tovec(x) , fn=function(x){func(tomat(x))} , gr=function(x){tovec(candidate(tomat(x), n=n))} , method="SANN" , ...)
  out$par <- tomat(out$par)
  rownames(out$par) <- rownames(x)
  colnames(out$par) <- colnames(x)
  out
}


"good"   <- function(x, method = "D", ...){
  jj <- marginals(x)
  N <- sum(x,na.rm=TRUE)
  B <- exp(
           sum(lchoose(jj$rs+ncol(x)-1,jj$rs))+
           sum(lchoose(jj$cs+nrow(x)-1,jj$cs))-
           lchoose(N+nrow(x)*ncol(x)-1,N)
           )
  
  if(any(is.na(x)) & !method=="A"){
    warning("Good's method is for matrices with no NA entries.  Answer supplied is indicative only (but should provide an upper bound)")
  }
  
  return(
         switch(method,
                A = no.of.boards(x, ...),
                B = B,
                C = 1.3*N^2*B/(nrow(x)*sum(jj$rs^2)),
                D = 1.3*N^4*B/(nrow(x)*ncol(x)*sum(outer(jj$rs^2,jj$cs^2))),
                "method must be one of A-D"
                )
         )
}

"candidate" <- function(x, n=100, give=FALSE){
  stopifnot(is.matrix(x))
  m <- marginals(x)
  cx <- .Cargs(x)
  x[is.na(x)] <- 0
  flash <- c("randpath", cx, list(ans=as.integer(as.vector(x))), n=as.integer(n), PACKAGE="aylmer")
  jj <- do.call(".C",flash)
  n <- jj$n
  if(give){
    return(n)
  }
  if(n==0){
    print(x)
    stop("no acceptable candidates found.  Consider increasing n")
  }
  out <- jj$ans
  dim(out) <- c(cx$nrow,cx$ncol)
  out[m$na] <- NA
  rownames(out) <- rownames(x)
  colnames(out) <- colnames(x)
  return(out)
}  
