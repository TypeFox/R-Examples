## CreatePACS object
setClass(
         Class = "Pacs",
         representation = representation(
           y          = "numeric",
           x          = "matrix",
           n          = "integer",
           p          = "integer",
           nItMax     = "integer",
           lambda     = "numeric",
           epsPACS    = "numeric",
           betaInput  = "numeric",
           betaOutput = "numeric",
           a0         = "numeric",
           K          = "integer"
           ),
         prototype = prototype(
           y          = numeric(),
           x          = matrix(),
           n          = integer(),
           p          = integer(),
           nItMax     = integer(),
           lambda     = numeric(),
           epsPACS    = numeric(),
           betaInput  = numeric(),
           betaOutput = numeric(),
           a0         = numeric(),
           K          = integer()
           )
         )
fitPacs  <- function(Y = rnorm(10),
                     X = matrix(rnorm(50), nrow = 10), 
                     lambda=0.5,
                     betaInput=rnorm(10),
                     epsPACS=1e-5,
                     nItMax=1000){
  if (lambda <= 0) {
    stop("[Clere:fitPacs] Non negative (or  = 0) values for tuning parameter lambda are not allowed!\n",call. = FALSE)
  }
  if (epsPACS <= 0) {
    stop("[Clere:fitPacs] Non negative (or  = 0) values for tolerance parameter epsPACS are not allowed!\n",call. = FALSE)
  }
  
  n <- nrow(X)
  p <- ncol(X)
  x <- scale(X)
  y <- Y-mean(Y)  
  
  pacsObj <- new("Pacs", y = y, x = x,
                 n = as.integer(n), p = as.integer(p), nItMax = as.integer(nItMax),
                 lambda = lambda, epsPACS = epsPACS, betaInput = betaInput)  
  .Call("pacs", pacsObj, PACKAGE = "clere")
  littleeps <- 1e-7
  nround    <- round(-log10(littleeps))
  rB        <- round(pacsObj@betaOutput,nround)
  if(sum(rB==0,na.rm=TRUE)>0){
    K <- length(unique(abs(rB[which(rB!=0)])))
  }else{
    K <- length(unique(abs(rB)))
  }
  pacsObj@a0 <- mean(Y-X%*%pacsObj@betaOutput)
  pacsObj@K  <- K
  return( pacsObj )
}


