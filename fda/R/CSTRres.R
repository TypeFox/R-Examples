CSTRsse <- function(par, datstruct, fitstruct, CSTRbasis, lambda){
#%  2007.07.06 by Spencer Graves
#
  estimate = fitstruct$estimate;
  if(is.null(estimate))
    stop("fitstruct has no 'estimate' component")
#  Compare estimate with par
  if(sum(estimate) != length(par)){
    cat("ERROR in CSTRfn:  sum(estimate) != length(par)\n")
    cat("par = ", par, "; \n")
    stop("fitstruct$estimate = ",
         paste(estimate, collapse=" "))
  }
  m <- 0
  {
#  1.1.  Estimate kref starting from par or use fitstruct?
    if(estimate[1]){
      m <- m+1
      kref <- par[m]
      fitstruct$kref <- kref
    }
    else
      kref <- fitstruct$kref
  }
  {
#  1.2.  Estimate EoverR starting from parvec or use fitstruct?
    if(estimate[2]){
      m <- m+1
      EoverR <- par[m]
      fitstruct$EoverR <- EoverR
    }
    else
      EoverR<- fitstruct$EoverR
  }
  {
#  1.3.  Estimate 'a' starting from parvec or use fitstruct?
    if(estimate[3]){
      m <- m+1
      a <- par[m]
      fitstruct$a <- a
    }
    else
      a <- fitstruct$a
  }
  {
#  1.4.  Estimate b starting from parvec or use fitstruct?
    if(estimate[4]){
      m <- m+1
      b <- par[m]
      fitstruct$b <- b
    }
    else
      b<- fitstruct$b
  }
#
  res <- CSTRres(kref=kref, EoverR=EoverR, a=a, b=b,
               datstruct=datstruct, fitstruct=fitstruct,
               CSTRbasis=CSTRbasis, lambda=lambda,
               gradwrd=FALSE)
  sse <- sum(res^2)
  if(is.na(sse)){
    parch <- paste(par, collapse=", ")
    warning('Missing values returned by CSTRres with par = ',
        parch, '; ', sum(is.na(res)), ' NAs out of ', length(res))
  }
  sse
}

CSTRres <- function(kref=NULL, EoverR=NULL, a=NULL, b=NULL,
               datstruct, fitstruct, CSTRbasis, lambda,
               gradwrd=FALSE){
#%  2007.06.02 by Spencer Graves
##
## 1.  Construct 'parvec'
##
#  parvec <- vector("list", 4)
#  names(parvec) <- c("kref", "EoverR", "a", "b")
#  parvec$kref  <- kref
#  parvec$Eover <- EoverR
#  parvec$a     <- a
#  parvec$b     <- b
  parvec <- list(kref=kref, EoverR=EoverR, a=a, b=b)
  pv     <- unlist(parvec)
##
## 2.  Call CSTRfn
##
  cstr. <- CSTRfn(parvec=pv, datstruct=datstruct, fitstruct=fitstruct,
                  CSTRbasis=CSTRbasis, lambda=lambda, gradwrd=gradwrd)
  Res <- as.vector(cstr.$res)
  if(length(d.r <- dim(Res))>1){
    ys <- dimnames(Res)[[2]]
    if(!is.null(ys)){
      resNames <- t(outer(ys, 1:d.r[1], paste, sep=""))
      Res      <- as.vector(Res)
      names(Res) <- resNames
    }
  }
  if(gradwrd)
    attr(Res, "gradient") <- cstr.$Dres
#
  Res
}
