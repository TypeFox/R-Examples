#------several benchmarks for mip testing
library(editrules)
library(lpSolveAPI)

genRandomEM <- function(nvar=10, nedits=10){
  i <- 1
  while(TRUE){
    A <- matrix(as.integer(rnorm(nvar*nedits)), ncol=nvar)
    ops <- rep("<=", nedits)
    b <- as.integer(rnorm(nedits))
    E <- as.editmatrix(A, ops=ops, b=b)
    cat("\r",i,": Checking feasibility")
    try(if (isFeasible(E)) {
      x <- as.integer(rnorm(nvar))
      names(x) <- getVars(E)
      return(list(E=E, x=x))
    })
    cat("\r")
    i <- i + 1
  }
}

genOrderedVars <- function(nvar=10){
  I <- diag(nvar)
  A <- cbind(I[,-1], 0) - I
  as.editmatrix(A, ops=rep("<=", nvar))
}

testOrdered <- function(nvars=10, method="localizer"){
  E <- genOrderedVars(nvars)
  res <- NULL
  data <- as.data.frame(matrix(0, ncol=nvars, nrow=nvars, dimnames=list(errors=NULL, vars=getVars(E))))
  for (errors in seq_len(nvars)){
    data[errors,seq_len(errors)] <- rev(seq_len(errors))  # creates errors...
  }
  
  
  le <- localizeErrors(E, data, method=method)
  
  cbind(le$status, nvars=nvars, method=method)
}

f <- file("benchordered.txt", open="wt")
writeLines(paste(c("weight", "degeneracy", "user", "system", "elapsed", "maxDurationExceeded", 
             "nvars", "method"), collapse="\t"), f)

m <- seq_len(100)
for (i in m){
  cat("nvars = ",i,"....")
  res <- testOrdered(nvars=i, method="mip")
  cat("\nWriting results\n")
  write.table(res, f, col.names=FALSE, row.names=FALSE)
  flush(f)
}
for (i in m){
  cat("nvars = ",i,"....")
  res <- testOrdered(nvars=i, method="localizer")
  cat("\nWriting results\n")
  write.table(res, f, col.names=FALSE, row.names=FALSE)
  flush(f)
}

close(f)

tab <- read.table("benchordered.txt", header=T)
str(tab)

#E <- genRandom2(nvar=100, nedits=20, r=4)
#E

# Ex <- genRandomEM(nvar=100, nedits=10)
# errorLocalizer.mip(Ex$E,Ex$x, verbose=TRUE)



