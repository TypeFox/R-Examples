# This is file ../spam/tests/demo_timing.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     








# We construct a few large matrices and we compare how much faster (slower)
# we are compared to the full matrix analysis.
# Since all the calculation are also done with full matrices, we do not
# exagerate with the sizes.

options( echo=FALSE)
library( spam, warn.conflict=FALSE)

set.seed(14)


# In the test function, we do not print out the actual times
# We would get too many differences pointed out!
compare <- function(expr1,expr2,tag=NULL)
  {
    if( !is.null(tag)) cat( "Comparing: ", tag, fill=TRUE)
    invisible(data.frame(full=system.time( expr1, TRUE)[1:3],
                     sparse=system.time( expr2, TRUE)[1:3],
                     row.names=c("user","system","elapsed")))
  }

xn <- 10
xm <- 12

# first start with a full matrix.
fmat1 <- matrix(rnorm(xn*xm),xn,xm)
smat1 <- as.spam(fmat1)

compare(fmat2 <- t(fmat1), smat2 <- t(smat1), "Transpose")

compare(ffmat <- fmat1 %*% fmat2,
        ssmat <- smat1 %*% smat2, "multiplication")

compare( solve(ffmat),  solve(ssmat), "solving")


compare(rbind(fmat1,fmat1),rbind(smat1,smat1))
compare(cbind(fmat1,fmat1),cbind(smat1,smat1))





# now create a sparse matrix.
fmat1[fmat1<3] <- 0
smat1 <- as.spam(fmat1)



compare(fmat2 <- t(fmat1), smat2 <- t(smat1), "Transpose")

compare(ffmat <- fmat1 %*% fmat2,
        ssmat <- smat1 %*% smat2, "multiplication")

compare(ffmat <- ffmat + diag(xn),
        ssmat <- ssmat + diag.spam(xn), "add identity")

compare(ffmat <- 1:xn %d+% ffmat,
        ssmat <- 1:xn %d+% ssmat, "add identity quicker")

compare( solve(ffmat),  solve(ssmat), "solving")

summary(ssmat)


# compare a few cbind/rbinds

compare(rbind(fmat1,fmat1),rbind(smat1,smat1))
compare(cbind(fmat1,fmat1),cbind(smat1,smat1))


options( echo=TRUE)





# illustrate the new spam x matrix multiply:
if (F){
  n <- 1000
  
  A <- spam(0,n,n)
  A[cbind(1:(n-1),2:n)] <- -c(2,rep(4,n-3),2)
  A[cbind(1:(n-2),3:n)] <- rep(1,n-2)
  A <- A + t( A)
  diag(A) <- c(1,5,rep(6,n-4),5,1)


  B <- array(rnorm(n*n),c(n,n))

  system.time(C1 <- .spam.matmul.mat(A,B))
  system.time(C2 <- .spam.matmul(A,B))
  norm(C1-C2)


}
