# ################################################
# ------------------------------------------------
# lm.fit
# ------------------------------------------------
# ################################################

library(pbdDMAT , quiet=T)

init.grid()



genmat <- function(n, p)
{
  x <<- matrix(rnorm(n*p, sd=10), n, p)
  y <<- matrix(rnorm(n, sd=10), n)
  
  invisible(0)
}


check.mdl <- function(x, y, bldim=2)
{
  dx <- as.ddmatrix(x, bldim=bldim)
  dy <- as.ddmatrix(y, bldim=bldim)

#  comm.cat("\n", quiet=T)

  mdl1 <- lm.fit(x, y)
  mdl2 <- lm.fit(dx, dy)

  comm.cat("--------------------------\n", quiet=T)
  comm.cat("size:\t\t", quiet=T)
  comm.cat(paste(n, "x", p, "\n", sep=""), quiet=T)

  tests <- logical(5)
  tests[1] <- all.equal(mdl1$rank, mdl2$rank)

  out1 <- mdl1$coefficients
  names(out1) <- NULL
  out2 <- as.vector(mdl2$coefficients)
  tests[2] <- all.equal(out1, out2)

#  out1 <- mdl1$effects
#  names(out1) <- NULL
#  out2 <- as.vector(mdl2$effects)
#  tests[2] <- all.equal(out1, out2)

  out1 <- mdl1$residuals
  names(out1) <- NULL
  out2 <- as.vector(mdl2$residuals)
  tests[3] <- all.equal(out1, out2)

  out1 <- mdl1$fitted.values
  names(out1) <- NULL
  out2 <- as.vector(mdl2$fitted.values)
  tests[4] <- all.equal(out1, out2)

  out1 <- mdl1$df.residual
  names(out1) <- NULL
  out2 <- mdl2$df.residual
  tests[5] <- all.equal(out1, out2)

  if (is.logical(tests)) {
    if (all(tests)==T)
      comm.cat("All ok:\t\tTRUE\n", quiet=T)
    else 
      comm.print(tests, quiet=T)
  } else {
    comm.cat("rank:\t\t", quiet=T)
    comm.cat(tests[1], quiet=T)
    comm.cat("\n", quiet=T) 

    comm.cat("coefficients:\t", quiet=T)
    comm.cat(tests[2], quiet=T)
    comm.cat("\n", quiet=T) 

    comm.cat("residuals:\t", quiet=T)
    comm.cat(tests[3], quiet=T)
    comm.cat("\n", quiet=T) 

    comm.cat("fitted.values:\t", quiet=T)
    comm.cat(tests[4], quiet=T)
    comm.cat("\n", quiet=T) 

    comm.cat("df.residual:\t", quiet=T)
    comm.cat(tests[5], quiet=T)
    comm.cat("\n", quiet=T) 
  }
  
  invisible(0)
}




comm.set.seed(seed=1234, diff=F)

n <- 8
p <- 6

genmat(n, p)
check.mdl(x, y)

genmat(n, p)
x[,4] <- x[,2]
check.mdl(x, y)

genmat(n, p)
x[,1] <- x[,2]
check.mdl(x, y)

genmat(n, p)
x[,6] <- x[,2]
check.mdl(x, y)

genmat(n, p)
x[,1] <- x[,2] <- x[,3] <- x[,6]
check.mdl(x, y)

genmat(n, p)
x[,1] <- x[,3] <- x[,4] <- x[,5]
check.mdl(x, y)


####genmat(n, p)
####x[,1] <- x[,2] <- x[,3] <- x[,4] <- x[,5] <- x[,6]
####check.mdl(x, y)
####comm.print("ok to disagree", quiet=T)




n <- 6
p <- 6


genmat(n, p)
check.mdl(x, y)

genmat(n, p)
x[,4] <- x[,2]
check.mdl(x, y)

genmat(n, p)
x[,1] <- x[,2]
check.mdl(x, y)

genmat(n, p)
x[,6] <- x[,2]
check.mdl(x, y)

genmat(n, p)
x[,1] <- x[,2] <- x[,3] <- x[,6]
check.mdl(x, y)

genmat(n, p)
x[,1] <- x[,3] <- x[,4] <- x[,5]
check.mdl(x, y)

#genmat(n, p)
#x[,1] <- x[,2] <- x[,3] <- x[,4] <- x[,5] <- x[,6]
#check.mdl(x, y)
#comm.print("ok to disagree", quiet=T)





n <- 6
p <- 8


genmat(n, p)
check.mdl(x, y)

genmat(n, p)
x[,4] <- x[,2]
check.mdl(x, y)

genmat(n, p)
x[,1] <- x[,2]
check.mdl(x, y)

genmat(n, p)
x[,1] <- x[,2] <- x[,3]
check.mdl(x, y)

genmat(n, p)
x[,1] <- x[,2] <- x[,3] <- x[,4]
check.mdl(x, y)



finalize()
