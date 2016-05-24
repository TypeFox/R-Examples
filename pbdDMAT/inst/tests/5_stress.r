### Parameter stress tests
library(pbdDMAT, quietly=TRUE)
init.grid(quiet=TRUE)

bldims <- c(2, 4, 5, 7, 12, 13, 16, 24, 64)
ictxts <- 0:2

m <- 10
n <- 6
p <- 12


comm.set.seed(1234, diff=FALSE)

x <- matrix(rnorm(m*n), m, n)
y <- matrix(rnorm(n*p), n, p)

dx <- as.ddmatrix(x)
dy <- as.ddmatrix(y)


stress <- function(x, y, dx, dy, f, bldims, ctxts)
{
  printfun <- gsub(capture.output(f)[2], pattern=" +", replacement=" ")
  comm.cat(paste("TESTING: ", printfun), "\n", quiet=TRUE)
  
  truth <- f(x, y)
  
  anyfail <- FALSE
  
  for (bldim in bldims)
  {
    for (ctxt in ictxts)
    {
      .BLDIM <- bldim
      .ICTXT <- ctxt
      
      test <- as.matrix(f(dx, dy))
      success <- comm.all(all.equal(test, truth))
      
      if (!success)
      {
        anyfail <- TRUE
        
        comm.cat("FAILURE IN FUNCTION: ", quiet=TRUE)
        comm.cat(printfun, quiet=TRUE)
        comm.cat(paste("\n.BLDIM=", .BLDIM, "\n.ICTXT=", .ICTXT), quiet=TRUE)
        comm.cat("\n", quiet=TRUE)
      }
    }
  }
  
  if (!anyfail)
    comm.cat("All passed!\n\n", quiet=TRUE)
  else
    comm.cat("\n", quiet=TRUE)
  
  invisible()
}



### %*%
f <- function(x, y) x %*% y
stress(x, y, dx, dy, f, bldims, ctxts)


### cov()
f <- function(x, y) cov(x)
stress(x, y, dx, dy, f, bldims, ctxts)

### svd()
f <- function(x, y) La.svd(x)$vt
stress(x, y, dx, dy, f, bldims, ctxts)



finalize()
