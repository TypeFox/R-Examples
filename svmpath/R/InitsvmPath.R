"InitsvmPath" <-
  function(Rmat, cvec, const)
  {
    n <- length(cvec)
    zsmall <- c(Rmat, cvec, const)
    lenz <- length(zsmall)
    sol <- .Fortran("qp",
                    xn = as.double(rep(0,n+1)),
                    n = as.integer(n),
                    zsmall = as.double(zsmall),
                    lenz = as.integer(lenz),
                    inform = as.integer(0),
                    PACKAGE="svmpath"
                    )
    if (sol$inform!=0) print("convergence warning in initialization\n")
    list(alpha=sol$xn[1:n], obj=sol$zsmall[lenz])
  }

