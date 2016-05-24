bwcontrol <- function(maxiter = 500, tol = 1e-05, prt = TRUE,
          posdiff = TRUE, converge = expression(diff<tol)){
    list(maxiter=maxiter, tol=tol, prt=prt, posdiff=posdiff,
         converge=converge)
}

