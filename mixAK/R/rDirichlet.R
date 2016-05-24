##
##  PURPOSE:   Dirichlet distribution
##             * random numbers generation
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   07/11/2007
##
##  FUNCTIONS:  rDirichlet
##
## ======================================================================

## *************************************************************
## rDirichlet
## *************************************************************
rDirichlet <- function(n, alpha=c(1, 1))
{
  thispackage <- "mixAK"

  if (any(alpha <= 0)) stop("All alpha's must be positive.")
  K <- length(alpha)

  SAMPLE <- .C("rDirichlet_R", x=double(K*n),
                               alpha=as.double(alpha),
                               K=as.integer(K),
                               npoints=as.integer(n),
               PACKAGE=thispackage)

  if (n == 1){
    names(SAMPLE$x) <- names(alpha)
  }else{
    SAMPLE$x <- matrix(SAMPLE$x, byrow=TRUE, ncol=K, nrow=n)
    colnames(SAMPLE$x) <- names(mean)
    rownames(SAMPLE$x) <- 1:n    
  }  

  return(SAMPLE$x)
}

