# Author: P. Poncet

vieu <-
function(x,                       # sample (the data)
         bw = NULL,               # bandwidth
         kernel = "gaussian",     # kernel used
         abc = FALSE,             # if FALSE, 'optim' is used
         ...)
{
###################################################
# Mode estimator based on the kernel estimate 
#  of the derivative of the density 
###################################################

  if (pmatch(tolower(kernel), "normal", nomatch = 0)) {
    kernel <- "gaussian"
  } else kernel <- match.arg(tolower(kernel), c(.kernelList, "uniform")) # '.kernelList' is defined in 'K.R'
    
  ## Initialization
  nx <- length(x)
  if (is.null(bw)) bw <- bw.SJ(x)
  
  #fn <-
  #function(z)
  #{
  #  z <- (rep(z,nx) - x)/rep(bw,nx)   #! enlever 'rep' : on ne peut pas appliquer 'fn' à un vecteur 'z' de longueur > 1 !
  #                                    #! cela permettrait d'accelerer la methode 'abc'  
  #  z <- do.call(paste(".kernel.d", kernel, sep = ""), list(z))$k
  #  return(sum(z))   #! à modifier
  #}

  fn <- 
  function(z)
  {
    mat <- z/bw - x/bw
    k <- do.call(paste(".kernel.d", kernel, sep = ""), list(mat))$k
    return(sum(k))
  }
  
  #FN <- 
  #function(z)
  #{
  #  mat <- kronecker(z/bw, t(-x/bw), FUN = "+")
  #  k <- do.call(paste(".kernel.d", kernel, sep = ""), list(mat))$k
  #  return(rowSums(k))
  #}
  
  if (!abc) {
    r <- uniroot(f = fn, interval = c(min(x), max(x)), ...)
    M <- r$root
  } else {
    FN <- Vectorize(fn)
    f <- abs(FN(x)) 
    M <- x[f == min(f)]
  }
  
  ## Output
  return(M)   
}
