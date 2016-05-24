waldtest.asreml <- function(object, cc)
{
      call <- match.call()
      if(oldClass(object) != "asreml")
        stop("Requires an object of class asreml\n")
      if(is.null(object$Cfixed)) {
        warning("Requires C matrix from model object. Refitting test model with argument \"Cfixed = TRUE\"\n")
        object <- update(object, Cfixed = TRUE)
      }  
      temp <- object$Cfixed
      tau <- object$coefficients$fixed
      nc <- length(tau)
      vrb <- matrix(0, nc, nc)
      vrb[!lower.tri(vrb)] <- temp[1:((nc * (nc + 1))/2)]
      vrb <- vrb + t(vrb) - diag(diag(vrb))
      sigma2 <- object$sigma2
      vrb <- vrb/sigma2
      res <- cintern(cc = cc, tau = tau, vrb = vrb, sigma2 = sigma2) 
      oldClass(res) <- "wald"
      res
}


