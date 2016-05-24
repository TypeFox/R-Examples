ci.boot.lmf <-
function(x, clevel = 0.05)
{
  #Check class
  if(class(x) != "boot.lmf")
    stop("This is not bootstrap replicates from a 'lmf' model")
  #Keep the call
  ret <- list(call = match.call())
  #Keep the number of bootstrap replicates
  ret$nboot <- x$nboot
  #Keep the what argument from x$call
  ret$what <- x$call$what
  #Keep the level argument
  ret$clevel <- clevel
  #Keep the unique ages from x$uage
  ret$uage <- x$uage
  #Keep the number of age classes from x$nage
  ret$nage <- x$nage
  #Set up quantile function
  qfn <- function(a, clevel) {paste("(", paste(format(
    quantile(a, c(clevel/2, 1 - clevel/2)), digits = getOption("digits") - 3),
                                               collapse = ","), ")", sep = "")}
  #Calculate confidence intervals after checking which parameters
  #has been bootstrapped
  if(!is.null(x$lboot))
  {
    #lboot - CI for transition matrix
    lt <- array(Reduce(cbind, x$lboot),
                dim = c(dim(x$l)[1], dim(x$l)[2], x$nboot))
    ret$l <- apply(lt, 1:2, qfn, clevel = clevel)
    #luvboot - CI for lambda, u and v
    ret$luv <- apply(x$luvboot, 2, qfn, clevel = clevel)
  }
  if(!is.null(x$djboot))
  {
    #djboot - CI for sigma2.dj
    ret$sigma2.dj <- apply(x$djboot, 2, qfn, clevel = clevel)
    #dboot - CI for sigma2.d
    ret$sigma2.d <- qfn(x$dboot, clevel = clevel)
    #atboot
    #No CI possible
    #Atboot
    #No CI possible
    #Mboot - CI for temporal covariance matrix
    Mt <- array(Reduce(cbind, x$Mboot),
                dim = c(x$npar, x$npar, x$nboot))
    ret$M <- apply(Mt, 1:2, qfn, clevel = clevel)
    #Inherit names
    dimnames(ret$M) <- dimnames(x$Mboot[[1]])
    #aMboot - CI for temporal alpha estimates
    aMt <- do.call(rbind, x$aMboot)
    ret$aM <- apply(aMt, 2, qfn, clevel = clevel)
    #atCboot
    #No CI possible
    #eboot - CI for environmental variance (sigma2.e)
    ret$sigma2.e <- qfn(x$eboot, clevel = clevel)
    #Anfboot - CI for temporal covariance matrix
    #assuming M = 0 (no fluctuating selection)
    Anft <- array(Reduce(cbind, x$Anfboot),
                  dim = c(x$npar, x$npar, x$nboot))
    ret$Anf <- apply(Anft, 1:2, qfn, clevel = clevel)
    #Inherit names
    dimnames(ret$Anf) <- dimnames(x$Anfboot[[1]])
    #anfboot
    anft <- do.call(rbind, x$anfboot)
    ret$anf <- apply(anft, 2, qfn, clevel = clevel)
  }
  class(ret) <- "ci.boot.lmf"
  ret
}
