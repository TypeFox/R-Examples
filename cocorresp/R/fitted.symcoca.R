"fitted.symcoca" <-
function(object, ...)
  {
    rsum1 <- object$rowsum[[1]]
    csum1 <- object$colsum[[1]]
    tot1 <- sum(rsum1)
    rsum2 <- object$rowsum[[2]]
    csum2 <- object$colsum[[2]]
    tot2 <- sum(rsum2)
    exp1 <- rsum1 %*% t(csum1) / tot1
    exp2 <- rsum2 %*% t(csum2) / tot2
    Yhat1 <- exp1 * (1 + (object$scores$site$Y %*% t(object$scores$species$Y)))
    Yhat2 <- exp2 * (1 + (object$scores$site$X %*% t(object$scores$species$X)))
    rownames(Yhat1) <- rownames(Yhat2) <- rownames(object$scores$site$Y)
    retval <- list(Yhat1 = Yhat1, Yhat2 = Yhat2,
                   nam.dat = c(object$nam.dat$namY, object$nam.dat$namX))
    names(retval$nam.dat) <- c("namY","namX")
    class(retval) <- "fitted.symcoca"
    retval
  }

