MgFormule632 <- function(App,BCV,NoInf,SmallerBetter=FALSE){
  ## ---------------------------------------------------------
  ## for a non-overfitting model the diff performance
  ##
  ##                 Apparent - NoInf
  ## 
  ## should be large relative to the diff performance
  ## 
  ##                 Apparent - BCV 
  ## ---------------------------------------------------------
  if (SmallerBetter==TRUE){
    overfit <- (BCV - App) / (NoInf - App)
    overfit[BCV>=NoInf] <- 1
    # the random forest can have NoInf <= App 
    # overfit[BCV<=App|NoInf<=App] <- 0
    overfit[BCV<=App] <- 0
    w <- .632 / (1 - .368 * overfit)
    w[BCV<=App] <- 0
    B632Plus <- (1-w) * App + w * BCV
  }
  else{
    overfit <- (App - BCV) / (App - NoInf)
    overfit[NoInf>=BCV] <- 1
    # the random forest can have NoInf <= App
    # overfit[App<=BCV|App<=NoInf] <- 0
    overfit[App<=BCV] <- 0
    w <- .632 / (1 - .368 * overfit)
    w[App<=BCV] <- 0
    #    B632 <- (1-.632) * App + .632 * BCV
    B632Plus <- (1-w) * App + w * BCV
  }
  list(B632Plus=B632Plus,
       overfit=overfit,
       weight=w)
}
