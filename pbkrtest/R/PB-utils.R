
##########################################################
###
### Likelihood ratio statistic
###
##########################################################

getLRT <- function(largeModel, smallModel){
  UseMethod("getLRT")
}

getLRT.merMod <-
    getLRT.mer <-
        function(largeModel, smallModel){
    ll.small <- logLik(smallModel, REML=FALSE)
    ll.large <- logLik(largeModel, REML=FALSE)
    tobs     <- 2*(ll.large-ll.small)
    df11     <- attr(ll.large, "df") - attr(ll.small, "df")
    p.X2     <- 1-pchisq(tobs, df11)
    c(tobs=tobs, df=df11, p.value=p.X2)
}


getLRT.lm <- function(largeModel, smallModel){
  ll.small <- logLik(smallModel)
  ll.large <- logLik(largeModel)
  tobs     <- 2*(ll.large-ll.small)
  df11     <- attr(ll.large, "df") - attr(ll.small, "df")
  p.X2     <- 1-pchisq(tobs, df11)
  c(tobs=tobs, df=df11, p.value=p.X2)
}







