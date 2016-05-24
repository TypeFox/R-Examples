multistate <- function(Ri = Ri, Ci = Ci){
    ## #########################################################
    ## negative log likelihood of a cox type multistate model ##
    ## a.k.a. partial likelihood loss function                ##
    ## #########################################################
    plloss <- function(y, f, w) {
        event <- y[, 3]
        n <- length(event)
        if (length(f) == 1) {f <- rep(f, n)}
        ef <- exp(f)
        risk <- do.call(c, lapply(X = Ri, FUN = helpfunctionmultistate1, ef = ef))
        lpl <- sum(event * (f - log(risk)))
        return(lpl)}
    ## #################################
    ## construct new boosting family ##
    ## #################################
    Family(
    ## #########################################
    ## negative gradient of the loss function ##
    ## #########################################
    ngradient = function(y, f, w) {
      event <- y[, 3]
      n <- length(event)
      if(length(f) == 1){f <- rep(f, n)}
      ef <- exp(f)
      dummy <- secondpart <- rep(0, n)
      dummy <- do.call(c, lapply(X = Ri, FUN = helpfunctionmultistate1, ef = ef))
      dummy[which(dummy==0)] <- 1e-05
      secondpart <- do.call(c, lapply(X = Ci, FUN = helpfunctionmultistate2, dummy = dummy))
      gradients <- event - (ef*secondpart)
      return(gradients)},
    ## ################
    ## risk function ##
    ## ################
    risk = risk <- function(y, f, w) -sum(plloss(y, f, w)),
    ## ##################
    ## offset function ##
    ## ##################
    offset = function(y, w) 0,
    ## ################
    ## loss function ##
    ## ################
    #loss = loss <- function(y, f) {
    #    event <- y[, 3]
    #    n <- length(event)
    #    if (length(f) == 1) {f <- rep(f, n)}
    #    ef <- exp(f)
    #    risk <- rep(0, n)
    #    for (i in 1:n) {
    #        risk[i] <- sum(ef[riskset[[i]]])
    #    }
    #    lpl <- sum(event * (f - log(risk)))
    #    return(-mean(lpl))},
    name = "family for multistate models (2014-10-21).")}
