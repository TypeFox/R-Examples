CaoKosorok <- function(Tstar, 
                       gammas,
                       cSequence,
                       tSequence){

  if(is.null(cSequence)) cSequence <- seq(0.01, 6, 0.01) 
  if(is.null(tSequence)) tSequence <- seq(0.01, 6, 0.01)

  tst <- Tstar < 0
  Tstar[tst] <- -Tstar[tst]

  CKfunc <- function(x, Tstar){
    g_c.Tstar <- pmin(Tstar,x)/x
    g_c.Tstar.hat <- mean(g_c.Tstar)

    E.g_c <- 2 * (1 - exp(-x*x/2))/(x * sqrt(2*pi)) + 
             2 * (1 - pnorm(x))

    if(E.g_c >= 1){
      return(NA)
    } else {
      return( {g_c.Tstar.hat - E.g_c}/{1 - E.g_c} )
    }

  }

  pi1 <- sapply(X = cSequence, 
                FUN = CKfunc, 
                Tstar = Tstar)

  pi1.hat <- max(pi1)

  ip <- which.max(pi1)

  if(ip == length(cSequence)){
    warning("Maximum is at upper bound of cSequence.")
  } else if(ip == 1){
    warning("Maximum is at lower bound of cSequence.")
  }

  fdr <- sapply(X = tSequence, 
                FUN = function(x,phat){
                        meanT <- mean(Tstar >= x)
                        if(meanT != 0){
                          return(2 * {1-phat} * {1-pnorm(x)} / meanT)
                        } else {
                          return(NA)
                        }
                      }, 
                phat = pi1.hat)

  indicator.ck <- sapply(X = gammas, 
                         FUN = function(x){
                                 tst <- which(fdr <= x + 1e-8)
                                 if(length(tst) > 0){
                                   return(Tstar >= tSequence[tst[1]])
                                 } else {
                                   return(array(0,length(Tstar)))
                                 }
                               })

  if(any(is.na(indicator.ck))){
    warning(paste("NA's returned by Cao-Kosork method.",
                  "Verify tSequence input.\n", sep=""))
  }

  if(is.matrix(indicator.ck)){
    colnames(indicator.ck) <- round(gammas,3)
  } else {
    indicator.ck <- matrix(indicator.ck,nrow=1)
    colnames(indicator.ck) <- round(gammas,3)
  }

  return(list('ind' = indicator.ck,
              'pi1' = pi1.hat))

}

