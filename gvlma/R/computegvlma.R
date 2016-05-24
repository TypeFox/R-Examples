"computegvlma" <-
function(lmobj, alphalevel, v)
{
  y <- resid(lmobj) + fitted(lmobj)
  X <- model.matrix(lmobj)
  n <- length(y)
  nr <- nrow(X)
  nc <- ncol(X)
  sighat <- sqrt(sum(resid(lmobj)^2)/n)
  stdresids <- resid(lmobj)/sighat
  fitted <- fitted(lmobj)
  betahat <- coef(lmobj)
  w <- X[,-1]
  wbar <- t(w) %*% rep(1, n)/n
  wbarmat <- matrix(rep(wbar, n), nrow = n, byrow = TRUE)
  wmwbar <- w - wbarmat
  ## wmwbar <- sweep(w, MARGIN = 2, STATS = wbar, FUN = "-")    
  sigwhat <- t(wmwbar) %*% wmwbar/n
  sigwhatinv <- solve(sigwhat)
  fittedmmean <- fitted - mean(y)
  w2 <- t(betahat[-1]) %*% sigwhat %*% betahat[-1]
  gam <- t(fittedmmean^2) %*% wmwbar/n
  w3 <- (gam) %*% sigwhatinv %*% t(gam)
  w4 <- sum(fittedmmean^4)/n
  sig2q3 <- w4 - w2^2 - w3
  vmvbar <- v - mean(v)
  v2 <- sum(vmvbar^2)/n
  nsq <- sqrt(n)
  Q <- (1/nsq) * c(sum(stdresids^3),
                   sum(stdresids^4 - 3),
                   (fittedmmean^2) %*% stdresids,
                   sum(vmvbar * (stdresids^2 - 1)))
  stat1 <- Q[1]^2/6
  stat1pvalue <- 1 - pchisq(stat1, 1)
  stat2 <- Q[2]^2/24
  stat2pvalue <- 1 - pchisq(stat2, 1)
  stat3 <- Q[3]^2/(w4 - w2^2 - w3)
  stat3pvalue <- 1 - pchisq(stat3, 1)
  stat4 <- Q[4]^2/(2 * v2)
  stat4pvalue <- 1 - pchisq(stat4, 1)
  Gstat4 <- stat1 + stat2 + stat3 + stat4
  Gstat4pvalue <- 1 - pchisq(Gstat4, 4)
  decision1 <- 1
  if(stat1pvalue > alphalevel)
    decision1 <- 0
  decision2 <- 1
  if(stat2pvalue > alphalevel)
    decision2 <- 0
  decision3 <- 1
  if(stat3pvalue > alphalevel)
    decision3 <- 0
  decision4 <- 1
  if(stat4pvalue > alphalevel)
    decision4 <- 0
  Gdecision <- 1
  if(Gstat4pvalue > alphalevel)
    Gdecision <- 0
  z <- lmobj
  z$GlobalTest <-
    list(LevelOfSignificance = alphalevel,
         GlobalStat4 = list(Value = Gstat4,
           pvalue = Gstat4pvalue,
           Decision = Gdecision),
         DirectionalStat1 = list(Value = stat1, 
           pvalue = stat1pvalue,
           Decision = decision1), 
         DirectionalStat2 = list(Value = stat2,
           pvalue = stat2pvalue,
           Decision = decision2),
         DirectionalStat3 = list(Value = stat3,
           pvalue = stat3pvalue,
           Decision = decision3),
         DirectionalStat4 = list(Value = stat4, 
           pvalue = stat4pvalue,
           Decision = decision4),
         timeseq = v
         )
  class(z) <- c("gvlma", "lm")
  z
}

