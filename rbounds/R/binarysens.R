binarysens <- function(x, y=NULL, Gamma=6, GammaInc=1) {
  if (length(x)==1) {
    ctrl <- x
    trt <- y
  }
  else {
    y.c <- x$mdata$Y[x$mdata$Tr==0]
    y.t <- x$mdata$Y[x$mdata$Tr==1]
    table(y.t, y.c)
    y.tmp1 <- table(y.t, y.c)[2]
    y.tmp2 <- table(y.t, y.c)[3]
    
    if(y.tmp1 >= y.tmp2) {
      trt <- y.tmp1
      ctrl <- y.tmp2
    }
    else {
      trt <- y.tmp2
      ctrl <- y.tmp1
    }
  }
  gamma <- seq(1, Gamma, by=GammaInc)
  mx <- ctrl + trt
  up <- c()
  lo <- c()
  series <- seq(trt, mx, by=1)
  n.it <- length(gamma)
  for(i in 1:n.it) {
    p.plus <- gamma[i]/(1 + gamma[i])
    p.minus <- 1/(1 + gamma[i])
    up.tmp <- sum(dbinom(series, mx, prob=p.plus))
    lo.tmp <- sum(dbinom(series, mx, prob=p.minus))
    up <- c(up, up.tmp)
    lo <- c(lo, lo.tmp)
  }
  
  pval <- lo[1]
  bounds <- data.frame(gamma, round(lo, 5), round(up, 5))
  colnames(bounds) <- c("Gamma", "Lower bound", "Upper bound")
  
  msg <- "Rosenbaum Sensitivity Test \n"
  note <- "Note: Gamma is Odds of Differential Assignment To
 Treatment Due to Unobserved Factors \n"
  
  Obj <- list(Gamma = Gamma, GammaInc = GammaInc, pval = pval,
              msg = msg, bounds = bounds, note = note)
  class(Obj) <- c("rbounds", class(Obj))
  
  Obj
}

