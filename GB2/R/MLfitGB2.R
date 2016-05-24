#Empirical estimates

main.emp <- function(z, w=rep(1, length(z))){

  median.emp <- weightedMedian(z, w)
  mean.emp   <- sum(w*z)/sum(w)
  arpr.emp   <- arpr(z, w)$value
  rmpg.emp   <- rmpg(z, w)$value
  qsr.emp    <- qsr(z, w)$value
  gini.emp   <- (gini(z, w)$value)/100
  main    <- c(median.emp, mean.emp, arpr.emp, rmpg.emp, qsr.emp, gini.emp)
  names(main) <- c("median", "mean", "arpr", "rmpg", "qsr", "gini")
  return(main)
}

# Fit by maximum likelihood, sample of persons

mlfit.gb2 <- function(z, w=rep(1, length(z))){
	
  d     <- data.frame(inc=z,w=w)
  d     <- d[!is.na(d$inc),]
   
# Truncate at 0

  inc   <- d$inc[d$inc > 0]
  w     <- d$w[d$inc > 0]

# Full log-likelihood fit

  fitf    <- ml.gb2(inc, w)$opt1
  af      <- fitf$par[1]
  bf      <- fitf$par[2]
  pf      <- fitf$par[3]
  qf      <- fitf$par[4]
  flik    <- fitf$value
  indicf  <- main.gb2(0.6, af, bf, pf, qf)

# Profile log-likelihood fit

  fitp    <- profml.gb2(inc, w)$opt1
  ap      <- fitp$par[1]
  bp      <- fitp$par[2]
  pp      <- prof.gb2(inc, ap, bp, w)[3]
  qp      <- prof.gb2(inc, ap, bp, w)[4]
  plik    <- fitp$value
  indicp  <- main.gb2(0.6, ap, bp, pp, qp)

# Values of the empirical estimates 
  
  indicE  <- main.emp(inc, w)
  
  type=c("Emp. est","ML full","ML prof")
  results <- data.frame(type=type,
        median=round(c(indicE[1],indicf[1],indicp[1])),
        mean=round(c(indicE[2],indicf[2],indicp[2])),
        ARPR=round(c(indicE[3],indicf[3],indicp[3]), digits=2),
        RMPG=round(c(indicE[4],indicf[4],indicp[4]), digits=2),
        QSR=round(c(indicE[5],indicf[5],indicp[5]), digits=2),
        GINI=round(c(indicE[6],indicf[6],indicp[6]), digits=2),
        likelihood=round(c(NA,flik,plik), digits=3),
        a=round(c(NA,af,ap), digits=2), b=round(c(NA,bf,bp), digits=2),p=round(c(NA,pf,pp), digits=2), q=round(c(NA,qf,qp), digits=2))

return(list(data.frame(results), fitf, fitp))
}

