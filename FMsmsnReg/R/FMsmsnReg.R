FMsmsnReg <- function(y, x1, Abetas = NULL, medj= NULL, sigma2 = NULL, shape = NULL, pii = NULL, g = NULL, get.init = TRUE, criteria = TRUE, group = FALSE,
                      family = "Skew.normal", error = 0.00001, iter.max = 100, obs.prob= FALSE, kmeans.param = NULL, show.convergence="TRUE", cp=0.4)
{
  if(ncol(as.matrix(y)) > 1) stop("This function is only for univariate response y!")
  if((family != "t") && (family != "Skew.t") && (family != "Skew.cn") && (family != "Skew.slash") && (family != "Skew.normal") && (family != "Normal")) {stop(paste("Family",family,"not recognized.\n",sep=" "))}
  if((length(g) == 0) && ((length(medj)==0) || (length(sigma2)==0) || (length(shape)==0) || (length(pii)==0)))  stop("The model is not specified correctly.\n")
  if(get.init == FALSE)
  {
   g <- length(medj)
   if((length(medj) != length(sigma2)) || (length(medj) != length(pii))) {stop("The size of the initial values are not compatibles.\n")}
   if((family == "Skew.t" || family == "Skew.cn" || family == "Skew.slash" || family == "Skew.normal") & (length(medj) != length(shape))) { stop("The size of the initial values are not compatibles.\n") }
   if(sum(pii) != 1) {stop("probability of pii does not sum to 1.\n")}
  }

  if((length(g)!= 0) && (g < 1)) { stop("g must be greater than 0.\n") }

  #Running the algorithm
  out <- EMsmsn.mixR(y, x1, Abetas, medj, sigma2, shape, pii, g, get.init, criteria, group, family, error, iter.max, obs.prob, kmeans.param)
  cat('\n')
  cat('------------------------------------------------------------\n')
  cat('Finite Mixture of Scale Mixture Skew Normal Regression Model\n')
  cat('------------------------------------------------------------\n')
  cat('\n')
  cat('Observations =',length(y))
  cat('\n')
  cat('\n')
  cat('Family =',family)
  cat('\n')
  cat('-----------\n')
  cat('Estimates\n')
  cat('-----------\n')
  cat('\n')
  print(round(out$ttable,5))
  cat('\n')
  if(criteria==TRUE)
  {
   cat('------------------------\n')
   cat('Model selection criteria\n')
   cat('------------------------\n')
   cat('\n')
   critFin <- c(out$loglik, out$aic, out$bic, out$edc, out$icl)
   critFin <- round(t(as.matrix(critFin)),digits=3)
   dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC","EDC","ICL"))
   print(critFin)
  cat('\n')
  }
  cat('-------\n')
  cat('Details\n')
  cat('-------\n')
  cat('\n')
  cat("Convergence reached? =",(out$iter < iter.max))
  cat('\n')
  cat('EM iterations =',out$iter,"/",iter.max)
  cat('\n')
  cat('Criteria =',round(out$criterio,9))
  cat('\n')
  cat("Processing time =",out$time,units(out$time))
  cat('\n\n')

  q <- ncol(x1)

  if(show.convergence=="TRUE" && (family == "Skew.t" || family == "Skew.slash"))
  {
    cpl    <- cp*iter.max
    npar   <- 4*g + 1 + q

    labels <- list()
    for(i in 1:q){labels[[i]]         = bquote(beta[.(i-1)])}
    for(i in 1:g){labels[[i + q]]     = bquote(sigma[.(i)]) }
    for(i in 1:g){labels[[i + q + 2]] = bquote(mu[.(i)])    }
    for(i in 1:g){labels[[i + q + 4]] = bquote(lambda[.(i)])}
    for(i in 1:g){labels[[i + q + 6]] = bquote(p[.(i)])}
    labels[[q+8+1]] = bquote(nu)

    par(mar=c(4, 4.5, 1, 0.5))
    op <- par(mfrow=c(ifelse(npar%%3==0,npar%/%3,(npar%/%3)+1),3))

    for(i in 1:npar)
    {
      plot.ts(out$tetam[i,],xlab="Iteration",ylab=labels[[i]])
      abline(v=cpl,lty=2)
    }
  }

  if(show.convergence=="TRUE" && family == "Skew.cn")
  {
    cpl    <- cp*iter.max
    npar   <- 4*g + 2 + q

    labels <- list()
    for(i in 1:q){labels[[i]]         = bquote(beta[.(i-1)])}
    for(i in 1:g){labels[[i + q]]     = bquote(sigma[.(i)]^2) }
    for(i in 1:g){labels[[i + q + 2]] = bquote(mu[.(i)])    }
    for(i in 1:g){labels[[i + q + 4]] = bquote(lambda[.(i)])}
    for(i in 1:g){labels[[i + q + 6]] = bquote(p[.(i)])}
    labels[[q+8+1]] = bquote(nu)
    labels[[q+8+2]] = bquote(gamma)

    par(mar=c(4, 4.5, 1, 0.5))
    op <- par(mfrow=c(ifelse(npar%%3==0,npar%/%3,(npar%/%3)+1),3))

    for(i in 1:npar)
    {
      plot.ts(out$tetam[i,],xlab="Iteration",ylab=labels[[i]])
      abline(v=cpl,lty=2)
    }
  }

  if(show.convergence=="TRUE" && family == "Skew.normal")
  {
    cpl    <- cp*iter.max
    npar   <- 4*g + q

    labels <- list()
    for(i in 1:q){labels[[i]]         = bquote(beta[.(i-1)])}
    for(i in 1:g){labels[[i + q]]     = bquote(sigma[.(i)]) }
    for(i in 1:g){labels[[i + q + 2]] = bquote(mu[.(i)])    }
    for(i in 1:g){labels[[i + q + 4]] = bquote(lambda[.(i)])}
    for(i in 1:g){labels[[i + q + 6]] = bquote(p[.(i)])}

    par(mar=c(4, 4.5, 1, 0.5))
    op <- par(mfrow=c(ifelse(npar%%3==0,npar%/%3,(npar%/%3)+1),3))

    for(i in 1:npar)
    {
      plot.ts(out$tetam[i,],xlab="Iteration",ylab=labels[[i]])
      abline(v=cpl,lty=2)
    }
  }

  par(mfrow=c(1,1))
  par(mar= c(5, 4, 4, 2) + 0.1)

  if(family=="Skew.normal")
   obj.out <- list(tetam = out$tetam, parameters = out$ttable,Abetas = out$Abetas,  medj=out$medj, sigma2 = out$sigma2, shape = out$shape, pii = out$pii, aic=out$aic, bic=out$bic, edc=out$edc, icl=out$icl,convergence = out$iter < iter.max)
  else
   obj.out <- list(tetam = out$tetam, parameters = out$ttable,Abetas = out$Abetas,  medj=out$medj, sigma2 = out$sigma2, shape = out$shape, pii = out$pii, nu = out$nu,aic=out$aic, bic=out$bic, edc=out$edc, icl=out$icl,convergence = out$iter < iter.max)
  class(obj.out) <- family
  return(obj.out)
}

