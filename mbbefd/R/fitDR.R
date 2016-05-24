
fitDR <- function(x, dist, method="mle", start=NULL, ...)
{
  if(any(is.na(x)))
    x <- x[!is.na(x)]
  if(any(x < 0 | x > 1))
    stop("Values outside [0,1] are not supported in fitDR.")
  method <- match.arg(method, c("mle", "tlmme"))
  dist <- match.arg(dist, c("oiunif", "oistpareto", "oibeta", "oigbeta", "mbbefd", "MBBEFD"))
  
  if(dist == "mbbefd")
  {
    initparmbbefd <- list(list(a=Trans.m10(0), b=Trans.1Inf(0)), 
                          list(a=Trans.0Inf(0), b=Trans.01(0)))
    #print(initparmbbefd)
    prefit <- prefitDR.mle(x, "mbbefd")
    #domain : (a,b) in (-1, 0) x (1, +Inf)
    if(all(!is.na(prefit[[1]])))
    {
      initparmbbefd[[1]] <- prefit[[1]]
      if(initparmbbefd[[1]]["b"] == 1)
        initparmbbefd[[1]]["b"] <- 2
    }
    #domain : (a,b) in (0, +Inf) x (0, 1)
    if(all(!is.na(prefit[[2]])))
    {
      initparmbbefd[[2]] <- prefit[[2]]
      if(initparmbbefd[[2]]["b"] == 1)
        initparmbbefd[[2]]["b"] <- 1/2
    }
    #print(initparmbbefd)
    if(method == "mle")
    {
      
      #wrap gradient -LL to match the call by fitdist
      grLL <- function(x, fix.arg, obs, ddistnam) 
        -grLLfunc(obs=obs, theta=x, dist="mbbefd")
      
      #domain : (a,b) in (-1, 0) x (1, +Inf)
      alabama1 <- mledist(x, distr="mbbefd", start=initparmbbefd[[1]], 
                        custom.optim= constrOptim.nl, hin=constrmbbefd1, 
                        control.outer=list(trace= FALSE), gr=grLL)
      #domain : (a,b) in (0, +Inf) x (0, 1)
      alabama2 <- mledist(x, distr="mbbefd", start=initparmbbefd[[2]], 
                        custom.optim= constrOptim.nl, hin=constrmbbefd2, 
                        control.outer=list(trace= FALSE), gr=grLL)
      
      if(alabama1$convergence == 100 && alabama2$convergence == 100)
        f1 <- alabama1
      else if(alabama1$convergence == 100 && alabama2$convergence != 100) 
        f1 <- alabama2
      else if(alabama1$convergence != 100 && alabama2$convergence == 100) 
        f1 <- alabama1
      else
      {
        if(alabama1$loglik > alabama2$loglik)
          f1 <- alabama1
        else
          f1 <- alabama2 
        
        
#         #computes Hessian of -LL at estimate values
#         f1hess <- -heLLfunc(obs=x, theta=f1$estimate, dist="mbbefd")
#         
#         if(all(!is.na(f1hess)) && qr(f1hess)$rank == NCOL(f1hess)){
#           f1$vcov <- solve(f1hess)
#           f1$sd <- sqrt(diag(f1$vcov))
#           f1$cor <- cov2cor(f1$vcov)
#         }else{
#           f1$vcov <- NA
#           f1$sd <- NA
#           f1$cor <- NA                            
#         }
      }
      f1 <- fitDR.addcomp(x=x, theta=f1$estimate, dist="mbbefd", method="mle", f1$convergence)
      
    }else if(method == "tlmme")
    {
      DIFF2 <- function(par, obs) 
      {
        (mmbbefd(1, par[1], par[2]) - mean(obs))^2 + (tlmbbefd(par[1], par[2]) - etl(obs))^2
      }
      alabama1 <- constrOptim.nl(unlist(initparmbbefd[[1]]), fn=DIFF2, hin=constrmbbefd1, 
                                 hin.jac=constrmbbefd1jac, obs=x, control.outer=list(trace=FALSE))
      alabama2 <- constrOptim.nl(unlist(initparmbbefd[[2]]), fn=DIFF2, hin=constrmbbefd2, 
                                 hin.jac=constrmbbefd2jac, obs=x, control.outer=list(trace=FALSE))
      if(alabama1$convergence > 0 && alabama2$convergence > 0)
        f1 <- list(estimate=NA, convergence=100)
      else if(alabama1$convergence > 0 && alabama2$convergence == 0) 
        f1 <- list(estimate=alabama2$par, convergence=0)
      else if(alabama1$convergence == 0 && alabama2$convergence > 0) 
        f1 <- list(estimate=alabama1$par, convergence=0)
      else
      {
        if(alabama1$value < alabama2$value)
          f1 <- list(estimate=alabama1$par, convergence=0)
        else
          f1 <- list(estimate=alabama2$par, convergence=0)
      }
      f1 <- fitDR.addcomp(x=x, theta=f1$estimate, hessian=f1$hessian, dist="mbbefd", method="tlmme", f1$convergence)
      
    }else
      stop("not yet implemented")
  }else if(dist == "MBBEFD")
  {
    #starting values
    g <- 1/etl(x, na.rm=TRUE)
    if(is.infinite(g))
      g <- 2
    initparMBBEFD <- list(list(g=g, b=Trans.1Inf(0)), list(g=g, b=1/(2*g)))
    
    #try to improve initial value for b
    phalf <- mean(x <= 1/2)
    b <- Re(polyroot(c(phalf, (1-g)*(1-phalf), - 1 +g*(1-phalf))))
    if(any(b > 0 | b < 1))
      initparMBBEFD[[2]]["b"] <- max(b[b > 0 | b < 1])
    else if(any(b > 1))
      initparMBBEFD[[2]]["b"] <- max(b[b > 1])
      
#     cat("initial values")  
#     print(unlist(initparMBBEFD))
#     cat("constr > 0 : 1st set\n")
#     print(constrMBBEFD1(unlist(initparMBBEFD[[1]])))
#     cat("constr > 0 : 2nd set\n")
#     print(constrMBBEFD2(unlist(initparMBBEFD[[2]])))
#     
    
    prefit <- prefitDR.mle(x, "MBBEFD")
    if(all(!is.na(prefit[[1]])))
    { 
      initparMBBEFD[[1]] <- prefit[[1]]
      if(initparMBBEFD[[1]]["b"] == 1)
        initparMBBEFD[[1]]["b"] <- 2
    }
    if(all(!is.na(prefit[[2]])))
    {  
      initparMBBEFD[[2]] <- prefit[[2]]
      if(initparMBBEFD[[2]]["b"] == 1)
        initparMBBEFD[[2]]["b"] <- 1/2
    }
#     cat("after prefit\n")
#     print(unlist(initparMBBEFD))
#     
#     cat("constr > 0 : 1st set\n")
#     print(constrMBBEFD1(unlist(initparMBBEFD[[1]])))
#     cat("constr > 0 : 2nd set\n")
#     print(constrMBBEFD2(unlist(initparMBBEFD[[2]])))
    
    #cat("g", g, "b", b, "\n")
    if(method == "mle")
    {
      
      #domain : (g,b) in (1, +Inf) x (1, +Inf) with gb > 1
      alabama1 <- mledist(x, distr="MBBEFD", start=initparMBBEFD[[1]], 
                          custom.optim= constrOptim.nl, hin=constrMBBEFD1, 
                          control.outer=list(trace= FALSE), hin.jac=constrMBBEFD1jac, silent=TRUE)
      #domain : (g,b) in (1, +Inf) x (0, 1) with gb < 1
      alabama2 <- mledist(x, distr="MBBEFD", start=initparMBBEFD[[2]], 
                          custom.optim= constrOptim.nl, hin=constrMBBEFD2, 
                          control.outer=list(trace= FALSE), hin.jac=constrMBBEFD2jac, silent=TRUE)
      
      #print(summary(alabama1))
      #print(summary(alabama2))
      if(alabama1$convergence == 100 && alabama2$convergence == 100)
        f1 <- alabama1
      else if(alabama1$convergence == 100 && alabama2$convergence != 100) 
        f1 <- alabama2
      else if(alabama1$convergence != 100 && alabama2$convergence == 100) 
        f1 <- alabama1
      else
      {
        if(alabama1$loglik > alabama2$loglik)
          f1 <- alabama1
        else
          f1 <- alabama2 
        
#         if(all(!is.na(f1$hessian)) && qr(f1$hessian)$rank == NCOL(f1$hessian)){
#           f1$vcov <- solve(f1$hessian)
#           f1$sd <- sqrt(diag(varcovar))
#           f1$cor <- cov2cor(varcovar)
#         }else{
#           f1$vcov <- NA
#           f1$sd <- NA
#           f1$cor <- NA                            
#         }
#         
      }
      f1 <- fitDR.addcomp(x=x, theta=f1$estimate, hessian=f1$hessian, dist="MBBEFD", method="mle", f1$convergence)
      
      
    }else if(method == "tlmme")
    {
      DIFF2 <- function(par, obs) 
      {
        (mMBBEFD(1, par[1], par[2]) - mean(obs))^2 + (tlMBBEFD(par[1], par[2]) - etl(obs))^2
      }
      alabama1 <- constrOptim.nl(unlist(initparMBBEFD[[1]]), fn=DIFF2, hin=constrMBBEFD1, 
                                 obs=x, control.outer=list(trace=FALSE))
      alabama2 <- constrOptim.nl(unlist(initparMBBEFD[[2]]), fn=DIFF2, hin=constrMBBEFD2, 
                                 obs=x, control.outer=list(trace=FALSE))
      
      if(alabama1$convergence > 0 && alabama2$convergence > 0)
        f1 <- list(estimate=NA, convergence=100)
      else if(alabama1$convergence > 0 && alabama2$convergence == 0) 
        f1 <- list(estimate=alabama2$par, convergence=0)
      else if(alabama1$convergence == 0 && alabama2$convergence > 0) 
        f1 <- list(estimate=alabama1$par, convergence=0)
      else
      {
        if(alabama1$value < alabama2$value)
          f1 <- list(estimate=alabama1$par, convergence=0)
        else
          f1 <- list(estimate=alabama2$par, convergence=0)
      }
      f1 <- fitDR.addcomp(x=x, theta=f1$estimate, hessian=f1$hessian, dist="MBBEFD", method="tlmme", f1$convergence)
      
      
    }else
      stop("not yet implemented")  
      
  }else if(dist == "oiunif")
  {
    if(is.null(start))
      start <- list(p1=etl(x))
    
    #print(LLfunc(x, start$p1, dist))
    
    if(method %in% c("mle", "tlmme"))
    {
      if(method == "tlmme")
        method <- "mle"
      f1 <- fitdist(x, distr=dist, method=method, start=start,
                  lower=0, upper=1, ..., optim.method="Brent") #, control=list(trace=6, REPORT=1)
    
      
    }else
      stop("not yet implemented")  
    
  }else if(dist %in% c("oistpareto", "oibeta", "oigbeta")) #one-inflated distr
  {
    p1 <- etl(x, na.rm=TRUE)
    xneq1 <- x[x != 1]
    distneq1 <- substr(dist, 3, nchar(dist))
    
    uplolist <- list(upper=Inf, lower=0)
    if(is.null(start))
    {
      if(distneq1 == "stpareto")
      {
        start <- list(a=1)               
      }else if(distneq1 == "beta")
      {
        n <- length(xneq1)
        m <- mean(xneq1, na.rm=TRUE)
        v <- (n - 1)/n*var(xneq1, na.rm=TRUE)
        aux <- m*(1-m)/v - 1
        start <- list(shape1=m*aux, shape2=(1-m)*aux)
        
      }else if(distneq1 == "gbeta")
      {
        shape00 <- optimize(function(z) (Theil.emp(x, na.rm=TRUE) - Theil.theo.shape0(z, obs=x))^2, 
                            lower=0.01, upper=100)$minimum
        start <- c(list(shape0=shape00), as.list(fitdist(x^shape00, "beta", method="mme")$estimate))
      }else
        stop("wrong non-inflated distribution.")
    }else
    {
      if(distneq1 == "stpareto")
      {
        start <- start["a"]               
      }else if(distneq1 == "beta")
      {
        start <- start[c("shape1", "shape2")]
      }else if(distneq1 == "gbeta")
      {
        start <- start[c("shape0", "shape1", "shape2")]
      }else
        stop("wrong non-inflated distribution.")
    }
      
    if(method == "mle")
    {
      #check the initial value
      loglik0 <- LLfunc(xneq1, unlist(start), distneq1)
        
      if(is.infinite(loglik0))
        stop("initial value of the log-likelihood is infinite.")
    
        #improve initial parameters for GB1
      if(distneq1 == "gbeta")
      {
        
        prefit <- prefitDR.mle(x, "oigbeta")
        
        #f0 <- mledist(xneq1, dist="gbeta2", optim.method="BFGS", 
        #              control=list(trace=0, REPORT=1, maxit=100), start=lapply(start, log))
        #print(unlist(start))
        if(all(!is.na(prefit)))
          start <- as.list(prefit)
        #print(unlist(start))
        
        f1 <- fitdist(xneq1, distr=distneq1, method="mle", start=start, 
                      optim.method="Nelder-Mead", ...)
      }else
        f1 <- fitdist(xneq1, distr=distneq1, method="mle", start=start, 
                  lower=uplolist$lower, upper=uplolist$upper, ...)
      if(f1$convergence != 0)
      {
         stop("error in convergence when fitting data.")
      }else
      {
        f1$estimate <- c(f1$estimate, p1=p1) 
        f1$n <- length(x)
        f1$distname <- dist
        f1$data <- x
        
        #gof stat
        f1$loglik <- LLfunc(obs=x, theta=f1$estimate, dist=dist)
        npar <- length(f1$estimate)
        f1$aic <- -2*f1$loglik+2*npar
        f1$bic <- -2*f1$loglik+log(f1$n)*npar
        
        f1$vcov <- rbind(cbind(as.matrix(f1$vcov), rep(0, npar-1)), 
                         c(rep(0, npar-1), p1*(1-p1)))
        dimnames(f1$vcov) <- list(names(f1$estimate), names(f1$estimate))
        
        f1$sd <- sqrt(diag(f1$vcov))
        f1$cor <- cov2cor(f1$vcov)
        
      } 
      
    }else if(method == "tlmme")
    {
      start <- c(start, list(p1=p1))
      npar <- length(start)
      
      DIFF2 <- function(par, obs) 
      {
        PX1 <- do.call(paste0("tl", dist), as.list(par))
        EX <- do.call(paste0("m", dist), as.list(c(order=1, par)))
        if(npar <= 2)
          return( (EX - mean(obs))^2 + (PX1 - etl(obs))^2 )
        
        if(npar >= 3)
          EX2 <- do.call(paste0("m", dist), as.list(c(order=2, par)))
        if(npar >= 4)
          EX3 <- do.call(paste0("m", dist), as.list(c(order=3, par)))
        
        if(npar == 3)
          return( (EX - mean(obs))^2 + (EX2 - mean(obs^2))^2 + (PX1 - etl(obs))^2 )
        else if(npar == 4)
          return( (EX - mean(obs))^2 + (EX2 - mean(obs^2))^2 + (EX3 - mean(obs^3))^2 + (PX1 - etl(obs))^2 )
        else
          stop("not implemented")
      }
          
      res <- optim(par=unlist(start), fn=DIFF2, obs=x, method="L-BFGS-B", 
                   lower=uplolist$lower, upper=uplolist$upper)
      
      if(res$convergence > 0)
        f1 <- list(estimate=NA, convergence=100)
      else
      {
        f1 <- list(estimate=res$par, convergence=0)
      }
      f1$method <- "tlmme"
      f1$data <- x
      f1$n <- length(x)
      f1$distname <- dist
      f1$fix.arg <- f1$fix.arg.fun <- f1$dots <- f1$weights <- NULL
      f1$discrete <- FALSE
      #gof stat
      f1$loglik <- LLfunc(obs=x, theta=f1$estimate, dist=dist)
      
      f1$aic <- -2*f1$loglik+2*npar
      f1$bic <- -2*f1$loglik+log(f1$n)*npar
      
      f1$sd <- f1$vcov <- f1$cor <- NA
      
    }else
    {
      stop("not yet implemented.")
    }
    
  }else
    stop("Unknown distribution for destruction rate models.")
  
  #reorder components as fitdist object
  f1 <- f1[c("estimate", "method", "sd", "cor", "vcov", "loglik", "aic", "bic", "n", "data", 
             "distname", "fix.arg", "fix.arg.fun", "dots", "convergence", "discrete", "weights")]
  class(f1) <- c("DR", "fitdist")
  f1
}

