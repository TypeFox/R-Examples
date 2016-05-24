nca.ssd <- function(conc, time, n.tail=3, dose=0, method=c("z", "boott"), conf.level=0.95, nsample=1000, data){

  # function to define variance-covariance matrix
  sigma <- function(t, J, time, conc){

    # define objects
    sigma <- matrix(nrow=2*J, ncol=2*J)
    sigma[(1:(2*J)),(1:(2*J))] <- 0
    data <- cbind(conc, log(conc), time)
	
    # variances of observed and log-transformed data
    sx <- tapply(data[,1], data[,3], var)
    sy <- tapply(data[,2], data[,3], var)

    # variance-covariance matrix of observed and log-transformed data
    for(i in 1:J){      dintern <- subset(data, data[,3]==t[i])
      sigma[i+J,i] <- cov(dintern[,1], dintern[,2])      sigma[i,i+J] <- sigma[i+J,i] 
    }
    diag(sigma) <- c(sx, sy)
    sigma <- ifelse(is.na(sigma), 0, sigma)
    return(sigma)	
  }

  # function to calculate PK parameters work horse of SSD
  pkparms <- function(time, conc, d, k, w, t, J){	
		
    # define variance coveriance matrix
    M <- sigma(t=t, J=J, time=time, conc=conc)

    # number of animals per time point
    n <- tapply(conc, time, length)

    # means and variances of concentrations
    x <- conc
    x[x==0] <- 1E-256
    xq <- as.vector(tapply(x, time, mean))	
					
    # log transformed values	
    y <- log(x)
    yq <- tapply(y, time, mean, na.rm=TRUE) 
    sy <- tapply(y, time, var, na.rm=TRUE)
	
    # elimination rate constant
    i <- c((k+1):J) 				
    a <- c(rep(NA, k), t[i]-mean(t[i]))			
    u <- c(rep(NA, k), a[i]/sum(a[i]^2))	
    lambda <- sum(u[i]*yq[i])*(-1)
	
    # define index for vectorizations
    u <- u[!is.na(u)]
    i <- c(1:(J-1))

    # AUC from 0 to infinity #
    vec <- c(w[i], w[J] + 1/lambda, rep(0,k), xq[J]*u/lambda^2)
    est <- sum(w*xq) + xq[J]/lambda
    var <- (vec%*%M%*%vec)
    auc <- c(est, var)

    # AUMC from 0 TO infinity #
    vec <- c(t[i]*w[i], t[J]*w[J] + t[J]/lambda + 1/lambda^2, rep(0,k), ((t[J]*lambda+2)*xq[J]*u)/lambda^3)
    est <- sum(t*w*xq) + (xq[J]/lambda)*(t[J] + 1/lambda)
    var <- (vec%*%M%*%vec)
    aumc <- c(est, var)
		
    # MRT #
    c1 <- (-aumc[1] - w[J]*lambda*aumc[1]) / (lambda*auc[1]^2)
    c2 <- (lambda*t[J] + lambda^2*t[J]*w[J]+1) / (lambda^2*auc[1])
    c3 <- (2*auc[1] - lambda*aumc[1] + lambda*t[J]*auc[1]) / (lambda^3*auc[1]^2)
    vec <- c(w[i]*(t[i]/auc[1] - aumc[1]/auc[1]^2), c1 + c2, rep(0,k), c3*xq[J]*u)
    est <- aumc[1] / auc[1]
    var <- (vec%*%M%*%vec)
    mrt <- c(est, var)		

    # half-life #
    hl <- c(log(2)*mrt[1], log(2)^2*mrt[2])

    # clearance #
    vec <- c(-w[i]*d/auc[1]^2, (-d-d*lambda*w[J])/(lambda*auc[1]^2), rep(0,k), (-1)*(xq[J]*u*d)/(lambda^2*auc[1]^2))
    est <- d / auc[1]
    var <- (vec%*%M%*%vec)
    cls <- c(est, var)
	
    # Volume of distribution at steady state
    A <- 2*xq[J] + 2*lambda*xq[J]*t[J] + 2*lambda^2*sum(t*w*xq)
    B1 <- auc[1]*lambda - 2*xq[J] - 2*lambda*xq[J]*t[J] - 2*lambda*xq[J]*w[J] + auc[1]*lambda^2*t[J] 
    B2 <- auc[1]*lambda^3*t[J]*w[J] - 2*lambda^2*sum(t[i]*w[i]*xq[i]) - 4*lambda^2*xq[J]*t[J]*w[J] 
    B3 <- -2*w[J]*lambda^3*sum(t*w*xq)
    B <- B1 + B2 + B3
    C <- 2*auc[1]*lambda - 2*xq[J] - 2*lambda*xq[J]*t[J] + auc[1]*lambda^2*t[J] - 2*lambda^2*sum(t*w*xq)
    c1 <- (d*w[i]*t[i])/auc[1]^2 - (d*w[i]*A)/(auc[1]^3*lambda^2) 
    c2 <- (d*B) / (auc[1]^3 * lambda^3)
    c3 <- (d*xq[J]*u*C) / (lambda^4 * auc[1]^3)
    vec <- c(c1, c2, rep(0,k), c3)
    est <- d*aumc[1] / auc[1]^2
    var <- (vec%*%M%*%vec)
    vss <- c(est, var)

    # summarize results #
    res <- rbind(auc, aumc, mrt, hl, cls, vss)
    rownames(res) <- c('AUC to infinity', 'AUMC to infinity', 'Mean residence time', 'Half-life', 'Clearance', 'Volume of Distribution')
    return(res)
  }


  get.confint <- function(conc=conc, time=time, k=k, d=d, alpha=alpha, nsample=0){
	
    # use mean over all n per time points as n
    n <- mean(tapply(conc, time, length))
    w <- .weight(unique(time))
    t <- unique(time)
    J <- length(t)

    # get estimates and variances
    obsv.parms <- pkparms(time=time, conc=conc, d=d, k=k, w=w, t=t, J=J)
    # get asympotic confidence intervals at level 1-alpha
    z <- qnorm(1-alpha/2)
    obsv.stderr <- sqrt(obsv.parms[,2]/n)
    asymp.lower <- obsv.parms[,1] - obsv.stderr*z
    asymp.upper <- obsv.parms[,1] + obsv.stderr*z
    asymp <- data.frame(est=obsv.parms[,1], stderr=obsv.stderr, lower=asymp.lower, upper=asymp.upper,method=rep('z',6))
    res <- asymp
    if(nsample>0){
      boot.stat <- matrix(nrow=nsample, ncol=6)
      for(i in 1:nsample){
        boot.conc <- as.vector(unlist(tapply(conc, time, sample, replace=TRUE)))
        boot.parms <- pkparms(conc=boot.conc, time=time, d=d, k=k, w=w, t=t, J=J)
        boot.stat[i,] <- (boot.parms[,1]-obsv.parms[,1]) / sqrt(boot.parms[,2]/n)
      }

      t.lb <- apply(boot.stat, 2, quantile, probs=c(alpha/2), method=5, na.rm=TRUE)
      t.ub <- apply(boot.stat, 2, quantile, probs=c(1-alpha/2), method=5, na.rm=TRUE)

      boott.lower <- obsv.parms[,1] - t.ub*obsv.stderr
      boott.upper <- obsv.parms[,1] - t.lb*obsv.stderr	
   
      res <- data.frame(est=rep(asymp$est,each=2), stderr=rep(asymp$stderr,each=2), 
              lower= as.vector(rbind(asymp.lower,boott.lower)), upper= as.vector(rbind(asymp.upper,boott.upper)),
              method=rep(c('z','boott'),6))
    }
 
    return(res)		
  }

  if(!missing(data)){
    cnames <- colnames(data)
    if(!any(cnames=='conc')){stop("data does not contain a variable conc")}
    if(!any(cnames=='time')){stop("data does not contain a variable time")}
    conc <- data$conc
    time <- data$time
    if(any(cnames=='group')){
      group <- data$group
    }else{
      group <- NULL
    }
  }

  if (!is.vector(time) || !is.vector(conc)) {stop('Argument time and/or conc invalid')}
  method <- match.arg(method,several.ok=TRUE)
  if(!(any(method=='boott'))) nsample <- 0 ## prevent bootstrapping if not required
  if(any(method=='boott') & nsample ==0) stop('boott called with zero resamples')
  # check input data
  data <- data.frame(conc,time)
  data <- data[order(data$time),]
  data <- na.omit(data)
  data <- subset(data, data$conc >= 0)
  conc <- data$conc
  time <- data$time

  k <- length(unique(time)) - n.tail
  alpha <- 1-conf.level

  res <- get.confint(conc=conc, time=time, k=k, d=dose, alpha=alpha, nsample=nsample)
  auc.obs <- auc.ssd(conc=conc, time=time, method=levels(res$method),conf.level=conf.level, nsample=nsample)$CIs

  if(any(method=='boott')){
    res <- rbind(auc.obs[2:1,-5],res)
  }else{
    res <- rbind(auc.obs[1,-5],res)
  }

  r.names <- c('AUC to tlast', 'AUC to infinity', 'AUMC to infinity', 'Mean residence time', 'non-compartmental half-life', 'Clearance', 'Volume of distribution at steady state')
  rownames(res) <- paste(conf.level*100,'% CI for the ', rep(r.names,each=length(levels(res$method))), ' using a ', sort(levels(res$method),decreasing=TRUE),' distribution', sep='')

  if(!any(method=='z')){
    res <- res[seq(2,14,2),]
  }
  out <- NULL
  out$CIs <- res
  out$design<-"ssd"
  out$est <- matrix(split(res,res$method)[[1]][,1],ncol=1)
  rownames(out$est) <- r.names
  colnames(out$est) <- 'est'
  out$conf.level <- conf.level
  class(out)<-"PK"	
  out$conc <- conc
  out$time <- time
  out$group <- NULL
  out$dose <- dose
  return(out)
}

