nca.batch <- function(conc, time, n.tail=3, dose=0, method=c("z", "boott"), conf.level=0.95, nsample=1000, data){

  # function to define variance-covariance matrix
  sigma <- function(tp, Inc, y){
    # define objects
    J <- length(tp)
    sig <- matrix(nrow=2*J, ncol=2*J)
    sig[(1:(2*J)),(1:(2*J))] <- 0
  
    #### data <- cbind(conc, log(conc), tp)

    ## number of observations for each combination of observations 
    R <- Inc %*% t(Inc)
    rvec <- rowSums(Inc)

    ## average per timepoint
    y.mean <- matrix(0,ncol=1,nrow=nrow(y))
    y.mean[,1] <- rowMeans(y, na.rm=TRUE)
    ly.mean <- matrix(0,ncol=1,nrow=nrow(y))
    ly.mean[,1] <- rowMeans(log(y), na.rm=TRUE)
    ## deviation from mean    

    y.dif <- y - kronecker(matrix(1,ncol=ncol(Inc),nrow=1),y.mean)
    ly.dif <- log(y) - kronecker(matrix(1,ncol=ncol(Inc),nrow=1),ly.mean)

    for (j1 in 1:(nrow(y))) {
      for(j2 in (1):nrow(y)) {
        r <- R[j1,j2]
        sig[j1,j2] <- sum(y.dif[j1,]*y.dif[j2,],na.rm=TRUE)/((r-1)+(1-r/rvec[j1])*(1-r/rvec[j2]))
        sig[j1,J+j2] <- sum(y.dif[j1,]*ly.dif[j2,],na.rm=TRUE)/((r-1)+(1-r/rvec[j1])*(1-r/rvec[j2]))
        sig[J+j1,j2] <- sum(ly.dif[j1,]*y.dif[j2,],na.rm=TRUE)/((r-1)+(1-r/rvec[j1])*(1-r/rvec[j2]))
        sig[J+j1,J+j2] <- sum(ly.dif[j1,]*ly.dif[j2,],na.rm=TRUE)/((r-1)+(1-r/rvec[j1])*(1-r/rvec[j2]))
      }
    }

    ## enforce Cauchy-Schwarz inequality
    for (j1 in 1:(2*nrow(y))) {
      for(j2 in (1):(2*nrow(y))) {
        if( !is.nan(sig[j1,j2]) & abs(sig[j1,j2]) > sqrt(sig[j1,j1]*sig[j2,j2])){
            sig[j1,j2] <- sign(sig[j1,j2])*sqrt(sig[j1,j1]*sig[j2,j2])
        }
      }
    }
    sig <- ifelse(is.na(sig), 0, sig)

    return(sig)
  }

  # function to calculate PK parameters work horse of SSD
  pkparms <- function(time, conc, d, k, w, t, J, sigma){	
		
    # define variance coveriance matrix
    M <- sigma

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

  get.confint <- function(conc=conc, time=time, k=k, d=d, alpha=alpha, nsample=0, sigma){
	
    # use mean over all n per time points as n
    n <- min(tapply(conc, time, length))

    w <- .weight(unique(time))
    t <- unique(time)
    J <- length(t)

    # get estimates and variances
    obsv.parms <- pkparms(time=time, conc=conc, d=d, k=k, w=w, t=t, J=J, sigma=sigma)
    # get asympotic confidence intervals at level 1-alpha
    z <- qnorm(1-alpha/2)
    # correct rounding error. Negative variance will be set almost zero
    obsv.parms[is.nan(sqrt(obsv.parms[,2])),2]<-1E-8
    obsv.stderr <- sqrt(obsv.parms[,2]/n)
    asymp.lower <- obsv.parms[,1] - obsv.stderr*z
    asymp.upper <- obsv.parms[,1] + obsv.stderr*z
    asymp <- data.frame(est=obsv.parms[,1], stderr=obsv.stderr, lower=asymp.lower, upper=asymp.upper,method=rep('z',6))
    res <- asymp
    if(nsample>0){
      boot.stat <- matrix(nrow=nsample, ncol=6)
      for(i in 1:nsample){
        boot.conc <- as.vector(unlist(tapply(conc, time, sample, replace=TRUE)))
        boot.parms <- pkparms(conc=boot.conc, time=time, d=d, k=k, w=w, t=t, J=J, sigma=sigma)
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
    if(!any(cnames=='id')){stop("data does not contain a variable id")}
    if(!any(cnames=='conc')){stop("data does not contain a variable conc")}
    if(!any(cnames=='time')){stop("data does not contain a variable time")}
    temp <- .formattobatch(data)
    conc <- time <- group <- NULL
    for(i in 1:length(temp)){
      conc[[i]] <- temp[[i]]$conc 
      time[[i]] <- temp[[i]]$time
      if(any(names(temp[[1]])=='group')){
        group[[i]] <- temp[[i]]$group
      }
    }
  }  

#  if (!is.vector(time) || !is.vector(conc)) {stop('Argument time and/or conc invalid')}
  method <- match.arg(method,several.ok=TRUE)
  if(!(any(method=='boott'))) nsample <- 0 ## prevent bootstrapping if not required
  if(any(method=='boott') & nsample ==0) stop('boott called with zero resamples')
  # check input data

  if(!missing(data)){
    cnames <- colnames(data)
    if(!any(cnames=='id')){stop("data does not contain a variable id")}
    if(!any(cnames=='conc')){stop("data does not contain a variable conc")}
    if(!any(cnames=='time')){stop("data does not contain a variable time")}
    temp <- .formattobatch(data)
    conc <- time <- group <- NULL
    for(i in 1:length(temp)){
      conc[[i]] <- temp[[i]]$conc 
      time[[i]] <- temp[[i]]$time
      if(any(names(temp[[1]])=='group')){
        group[[i]] <- temp[[i]]$group
      }
    }
  }  

  alpha <- 1-conf.level

###CHECK FIRST IF DATAMATRIX IS SUPPLIED
  times <- sort(unique(unlist(time)))
  y <- .batchtomatrix(conc,time)
  J <- nrow(y)
  n <- ncol(y)
  #incidence matrix
  Inc <-matrix(0,nrow=J,ncol=n)
  Inc[!is.na(y)] <- 1

  tvec <- rep(times,n)
  tvec <- tvec[!is.na(y)]
  yvec <- as.vector(y)
  yvec <- yvec[!is.na(yvec)]
  d<-cbind(tvec,yvec)
  d<-d[order(d[,1]),]
  tvec<-d[,1]
  yvec<-d[,2]

  sig <- sigma(tp=times,Inc=Inc,y)
  res <- get.confint(conc=yvec, time=tvec, k=J-n.tail, d=dose, alpha=alpha, nsample=nsample, sigma=sig)

  auc.obs <- auc.batch(conc=conc, time=time, method=levels(res$method),conf.level=conf.level, nsample=nsample)$CIs

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
  out$design<-"batch"
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

