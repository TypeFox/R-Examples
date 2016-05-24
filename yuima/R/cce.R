
# cce() Cumulative Covariance Estimator

#
# CumulativeCovarianceEstimator
#

# returns a matrix of var-cov estimates

########################################################
#         functions for the synchronization
########################################################

# refresh time
## data: a list of zoos

refresh_sampling <- function(data){
  
  d.size <- length(data)
  
  if(d.size>1){
    
    #ser.times <- vector(d.size, mode="list")
    #ser.times <- lapply(data,time)
    ser.times <- lapply(lapply(data,"time"),"as.numeric")
    ser.lengths <- sapply(data,"length")
    ser.samplings <- vector(d.size, mode="list")
    #refresh.times <- c()
    
    if(0){
      tmp <- sapply(ser.times,"[",1)    
      refresh.times[1] <- max(tmp)
      
      I <- rep(1,d.size)
      
      for(d in 1:d.size){
        #ser.samplings[[d]][1] <- max(which(ser.times[[d]]<=refresh.times[1]))
        while(!(ser.times[[d]][I[d]+1]>refresh.times[1])){
          I[d] <- I[d]+1
          if((I[d]+1)>ser.lengths[d]){
            break
          }
        }
        ser.samplings[[d]][1] <- I[d]
      }
      
      i <- 1
      
      #while(all(sapply(ser.samplings,"[",i)<ser.lengths)){
      while(all(I<ser.lengths)){
        
        J <- I
        tmp <- double(d.size)
        
        for(d in 1:d.size){
          #tmp[d] <- ser.times[[d]][min(which(ser.times[[d]]>refresh.times[i]))
          repeat{
            J[d] <- J[d] + 1
            tmp[d] <- ser.times[[d]][J[d]]
            if((!(J[d]<ser.lengths))||(tmp[d]>refresh.times[i])){
              break
            }
          }
        }
        
        i <- i+1
        
        refresh.times[i] <- max(tmp)
        
        for(d in 1:d.size){
          #ser.samplings[[d]][i] <- max(which(ser.times[[d]]<=refresh.times[i]))
          while(!(ser.times[[d]][I[d]+1]>refresh.times[i])){
            I[d] <- I[d]+1
            if((I[d]+1)>ser.lengths[d]){
              break
            }
          }
          ser.samplings[[d]][i] <- I[d]
        }
        
      }
    }
    
    MinL <- min(ser.lengths)
    
    refresh.times <- double(MinL)
    refresh.times[1] <- max(sapply(ser.times,"[",1))
    
    #idx <- matrix(.C("refreshsampling",
    #                 as.integer(d.size),
    #                 integer(d.size),
    #                 as.double(unlist(ser.times)),
    #                 as.double(refresh.times),
    #                 as.integer(append(ser.lengths,0,after=0)),
    #                 min(sapply(ser.times,FUN="tail",n=1)),
    #                 as.integer(MinL),
    #                 result=integer(d.size*MinL))$result,ncol=d.size)
    idx <- matrix(.C("refreshsampling",
                     as.integer(d.size),
                     integer(d.size),
                     as.double(unlist(ser.times)),
                     as.double(refresh.times),
                     as.integer(ser.lengths),
                     as.integer(diffinv(ser.lengths[-d.size],xi=0)),
                     min(sapply(ser.times,FUN="tail",n=1)),
                     as.integer(MinL),
                     result=integer(d.size*MinL))$result,ncol=d.size)
    
    result <- vector(d.size, mode="list")
    
    for(d in 1:d.size){
      #result[[d]] <- data[[d]][ser.samplings[[d]]]
      result[[d]] <- data[[d]][idx[,d]]
    }
    
    return(result)
    
  }else{
    return(data)
  }
  
}

# refresh sampling for PHY
## data: a list of zoos

refresh_sampling.PHY <- function(data){
  
  d.size <- length(data)
  
  if(d.size>1){
    
    ser.times <- lapply(lapply(data,"time"),"as.numeric")
    ser.lengths <- sapply(data,"length")    
    #refresh.times <- max(sapply(ser.times,"[",1))
    ser.samplings <- vector(d.size,mode="list")
    
    if(0){
      for(d in 1:d.size){
        ser.samplings[[d]][1] <- 1
      }
      
      I <- rep(1,d.size)
      i <- 1
      
      while(all(I<ser.lengths)){
        
        flag <- integer(d.size)
        
        for(d in 1:d.size){
          while(I[d]<ser.lengths[d]){
            I[d] <- I[d]+1
            if(ser.times[[d]][I[d]]>refresh.times[i]){
              flag[d] <- 1
              ser.samplings[[d]][i+1] <- I[d]
              break
            }
          }
        }
        
        if(any(flag<rep(1,d.size))){
          break
        }else{
          i <- i+1
          tmp <- double(d.size)
          for(d in 1:d.size){
            tmp[d] <- ser.times[[d]][ser.samplings[[d]][i]]
          }
          refresh.times <- append(refresh.times,max(tmp))
        }
        
      }
    }
    
    MinL <- min(ser.lengths)
    
    refresh.times <- double(MinL)
    refresh.times[1] <- max(sapply(ser.times,"[",1))
    
    #obj <- .C("refreshsamplingphy",
    #          as.integer(d.size),
    #          integer(d.size),
    #          as.double(unlist(ser.times)),
    #          rtimes=as.double(refresh.times),
    #          as.integer(append(ser.lengths,0,after=0)),
    #          min(sapply(ser.times,FUN="tail",n=1)),
    #          as.integer(MinL),
    #          Samplings=integer(d.size*(MinL+1)),
    #          rNum=integer(1))
    obj <- .C("refreshsamplingphy",
              as.integer(d.size),
              integer(d.size),
              as.double(unlist(ser.times)),
              rtimes=as.double(refresh.times),
              as.integer(ser.lengths),
              as.integer(diffinv(ser.lengths[-d.size],xi=0)),
              min(sapply(ser.times,FUN="tail",n=1)),
              as.integer(MinL),
              Samplings=integer(d.size*(MinL+1)),
              rNum=integer(1))
    
    refresh.times <- obj$rtimes[1:obj$rNum]
    idx <- matrix(obj$Samplings,ncol=d.size)
    
    result <- vector(d.size, mode="list")
    
    for(d in 1:d.size){
      #result[[d]] <- data[[d]][ser.samplings[[d]]]
      result[[d]] <- data[[d]][idx[,d]]
    }
    
    return(list(data=result,refresh.times=refresh.times))
    
  }else{
    return(list(data=data,refresh.times=as.numeric(time(data[[1]]))))
  }
  
}

# Bibinger's synchronization algorithm
## x: a zoo  y: a zoo

Bibsynchro <- function(x,y){
  
  xtime <- as.numeric(time(x))
  ytime <- as.numeric(time(y))
  
  xlength <- length(xtime)
  ylength <- length(ytime)
  N.max <- max(xlength,ylength)
  
  if(0){
    mu <- integer(N.max)
    w <- integer(N.max)
    q <- integer(N.max)
    r <- integer(N.max)
    
    if(xtime[1]<ytime[1]){
      I <- 1
      while(xtime[I]<ytime[1]){
        I <- I+1
        if(!(I<xlength)){
          break
        }
      }
      #mu0 <- min(which(ytime[1]<=xtime))
      mu0 <- I
      if(ytime[1]<xtime[mu0]){
        #q[1] <- mu0-1
        q[1] <- mu0
      }else{
        #q[1] <- mu0
        q[1] <- mu0+1
      }
      r[1] <- 2
    }else if(xtime[1]>ytime[1]){
      I <- 1
      while(xtime[I]<ytime[1]){
        I <- I+1
        if(!(I<xlength)){
          break
        }
      }
      #w0 <- min(which(xtime[1]<=ytime))
      w0 <- I
      q[1] <- 2
      if(xtime[1]<ytime[w0]){
        #r[1] <- w0-1
        r[1] <- w0
      }else{
        #r[1] <- w0
        r[1] <- w0+1
      }
    }else{
      q[1] <- 2
      r[1] <- 2
    }
    
    i <- 1
    
    repeat{
      #while((q[i]<xlength)&&(r[i]<ylength)){
      if(xtime[q[i]]<ytime[r[i]]){
        #tmp <- which(ytime[r[i]]<=xtime)
        #if(identical(all.equal(tmp,integer(0)),TRUE)){
        #  break
        #}
        #mu[i] <- min(tmp)
        if(xtime[xlength]<ytime[r[i]]){
          break
        }
        I <- q[i]
        while(xtime[I]<ytime[r[i]]){
          I <- I+1
        }
        mu[i] <- I
        w[i] <- r[i]
        if(ytime[r[i]]<xtime[mu[i]]){
          #q[i+1] <- mu[i]-1
          q[i+1] <- mu[i]
        }else{
          #q[i+1] <- mu[i]
          q[i+1] <- mu[i]+1
        }
        r[i+1] <- r[i]+1
      }else if(xtime[q[i]]>ytime[r[i]]){
        #tmp <- which(xtime[q[i]]<=ytime)
        #if(identical(all.equal(tmp,integer(0)),TRUE)){
        #  break
        #}
        if(xtime[q[i]]>ytime[ylength]){
          break
        }
        mu[i] <- q[i]
        #w[i] <- min(tmp)
        I <- r[i]
        while(xtime[q[i]]>ytime[I]){
          I <- I+1
        }
        w[i] <- I
        q[i+1] <- q[i]+1
        if(xtime[q[i]]<ytime[w[i]]){
          #r[i+1] <- w[i]-1
          r[i+1] <- w[i]
        }else{
          #r[i+1] <- w[i]
          r[i+1] <- w[i]+1
        }
      }else{
        mu[i] <- q[i]
        w[i] <- r[i]
        q[i+1] <- q[i]+1
        r[i+1] <- r[i]+1
      }
      
      i <- i+1
      
      if((q[i]>=xlength)||(r[i]>=ylength)){
        break
      }
      
    }
    
    mu[i] <- q[i]
    w[i] <- r[i]
  }
  
  sdata <- .C("bibsynchro",
              as.double(xtime),
              as.double(ytime),
              as.integer(xlength),
              as.integer(ylength),
              mu=integer(N.max),
              w=integer(N.max),
              q=integer(N.max),
              r=integer(N.max),
              Num=integer(1))
  
  Num <- sdata$Num
  
  #return(list(xg=as.numeric(x)[mu[1:i]],
  #            xl=as.numeric(x)[q[1:i]-1],
  #            ygamma=as.numeric(y)[w[1:i]],
  #            ylambda=as.numeric(y)[r[1:i]-1],
  #            num.data=i))
  return(list(xg=as.numeric(x)[sdata$mu[1:Num]+1],
              xl=as.numeric(x)[sdata$q[1:Num]],
              ygamma=as.numeric(y)[sdata$w[1:Num]+1],
              ylambda=as.numeric(y)[sdata$r[1:Num]],
              num.data=Num))
}


##############################################################
#         functions for tuning parameter
##############################################################

# Barndorff-Nielsen et al. (2009)

RV.sparse <- function(zdata,frequency=1200,utime){
  
  znum <- as.numeric(zdata)
  ztime <- as.numeric(time(zdata))*utime
  n.size <- length(zdata)
  end <- ztime[n.size]
  
  grid <- seq(ztime[1],end,by=frequency)
  n.sparse <- length(grid)
  
  #result <- double(frequency)
  #result <- 0
  #I <- rep(1,n.sparse)
  
  #for(t in 1:frequency){
  #  for(i in 1:n.sparse){
  #    while((ztime[I[i]+1]<=grid[i])&&(I[i]<n.size)){
  #      I[i] <- I[i]+1
  #    }
  #  }
  #  result[t] <- sum(diff(znum[I])^2)
    #result <- result+sum(diff(znum[I])^2)
  #  grid <- grid+rep(1,n.sparse)
  #  if(grid[n.sparse]>end){
  #    grid <- grid[-n.sparse]
  #    I <- I[-n.sparse]
  #    n.sparse <- n.sparse-1
  #  }
  #}
  
  K <- floor(end-grid[n.sparse]) + 1
  
  zmat <- matrix(.C("ctsubsampling",
                    as.double(znum),
                    as.double(ztime),
                    as.integer(frequency),
                    as.integer(n.sparse),
                    as.integer(n.size),
                    as.double(grid),
                    result=double(frequency*n.sparse))$result,
                 n.sparse,frequency)
  
  result <- double(frequency)
  result[1:K] <- colSums(diff(zmat[,1:K])^2)
  result[-(1:K)] <- colSums(diff(zmat[-n.sparse,-(1:K)])^2)
  
  return(mean(result))
  #return(result/frequency)
  #return(znum[I])
}

Omega_BNHLS <- function(zdata,sec=120,utime){
  
  #q <- ceiling(sec/mean(diff(as.numeric(time(zdata))*utime)))
  #q <- ceiling(sec*(length(zdata)-1)/utime)
  ztime <- as.numeric(time(zdata))
  q <- ceiling(sec*(length(zdata)-1)/(utime*(tail(ztime,n=1)-head(ztime,n=1))))
  obj <- diff(as.numeric(zdata),lag=q)
  n <- length(obj)
  
  #result <- 0
  result <- double(q)
  
  for(i in 1:q){
    tmp <- obj[seq(i,n,by=q)]
    n.tmp <- sum(tmp!=0)
    if(n.tmp>0){
      result[i] <- sum(tmp^2)/(2*n.tmp)
    }
  }
  
  #return(result/q)
  return(mean(result))
}

NSratio_BNHLS <- function(zdata,frequency=1200,sec=120,utime){
  
  IV <- RV.sparse(zdata,frequency,utime)
  Omega <- Omega_BNHLS(zdata,sec,utime)
  
  return(Omega/IV)
}

# Pre-averaging estimator

selectParam.pavg <- function(data,utime,frequency=1200,sec=120,
                             a.theta,b.theta,c.theta){
  
  if(missing(utime)) utime <- ifelse(is.numeric(time(data[[1]])),23400,1)
  
  NS <- sapply(data,FUN=NSratio_BNHLS,frequency=frequency,sec=sec,utime=utime)
  coef <- (b.theta+sqrt(b.theta^2+3*a.theta*c.theta))/a.theta
  
  return(sqrt(coef*NS))
}

selectParam.TSCV <- function(data,utime,frequency=1200,sec=120){
  
  if(missing(utime)) utime <- ifelse(is.numeric(time(data[[1]])),23400,1)
  
  NS <- sapply(data,FUN=NSratio_BNHLS,frequency=frequency,sec=sec,
               utime=utime)
  
  return((12*NS)^(2/3))
}

selectParam.GME <- function(data,utime,frequency=1200,sec=120){
  
  if(missing(utime)) utime <- ifelse(is.numeric(time(data[[1]])),23400,1)
  
  NS <- sapply(data,FUN=NSratio_BNHLS,frequency=frequency,sec=sec,
               utime=utime)
  
  a.theta <- 52/35
  #b.theta <- (12/5)*(NS^2+3*NS)
  b.theta <- (24/5)*(NS^2+2*NS)
  c.theta <- 48*NS^2
  
  return(sqrt((b.theta+sqrt(b.theta^2+12*a.theta*c.theta))/(2*a.theta)))
}

selectParam.RK <- function(data,utime,frequency=1200,sec=120,cstar=3.5134){
  
  if(missing(utime)) utime <- ifelse(is.numeric(time(data[[1]])),23400,1)
  
  NS <- sapply(data,FUN=NSratio_BNHLS,frequency=frequency,sec=sec,
               utime=utime)
  
  return(cstar*(NS^(2/5)))
}

if(0){
selectParam_BNHLS <- function(yuima,method,utime=1,frequency=1200,sec=120,kappa=7585/1161216,
                              kappabar=151/20160,kappatilde=1/24,cstar=3.51){
  data <- get.zoo.data(yuima)
  
  switch(method,
         "PHY"=selectParam.PHY(data,utime,frequency,sec,kappa,kappabar,kappatilde),
         "GME"=selectParam.GME(data,utime,frequency,sec),
         "RK"=selectParam.RK(data,utime,frequency,sec,cstar),
         "PTHY"=selectParam.PHY(data,utime,frequency,sec,kappa,kappabar,kappatilde))
  
}
}

###############################################################
#      functions for selecting the threshold functions
###############################################################

# local universal thresholding method

## Bipower variation
### x: a zoo lag: a positive integer

BPV <- function(x,lag=1){
  
  n <- length(x)-1
  dt <- diff(as.numeric(time(x)))
  obj <- abs(diff(as.numeric(x)))/sqrt(dt)
  
  result <- (pi/2)*(obj[1:(n-lag)]*dt[1:(n-lag)])%*%obj[(1+lag):n]
  
  return(result)
  
}

## local univaersal thresholding method
### data: a list of zoos  coef: a positive number

local_univ_threshold <- function(data,coef=5,eps=0.1){
  
  d.size <- length(data)
  
  result <- vector(d.size,mode="list") 
  
  for(d in 1:d.size){
    
    x <- data[[d]]
    #x <- as.numeric(data[[d]])
    n <- length(x)-1
    dt <- diff(as.numeric(time(x)))
    obj <- abs(diff(as.numeric(x)))/sqrt(dt)
    #dx <- diff(as.numeric(x))
    
    #xtime <- time(x)
    #xtime <- time(data[[d]])
    #if(is.numeric(xtime)){
    #  r <- max(diff(xtime))
    #}else{
    #  xtime <- as.numeric(xtime)
      #r <- max(diff(xtime))/(length(x)-1)
    #  r <- max(diff(xtime))/(tail(xtime,n=1)-head(xtime,n=1))
    #}
    #K <- ceiling(sqrt(1/r))
    #K <- max(ceiling(sqrt(1/r)),3)
    K <- max(ceiling(n^0.5),3)
    
    coef2 <- coef^2
    rho <- double(n+1)
    
    #tmp <- coef2*sum(dx^2)*n^(eps-1)
    #tmp <- double(n+1)
    #tmp[-(1:(K-1))] <- coef2*n^eps*rollmean(dx^2,k=K-1,align="left")
    #tmp[1:(K-1)] <- tmp[K]
    #dx[dx^2>tmp[-(n+1)]] <- 0
    
    if(K<n){
      #rho[1:K] <- coef2*(mad(diff(as.numeric(x)[1:(K+1)]))/0.6745)^2
      #rho[1:K] <- coef2*2*log(K)*(mad(diff(x[1:(K+1)]))/0.6745)^2
      #for(i in (K+1):n){
      #  rho[i] <- coef2*2*log(length(x))*BPV(x[(i-K):(i-1)])/(K-2)
      #}
      #rho[-(1:K)] <- coef2*2*log(n)^(1+eps)*
      # #rollapply(x[-n],width=K,FUN=BPV,align="left")/(K-2)
      rho[-(1:(K-1))] <- coef2*n^eps*(pi/2)*
        rollmean((obj*dt)[-n]*obj[-1],k=K-2,align="left")
      #rho[-(1:(K-1))] <- coef2*n^eps*rollmean(dx^2,k=K-1,align="left")
      rho[1:(K-1)] <- rho[K]
    }else{
      #rho <- coef2*(mad(diff(as.numeric(x)))/0.6745)^2
      #rho <- coef2*(mad(diff(x))/0.6745)^2
      #rho <- coef2*2*log(n)^(1+eps)*BPV(x)
      rho <- coef2*n^(eps-1)*BPV(x)
      #rho <- coef2*sum(dx^2)*n^(eps-1)
    }
    
    result[[d]] <- rho[-(n+1)]
    
  }
  
  return(result)
  
}

#################################################################
#       functions for calculating the estimators
#################################################################

# Hayashi-Yoshida estimator
## data: a list of zoos

HY <- function(data) {  
  
  n.series <- length(data)
  
  # allocate memory
  ser.X <- vector(n.series, mode="list")     # data in 'x'
  ser.times <- vector(n.series, mode="list") # time index in 'x'
  ser.diffX <- vector(n.series, mode="list") # difference of data
  
  for(i in 1:n.series){
    # set data and time index
    ser.X[[i]] <- as.numeric(data[[i]]) # we need to transform data into numeric to avoid problem with duplicated indexes below
    ser.times[[i]] <- as.numeric(time(data[[i]]))
    # set difference of the data 
    ser.diffX[[i]] <- diff( ser.X[[i]] )
  }
  
  # core part of cce
  
  cmat <- matrix(0, n.series, n.series)  # cov
  for(i in 1:n.series){
    for(j in i:n.series){ 
      #I <- rep(1,n.series)
      #Checking Starting Point
      #repeat{
      #while((I[i]<length(ser.times[[i]])) && (I[j]<length(ser.times[[j]]))){
      #  if(ser.times[[i]][I[i]] >= ser.times[[j]][I[j]+1]){
      #    I[j] <- I[j]+1   
      #  }else if(ser.times[[i]][I[i]+1] <= ser.times[[j]][I[j]]){
      #    I[i] <- I[i]+1   
      #  }else{
      #    break
      #  }
      #}
      
      
      #Main Component
      if(i!=j){
        #while((I[i]<length(ser.times[[i]])) && (I[j]<length(ser.times[[j]]))) {
        #  cmat[j,i] <- cmat[j,i] + (ser.diffX[[i]])[I[i]]*(ser.diffX[[j]])[I[j]]
        #  if(ser.times[[i]][I[i]+1]>ser.times[[j]][I[j]+1]){
        #    I[j] <- I[j] + 1
        #  }else if(ser.times[[i]][I[i]+1]<ser.times[[j]][I[j]+1]){
        #    I[i] <- I[i] +1
        #  }else{
        #    I[i] <- I[i]+1
        #    I[j] <- I[j]+1
        #  }
        #}
        cmat[j,i] <- .C("HayashiYoshida",as.integer(length(ser.times[[i]])),
                        as.integer(length(ser.times[[j]])),as.double(ser.times[[i]]),
                        as.double(ser.times[[j]]),as.double(ser.diffX[[i]]),
                        as.double(ser.diffX[[j]]),value=double(1))$value
      }else{
        cmat[i,j] <- sum(ser.diffX[[i]]^2)
      }
      cmat[i,j] <- cmat[j,i]  
    }
    
  }
  
  return(cmat)
}


#########################################################

# Modulated realized covariance (Cristensen et al.(2010))
## data: a list of zoos  theta: a positive number
## g: a real-valued function on [0,1] (a weight function)
## delta: a positive number

MRC <- function(data,theta,kn,g,delta=0,adj=TRUE){
  
  n.series <- length(data)
  
  #if(missing(theta)&&missing(kn))
  #  theta <- selectParam.pavg(data,utime=utime,a.theta=151/80640,
  #                            b.theta=1/96,c.theta=1/6)
  if(missing(theta)) theta <- 1
  
  cmat <- matrix(0, n.series, n.series)  # cov
  
  # synchronize the data and take the differences
  #diffX <- lapply(lapply(refresh_sampling(data),"as.numeric"),"diff")
  #diffX <- do.call("rbind",diffX) # transform to matrix
  diffX <- diff(do.call("cbind",lapply(refresh_sampling(data),"as.numeric")))
  
  #Num <- ncol(diffX)
  Num <- nrow(diffX)
  
  if(missing(kn)) kn <- floor(mean(theta)*Num^(1/2+delta))
  
  kn <- min(max(kn,2),Num+1)
  
  weight <- sapply((1:(kn-1))/kn,g)
  psi2.kn <- sum(weight^2)
  
  # pre-averaging
  #myfunc <- function(dx)rollapplyr(dx,width=kn-1,FUN="%*%",weight)
  #barX <- apply(diffX,1,FUN=myfunc)
  barX <- filter(diffX,filter=rev(weight),method="c",sides=1)[(kn-1):Num,]
  
  cmat <- (Num/(Num-kn+2))*t(barX)%*%barX/psi2.kn
  
  if(delta==0){
    psi1.kn <- kn^2*sum(diff(c(0,weight,0))^2)
    scale <- psi1.kn/(theta^2*psi2.kn*2*Num)
    #cmat <- cmat-scale*diffX%*%t(diffX)
    cmat <- cmat-scale*t(diffX)%*%diffX
    if(adj) cmat <- cmat/(1-scale)
  }
  
  return(cmat)
}

#########################################################

# Pre-averaged Hayashi-Yoshida estimator
## data: a list of zoos  theta: a positive number
## g: a real-valued function on [0,1] (a weight function)
## refreshing: a logical value (if TRUE we use the refreshed data)

PHY <- function(data,theta,kn,g,refreshing=TRUE,cwise=TRUE){
  
  n.series <- length(data)
  
  #if(missing(theta)&&missing(kn))
  #  theta <- selectParam.pavg(data,a.theta=7585/1161216,
  #                            b.theta=151/20160,c.theta=1/24)
  if(missing(theta)) theta <- 0.15
    
  cmat <- matrix(0, n.series, n.series)  # cov
  
  if(refreshing){
    
    if(cwise){
      
      if(missing(kn)){# refreshing,cwise,missing kn
        
        theta <- matrix(theta,n.series,n.series) 
        theta <- (theta+t(theta))/2
        
        for(i in 1:n.series){
          for(j in i:n.series){ 
            if(i!=j){
              
              resample <- refresh_sampling.PHY(list(data[[i]],data[[j]]))
              dat <- resample$data
              ser.numX <- length(resample$refresh.times)-1
              
              # set data and time index
              ser.X <- lapply(dat,"as.numeric")
              ser.times <- lapply(lapply(dat,"time"),"as.numeric")
              
              # set difference of the data
              ser.diffX <- lapply(ser.X,"diff")
              
              # pre-averaging
              kn <- min(max(ceiling(theta[i,j]*sqrt(ser.numX)),2),ser.numX)
              weight <- sapply((1:(kn-1))/kn,g)
              psi.kn <- sum(weight)
              
              #ser.barX <- list(rollapplyr(ser.diffX[[1]],width=kn-1,FUN="%*%",weight),
              #                 rollapplyr(ser.diffX[[2]],width=kn-1,FUN="%*%",weight))
              ser.barX <- list(filter(ser.diffX[[1]],rev(weight),method="c",
                                      sides=1)[(kn-1):length(ser.diffX[[1]])],
                               filter(ser.diffX[[2]],rev(weight),method="c",
                                      sides=1)[(kn-1):length(ser.diffX[[2]])])
              
              ser.num.barX <- sapply(ser.barX,"length")-1
              
              # core part of cce
              #start <- kn+1
              #end <- 1
              #for(k in 1:ser.num.barX[1]){
              #  while(!(ser.times[[1]][k]<ser.times[[2]][start])&&((start-kn)<ser.num.barX[2])){
              #    start <- start + 1
              #  }
              #  while((ser.times[[1]][k+kn]>ser.times[[2]][end+1])&&(end<ser.num.barX[2])){
              #    end <- end + 1
              #  }
              #  cmat[i,j] <- cmat[i,j] + ser.barX[[1]][k]*sum(ser.barX[[2]][(start-kn):end])
              #}
              cmat[i,j] <- .C("pHayashiYoshida",
                              as.integer(kn),
                              as.integer(ser.num.barX[1]),
                              as.integer(ser.num.barX[2]),
                              as.double(ser.times[[1]]),
                              as.double(ser.times[[2]]),
                              as.double(ser.barX[[1]]),
                              as.double(ser.barX[[2]]),
                              value=double(1))$value
              
              cmat[i,j] <- cmat[i,j]/(psi.kn^2)
              cmat[j,i] <- cmat[i,j]
              
            }else{
              
              diffX <- diff(as.numeric(data[[i]]))
              
              kn <- min(max(ceiling(theta[i,i]*sqrt(length(diffX))),2),
                        length(data[[i]]))
              weight <- sapply((1:(kn-1))/kn,g)
              psi.kn <- sum(weight)
              
              #barX <- rollapplyr(diffX,width=kn-1,FUN="%*%",weight)
              barX <- filter(diffX,rev(weight),method="c",
                             sides=1)[(kn-1):length(diffX)]
              tmp <- barX[-length(barX)]
              #cmat[i,j] <- tmp%*%rollapply(tmp,width=2*(kn-1)+1,FUN="sum",partial=TRUE)/(psi.kn)^2
              cmat[i,j] <- tmp%*%rollsum(c(double(kn-1),tmp,double(kn-1)),
                                         k=2*(kn-1)+1)/(psi.kn)^2
              
            }
          }
        }
      }else{# refreshing,cwise,not missing kn
        
        kn <- matrix(kn,n.series,n.series)
        
        for(i in 1:n.series){
          for(j in i:n.series){
            
            if(i!=j){
              
              resample <- refresh_sampling.PHY(list(data[[i]],data[[j]]))
              dat <- resample$data
              ser.numX <- length(resample$refresh.times)-1
              
              # set data and time index
              ser.X <- lapply(dat,"as.numeric")
              ser.times <- lapply(lapply(dat,"time"),"as.numeric")
              
              # set difference of the data
              ser.diffX <- lapply(ser.X,"diff")
              
              # pre-averaging
              knij <- min(max(kn[i,j],2),ser.numX+1)
              weight <- sapply((1:(knij-1))/knij,g)
              psi.kn <- sum(weight)
              
              #ser.barX <- list(rollapplyr(ser.diffX[[1]],width=knij-1,FUN="%*%",weight),
              #                 rollapplyr(ser.diffX[[2]],width=knij-1,FUN="%*%",weight))
              ser.barX <- list(filter(ser.diffX[[1]],rev(weight),method="c",
                                      sides=1)[(knij-1):length(ser.diffX[[1]])],
                               filter(ser.diffX[[2]],rev(weight),method="c",
                                      sides=1)[(knij-1):length(ser.diffX[[2]])])
              
              ser.num.barX <- sapply(ser.barX,"length")-1
              
              # core part of cce
              #start <- knij+1
              #end <- 1
              #for(k in 1:ser.num.barX[1]){
              #  while(!(ser.times[[1]][k]<ser.times[[2]][start])&&((start-knij)<ser.num.barX[2])){
              #    start <- start + 1
              #  }
              #  while((ser.times[[1]][k+knij]>ser.times[[2]][end+1])&&(end<ser.num.barX[2])){
              #    end <- end + 1
              #  }
              #  cmat[i,j] <- cmat[i,j] + ser.barX[[1]][k]*sum(ser.barX[[2]][(start-knij):end])
              #}
              cmat[i,j] <- .C("pHayashiYoshida",
                              as.integer(knij),
                              as.integer(ser.num.barX[1]),
                              as.integer(ser.num.barX[2]),
                              as.double(ser.times[[1]]),
                              as.double(ser.times[[2]]),
                              as.double(ser.barX[[1]]),
                              as.double(ser.barX[[2]]),
                              value=double(1))$value
              
              cmat[i,j] <- cmat[i,j]/(psi.kn^2)
              cmat[j,i] <- cmat[i,j]
              
            }else{
              
              diffX <- diff(as.numeric(data[[i]]))
              # pre-averaging
              kni <- min(max(kn[i,i],2),length(data[[i]]))
              weight <- sapply((1:(kni-1))/kni,g)
              psi.kn <- sum(weight)
              
              #barX <- rollapplyr(diffX,width=kni-1,FUN="%*%",weight)
              barX <- filter(diffX,rev(weight),method="c",
                             sides=1)[(kni-1):length(diffX)]
              tmp <- barX[-length(barX)]
              
              # core part of cce
              #cmat[i,j] <- tmp%*%rollapply(tmp,width=2*(kni-1)+1,FUN="sum",partial=TRUE)/(psi.kn)^2
              cmat[i,j] <- tmp%*%rollsum(c(double(kni-1),tmp,double(kni-1)),
                                         k=2*(kni-1)+1)/(psi.kn)^2
            }
          }
        }
      }
    }else{# refreshing,non-cwise
      
      # synchronize the data
      resample <- refresh_sampling.PHY(data)
      data <- resample$data
      ser.numX <- length(resample$refresh.times)-1
      
      # if missing kn, we choose it following Barndorff-Nielsen et al. (2011)
      if(missing(kn)){
        kn <- min(max(ceiling(mean(theta)*sqrt(ser.numX)),2),ser.numX+1)
      }
      kn <- kn[1]
      
      weight <- sapply((1:(kn-1))/kn,g)
      
      # allocate memory
      ser.X <- vector(n.series, mode="list")     # data in 'x'
      ser.times <- vector(n.series, mode="list") # time index in 'x'
      ser.diffX <- vector(n.series, mode="list") # difference of data
      ser.barX <- vector(n.series, mode="list")  # pre-averaged data
      ser.num.barX <- integer(n.series)         
      
      for(i in 1:n.series){
        # set data and time index
        ser.X[[i]] <- as.numeric(data[[i]]) # we need to transform data into numeric to avoid problem with duplicated indexes below
        ser.times[[i]] <- as.numeric(time(data[[i]]))
        
        # set difference of the data 
        ser.diffX[[i]] <- diff( ser.X[[i]] )
        
        # pre-averaging
        #ser.barX[[i]] <- rollapplyr(ser.diffX[[i]],width=kn-1,FUN="%*%",weight)
        ser.barX[[i]] <- filter(ser.diffX[[i]],rev(weight),method="c",
                                sides=1)[(kn-1):length(ser.diffX[[i]])]
        ser.num.barX[i] <- length(ser.barX[[i]])-1
      }
      
      # core part of cce
      
      for(i in 1:n.series){
        for(j in i:n.series){ 
          if(i!=j){
            #start <- kn+1
            #end <- 1
            #for(k in 1:ser.num.barX[i]){
            #  while(!(ser.times[[i]][k]<ser.times[[j]][start])&&((start-kn)<ser.num.barX[j])){
            #    start <- start + 1
            #  }
            #  while((ser.times[[i]][k+kn]>ser.times[[j]][end+1])&&(end<ser.num.barX[j])){
            #    end <- end + 1
            #  }
            #  cmat[i,j] <- cmat[i,j] + ser.barX[[i]][k]*sum(ser.barX[[j]][(start-kn):end])
            #}
            cmat[i,j] <- .C("pHayashiYoshida",
                            as.integer(kn),
                            as.integer(ser.num.barX[i]),
                            as.integer(ser.num.barX[j]),
                            as.double(ser.times[[i]]),
                            as.double(ser.times[[j]]),
                            as.double(ser.barX[[i]]),
                            as.double(ser.barX[[j]]),
                            value=double(1))$value
            cmat[j,i] <- cmat[i,j]
          }else{
            tmp <- ser.barX[[i]][1:ser.num.barX[i]]
            #cmat[i,j] <- tmp%*%rollapply(tmp,width=2*(kn-1)+1,FUN="sum",partial=TRUE)
            cmat[i,j] <- tmp%*%rollsum(c(double(kn-1),tmp,double(kn-1)),
                                       k=2*(kn-1)+1) 
          }
        }
      }
      
      psi.kn <- sum(weight)
      cmat <- cmat/(psi.kn^2)
      
    }
  }else{# non-refreshing
    
    if(cwise){# non-refreshing,cwise
      
      theta <- matrix(theta,n.series,n.series) 
      theta <- (theta+t(theta))/2
      
      if(missing(kn)){
        ntmp <- matrix(sapply(data,"length"),n.series,n.series)
        ntmp <- ntmp+t(ntmp)
        diag(ntmp) <- diag(ntmp)/2
        kn <- ceiling(theta*sqrt(ntmp))
        kn[kn<2] <- 2 # kn must be larger than 1
      }
      kn <- matrix(kn,n.series,n.series)
      
      for(i in 1:n.series){
        for(j in i:n.series){ 
          if(i!=j){
            
            dat <- list(data[[i]],data[[j]])
            
            # set data and time index
            ser.X <- lapply(dat,"as.numeric")
            ser.times <- lapply(lapply(dat,"time"),"as.numeric")
            # set difference of the data
            ser.diffX <- lapply(ser.X,"diff")
            
            # pre-averaging
            knij <- min(kn[i,j],sapply(ser.X,"length")) # kn must be less than the numbers of the observations
            weight <- sapply((1:(knij-1))/knij,g)
            #ser.barX <- list(rollapplyr(ser.diffX[[1]],width=knij-1,FUN="%*%",weight),
            #                 rollapplyr(ser.diffX[[2]],width=knij-1,FUN="%*%",weight))
            ser.barX <- list(filter(ser.diffX[[1]],rev(weight),method="c",
                                    sides=1)[(knij-1):length(ser.diffX[[1]])],
                             filter(ser.diffX[[2]],rev(weight),method="c",
                                    sides=1)[(knij-1):length(ser.diffX[[2]])])
            
            ser.num.barX <- sapply(ser.barX,"length")-1
            psi.kn <- sum(weight)
            
            # core part of cce
            start <- knij+1
            end <- 1
            #for(k in 1:ser.num.barX[1]){
            #  while(!(ser.times[[1]][k]<ser.times[[2]][start])&&((start-knij)<ser.num.barX[2])){
            #    start <- start + 1
            #  }
            #  while((ser.times[[1]][k+knij]>ser.times[[2]][end+1])&&(end<ser.num.barX[2])){
            #    end <- end + 1
            #  }
            #  cmat[i,j] <- cmat[i,j] + ser.barX[[1]][k]*sum(ser.barX[[2]][(start-knij):end])
            #}
            cmat[i,j] <- .C("pHayashiYoshida",
                            as.integer(knij),
                            as.integer(ser.num.barX[1]),
                            as.integer(ser.num.barX[2]),
                            as.double(ser.times[[1]]),
                            as.double(ser.times[[2]]),
                            as.double(ser.barX[[1]]),
                            as.double(ser.barX[[2]]),
                            value=double(1))$value
            
            cmat[i,j] <- cmat[i,j]/(psi.kn^2)
            cmat[j,i] <- cmat[i,j]
            
          }else{
            
            diffX <- diff(as.numeric(data[[i]]))
            
            # pre-averaging
            kni <- min(kn[i,i],length(data[[i]])) # kn must be less than the number of the observations
            weight <- sapply((1:(kni-1))/kni,g)
            psi.kn <- sum(weight)
            
            #barX <- rollapplyr(diffX,width=kni-1,FUN="%*%",weight)
            barX <- filter(diffX,rev(weight),method="c",
                           sides=1)[(kni-1):length(diffX)]
            tmp <- barX[-length(barX)]
            #cmat[i,j] <- tmp%*%rollapply(tmp,width=2*(kni-1)+1,FUN="sum",partial=TRUE)/(psi.kn)^2
            cmat[i,j] <- tmp%*%rollsum(c(double(kni-1),tmp,double(kni-1)),
                                       k=2*(kni-1)+1)/(psi.kn)^2
            
          }
        }
      }
    }else{# non-refreshing, non-cwise
      
      # allocate memory
      ser.X <- vector(n.series, mode="list")     # data in 'x'
      ser.times <- vector(n.series, mode="list") # time index in 'x'
      ser.diffX <- vector(n.series, mode="list") # difference of data
      
      ser.numX <- integer(n.series)
      ser.barX <- vector(n.series, mode="list")
      
      for(i in 1:n.series){
        # set data and time index
        ser.X[[i]] <- as.numeric(data[[i]]) # we need to transform data into numeric to avoid problem with duplicated indexes below
        ser.times[[i]] <- as.numeric(time(data[[i]]))
        
        # set difference of the data 
        ser.diffX[[i]] <- diff( ser.X[[i]] )    
        ser.numX[i] <- length(ser.diffX[[i]])
      }
      
      # pre-averaging  
      if(missing(kn)){
        kn <- min(max(ceiling(mean(theta)*sqrt(sum(ser.numX))),2),
                  ser.numX+1)
      }
      kn <- kn[1]
      
      weight <- sapply((1:(kn-1))/kn,g)
      psi.kn <- sum(weight)
      
      for(i in 1:n.series){
        #ser.barX[[i]] <- rollapplyr(ser.diffX[[i]],width=kn-1,FUN="%*%",weight)
        ser.barX[[i]] <- filter(ser.diffX[[i]],rev(weight),method="c",
                                sides=1)[(kn-1):length(ser.diffX[[i]])]
      }
      
      ser.num.barX <- sapply(ser.barX,"length")-1
      
      # core part of cce
      
      cmat <- matrix(0, n.series, n.series)  # cov
      for(i in 1:n.series){
        for(j in i:n.series){ 
          if(i!=j){
            #start <- kn+1
            #end <- 1
            #for(k in 1:ser.num.barX[i]){
            #  while(!(ser.times[[i]][k]<ser.times[[j]][start])&&((start-kn)<ser.num.barX[j])){
            #    start <- start + 1
            #  }
            #  while((ser.times[[i]][k+kn]>ser.times[[j]][end+1])&&(end<ser.num.barX[j])){
            #    end <- end + 1
            #  }
            #  cmat[i,j] <- cmat[i,j] + ser.barX[[i]][k]*sum(ser.barX[[j]][(start-kn):end])
            #}
            cmat[i,j] <- .C("pHayashiYoshida",
                            as.integer(kn),
                            as.integer(ser.num.barX[i]),
                            as.integer(ser.num.barX[j]),
                            as.double(ser.times[[i]]),
                            as.double(ser.times[[j]]),
                            as.double(ser.barX[[i]]),
                            as.double(ser.barX[[j]]),
                            value=double(1))$value
            cmat[j,i] <- cmat[i,j]
          }else{
            tmp <- ser.barX[[i]][1:ser.num.barX[i]]
            #cmat[i,j] <- tmp%*%rollapply(tmp,width=2*(kn-1)+1,FUN="sum",partial=TRUE)
            cmat[i,j] <- tmp%*%rollsum(c(double(kn-1),tmp,double(kn-1)),
                                       k=2*(kn-1)+1)
          }
        }
      }
      
      cmat <- cmat/(psi.kn^2)
      
    }
  }
  
  return(cmat)
}


################################################################

# previous tick Two Scales realized CoVariance
## data: a list of zoos  c.two: a postive number

TSCV <- function(data,K,c.two,J=1,adj=TRUE,utime){
    
  #X <- do.call("rbind",lapply(refresh_sampling(data),"as.numeric"))
  #N <- ncol(X)
  X <- do.call("cbind",lapply(refresh_sampling(data),"as.numeric"))
  N <- nrow(X)
  
  if(missing(K)){
    if(missing(c.two)) c.two <- selectParam.TSCV(data,utime=utime)
    K <- ceiling(mean(c.two)*N^(2/3))
  }
  
  scale <- (N-K+1)*J/(K*(N-J+1))
  
  #diffX1 <- apply(X,1,"diff",lag=K)
  #diffX2 <- apply(X,1,"diff",lag=J)
  diffX1 <- diff(X,lag=K)
  diffX2 <- diff(X,lag=J)
  
  cmat <- t(diffX1)%*%diffX1/K-scale*t(diffX2)%*%diffX2/J
  if(adj) cmat <- cmat/(1-scale)
  
  return(cmat)
}

################################################################

# Generalized multiscale estimator
## data: a list of zoos  c.multi: a postive number

GME <- function(data,c.multi,utime){
  
  d.size <- length(data)
  
  if(missing(c.multi)) c.multi <- selectParam.GME(data,utime=utime)
    
  c.multi <- matrix(c.multi,d.size,d.size) 
  c.multi <- (c.multi+t(c.multi))/2
  
  cmat <- matrix(0,d.size,d.size) 
  for(i in 1:d.size){
    for(j in i:d.size){ 
      
      if(i!=j){
        
        sdata <- Bibsynchro(data[[i]],data[[j]])
        N <- sdata$num.data
        
        M <- ceiling(c.multi[i,j]*sqrt(N))
        M <- min(c(M,N)) # M must be smaller than N
        M <- max(c(M,2)) # M must be lager than 2
        
        alpha <- 12*(1:M)^2/(M^3-M)-6*(1:M)/(M^2-1)-6*(1:M)/(M^3-M)
        
        #tmp <- double(M)
        
        #for(m in 1:M){
        #  tmp[m] <- (sdata$xg[m:N]-sdata$xl[1:(N-m+1)])%*%
        #    (sdata$ygamma[m:N]-sdata$ylambda[1:(N-m+1)])
        #}
        tmp <- .C("msrc",
                  as.integer(M),
                  as.integer(N),
                  as.double(sdata$xg),
                  as.double(sdata$xl),
                  as.double(sdata$ygamma),
                  as.double(sdata$ylambda),
                  result=double(M))$result
        
        cmat[i,j] <- (alpha/(1:M))%*%tmp
        
      }else{
        
        X <- as.numeric(data[[i]])
        
        N <- length(X)-1
        
        M <- ceiling(c.multi[i,j]*sqrt(N))
        M <- min(c(M,N)) # M must be smaller than N
        M <- max(c(M,2)) # M must be lager than 2
        
        alpha <- 12*(1:M)^2/(M^3-M)-6*(1:M)/(M^2-1)-6*(1:M)/(M^3-M)
        
        #tmp <- double(M)
        
        #for(m in 1:M){
          #tmp[m] <- sum((X[m:N]-X[1:(N-m+1)])^2)
        #  tmp[m] <- sum(diff(X,lag=m)^2)
        #}
        tmp <- .C("msrc",
                  as.integer(M),
                  as.integer(N),
                  as.double(X[-1]),
                  as.double(X[1:N]),
                  as.double(X[-1]),
                  as.double(X[1:N]),
                  result=double(M))$result
        
        cmat[i,j] <- (alpha/(1:M))%*%tmp
        
      }
      cmat[j,i] <- cmat[i,j]  
    }
    
  }
  
  return(cmat)
}

################################################################

# Realized kernel

## data: a list of zoos  
## kernel: a real-valued function on [0,infty) (a kernel)
## H: a positive number  m: a postive integer

RK <- function(data,kernel,H,c.RK,eta=3/5,m=2,ftregion=0,utime){
  
  # if missing kernel, we use the Parzen kernel
  if(missing(kernel)){
    kernel <- function(x){
      if(x<=1/2){
        return(1-6*x^2+6*x^3)
      }else if(x<=1){
        return(2*(1-x)^3)
      }else{
        return(0)
      }
    }
  }
  
  d <- length(data)
  
  tmp <- lapply(refresh_sampling(data),"as.numeric")
  #tmp <- do.call("rbind",tmp)
  tmp <- do.call("cbind",tmp)
  
  #n <- max(c(ncol(tmp)-2*m+1,2))
  n <- max(c(nrow(tmp)-2*m+1,2))
  
  # if missing H, we select it following Barndorff-Nielsen et al.(2011)
  if(missing(H)){
    if(missing(c.RK)) c.RK <- selectParam.RK(data,utime=utime)
    H <- mean(c.RK)*n^eta
  }
  
  #X <- matrix(0,d,n+1)
  X <- matrix(0,n+1,d)
  
  #X[,1] <- apply(matrix(tmp[,1:m],d,m),1,"sum")/m
  #X[,1] <- apply(matrix(tmp[,1:m],d,m),1,"mean")
  X[1,] <- colMeans(matrix(tmp[1:m,],m,d))
  #X[,2:n] <- tmp[,(m+1):(m+n-1)]
  X[2:n,] <- tmp[(m+1):(m+n-1),]
  #X[,n+1] <- apply(matrix(tmp[,(n+m):(n+2*m-1)],d,m),1,"sum")/m
  #X[,n+1] <- apply(matrix(tmp[,(n+m):(n+2*m-1)],d,m),1,"mean")
  X[n+1,] <- colMeans(matrix(tmp[(n+m):(n+2*m-1),],m,d))
  
  cc <- floor(ftregion*H) # flat-top region
  Kh <- rep(1,cc)
  Kh <- append(Kh,sapply(((cc+1):(n-1))/H,kernel))
  h.size <- max(which(Kh!=0))
  
  #diffX <- apply(X,1,FUN="diff")
  #diffX <- diff(X)
  #Gamma <- array(0,dim=c(h.size+1,d,d))
  #for(h in 1:(h.size+1))
  #  Gamma[h,,] <- t(diffX)[,h:n]%*%diffX[1:(n-h+1),]
  Gamma <- acf(diff(X),lag.max=h.size,type="covariance",
               plot=FALSE,demean=FALSE)$acf*n

  cmat <- matrix(0,d,d)
  for (i in 1:d){
    cmat[,i] <- kernel(0)*Gamma[1,,i]+
      Kh[1:h.size]%*%(Gamma[-1,,i]+Gamma[-1,i,])
  }
  
  return(cmat)
}


#############################################################

# QMLE (Ait-Sahalia et al.(2010))
## Old sources

if(0){
ql.xiu <- function(zdata){
  
  diffX <- diff(as.numeric(zdata))
  inv.difft <- diff(as.numeric(time(zdata)))^(-1)
  
  n <- length(diffX)
  a <- 4*inv.difft*sin((pi/2)*seq(1,2*n-1,by=2)/(2*n+1))^2
  
  z <- double(n)
  for(k in 1:n){
    z[k] <- cos(pi*((2*k-1)/(2*n+1))*(1:n-1/2))%*%diffX
  }
  z <- sqrt(2/(n+1/2))*z
  #z <- sqrt(2/(n+1/2))*cos(pi*((2*(1:n)-1)/(2*n+1))%o%(1:n-1/2))%*%diffX
  z2 <- inv.difft*z^2
  
  n.ql <- function(v){
    V <- v[1]+a*v[2]
    return(sum(log(V)+V^(-1)*z2))
  }
  
  gr <- function(v){
    V <- v[1]+a*v[2]
    obj <- V^(-1)-V^(-2)*z2
    return(c(sum(obj),sum(a*obj)))
  }
  
  return(list(n.ql=n.ql,gr=gr))
}

cce.qmle <- function(data,opt.method="BFGS",vol.init=NA,
                     covol.init=NA,nvar.init=NA,ncov.init=NA,...,utime){
  
  d.size <- length(data)
  dd <- d.size*(d.size-1)/2
  
  vol.init <- matrix(vol.init,1,d.size)
  nvar.init <- matrix(vol.init,1,d.size)
  covol.init <- matrix(covol.init,2,dd)
  ncov.init <- matrix(ncov.init,2,dd)
  
  if(missing(utime)) utime <- ifelse(is.numeric(time(data[[1]])),23400,1)
  
  cmat <- matrix(0,d.size,d.size)
  
  for(i in 1:d.size){
    for(j in i:d.size){
      if(i!=j){
        
        idx <- d.size*(i-1)-(i-1)*i/2+(j-i)
        
        dat <- refresh_sampling(list(data[[i]],data[[j]]))
        dattime <- apply(do.call("rbind",
                                 lapply(lapply(dat,"time"),"as.numeric")),
                                 2,"max")
        dat[[1]] <- zoo(as.numeric(dat[[1]]),dattime)
        dat[[2]] <- zoo(as.numeric(dat[[2]]),dattime)
        
        dat1 <- dat[[1]]+dat[[2]]
        dat2 <- dat[[1]]-dat[[2]]
        
        ql1 <- ql.xiu(dat1)
        ql2 <- ql.xiu(dat2)
        
        Sigma1 <- covol.init[1,idx]
        Sigma2 <- covol.init[2,idx]
        Omega1 <- ncov.init[1,idx]
        Omega2 <- ncov.init[2,idx]
      
        if(is.na(Sigma1)) Sigma1 <- RV.sparse(dat1,frequency=1200,utime=utime)
        if(is.na(Sigma2)) Sigma2 <- RV.sparse(dat2,frequency=1200,utime=utime)
        if(is.na(Omega1)) Omega1 <- Omega_BNHLS(dat1,sec=120,utime=utime)
        if(is.na(Omega2)) Omega2 <- Omega_BNHLS(dat2,sec=120,utime=utime)
        
        #par1 <- optim(par1,fn=ql.ks(dat1),method=opt.method)
        #par2 <- optim(par2,fn=ql.ks(dat2),method=opt.method)
        par1 <- constrOptim(theta=c(Sigma1,Omega1),
                            f=ql1$n.ql,grad=ql1$gr,method=opt.method,
                            ui=diag(2),ci=0,...)$par[1]
        par2 <- constrOptim(theta=c(Sigma2,Omega2),
                            f=ql2$n.ql,grad=ql2$gr,method=opt.method,
                            ui=diag(2),ci=0,...)$par[1]
        cmat[i,j] <- (par1-par2)/4
        cmat[j,i] <- cmat[i,j]
      }else{
        
        ql <- ql.xiu(data[[i]])
        
        Sigma <- vol.init[i]
        Omega <- nvar.init[i]
        
        if(is.na(Sigma)) Sigma <- RV.sparse(data[[i]],frequency=1200,utime=utime)
        if(is.na(Omega)) Omega <- Omega_BNHLS(data[[i]],sec=120,utime=utime)
        
        #cmat[i,i] <- optim(par,fn=ql.ks(data[[i]]),method=opt.method)
        cmat[i,i] <- constrOptim(theta=c(Sigma,Omega),f=ql$n.ql,grad=ql$gr,
                                 method=opt.method,
                                 ui=diag(2),ci=0,...)$par[1]
        
      }
    }
  }
  
  return(cmat)
}
}

## New sources (2014/11/10, use arima0)

cce.qmle <- function(data,vol.init=NA,covol.init=NA,nvar.init=NA,ncov.init=NA){
  
  d.size <- length(data)
  dd <- d.size*(d.size-1)/2
  
  vol.init <- matrix(vol.init,1,d.size)
  nvar.init <- matrix(vol.init,1,d.size)
  covol.init <- matrix(covol.init,2,dd)
  ncov.init <- matrix(ncov.init,2,dd)
  
  cmat <- matrix(0,d.size,d.size)
  
  for(i in 1:d.size){
    for(j in i:d.size){
      if(i!=j){
        
        idx <- d.size*(i-1)-(i-1)*i/2+(j-i)
        
        dat <- refresh_sampling(list(data[[i]],data[[j]]))
        dat[[1]] <- diff(as.numeric(dat[[1]]))
        dat[[2]] <- diff(as.numeric(dat[[2]]))
        n <- length(dat[[1]])
        
        Sigma1 <- covol.init[1,idx]
        Sigma2 <- covol.init[2,idx]
        Omega1 <- ncov.init[1,idx]
        Omega2 <- ncov.init[2,idx]
        
        init <- (-2*Omega1-Sigma1/n+sqrt(Sigma1*(4*Omega1+Sigma1/n)/n))/(2*Omega1)
        obj <- arima0(dat[[1]]+dat[[2]],order=c(0,0,1),include.mean=FALSE,
                      init=init)
        par1 <- n*obj$sigma2*(1+obj$coef)^2
        
        init <- (-2*Omega2-Sigma2/n+sqrt(Sigma2*(4*Omega2+Sigma2/n)/n))/(2*Omega2)
        obj <- arima0(dat[[1]]-dat[[2]],order=c(0,0,1),include.mean=FALSE,
                      init=init)
        par2 <- n*obj$sigma2*(1+obj$coef)^2
        
        cmat[i,j] <- (par1-par2)/4
        cmat[j,i] <- cmat[i,j]
      }else{
        
        dx <- diff(as.numeric(data[[i]]))
        n <- length(dx)
        
        Sigma <- vol.init[i]
        Omega <- nvar.init[i]
        
        init <- (-2*Omega-Sigma/n+sqrt(Sigma*(4*Omega+Sigma/n)/n))/(2*Omega)
        obj <- arima0(dx,order=c(0,0,1),include.mean=FALSE,init=init)
        cmat[i,i] <- n*obj$sigma2*(1+obj$coef)^2
        
      }
    }
  }
  
  return(cmat)
}


##################################################################

# Separating information maximum likelihood estimator (Kunitomo and Sato (2008))

SIML <- function(data,mn,alpha=0.4){
  
  d.size <- length(data)
  
  data <- lapply(refresh_sampling(data),"as.numeric")
  data <- do.call("cbind",data)
  
  n.size <- nrow(data)-1
  
  if(missing(mn)){
    mn <- ceiling(n.size^alpha)
  }
  
  #C <- matrix(1,n.size,n.size)
  #C[upper.tri(C)] <- 0
  
  #P <- matrix(0,n.size,n.size)
  
  #for(j in 1:n.size){
  #  for(k in 1:n.size){
  #    P[j,k] <- sqrt(2/(n.size+1/2))*cos((2*pi/(2*n.size+1))*(j-1/2)*(k-1/2))
  #  }
  #}
  
  #Z <- data[-1,]-matrix(1,n.size,1)%*%matrix(data[1,],1,d.size)
  #Z <- sqrt(n.size)*t(P)%*%solve(C)%*%Z
  
  #diff.Y <- apply(data,2,"diff")
  diff.Y <- diff(data)
  
  #Z <- matrix(0,n.size,d.size)
  #for(j in 1:n.size){
  #  pj <- sqrt(2/(n.size+1/2))*cos((2*pi/(2*n.size+1))*(j-1/2)*(1:n.size-1/2))
  #  Z[j,] <- matrix(pj,1,n.size)%*%diff.Y
  #}
  #Z <- matrix(0,mn,d.size)
  #for(j in 1:mn){
  #  pj <- sqrt(2/(n.size+1/2))*cos((2*pi/(2*n.size+1))*(j-1/2)*(1:n.size-1/2))
  #  Z[j,] <- matrix(pj,1,n.size)%*%diff.Y
  #}
  #Z <- sqrt(n.size)*Z
  Z <- sqrt(n.size)*
    sqrt(2/(n.size+1/2))*cos((2*pi/(2*n.size+1))*(1:mn-1/2)%o%(1:n.size-1/2))%*%
    diff.Y
  
  cmat <- matrix(0,d.size,d.size)
  for(k in 1:mn){
    cmat <- cmat+matrix(Z[k,],d.size,1)%*%matrix(Z[k,],1,d.size)
  }
  
  return(cmat/mn)
}


#############################################################

# Truncated Hayashi-Yoshida estimator
## data: a list of zoos 
## threshold: a numeric vector or a list of numeric vectors or zoos

THY <- function(data,threshold) {  
  
  n.series <- length(data)
  
  if(missing(threshold)){
    threshold <- local_univ_threshold(data,coef=5,eps=0.1)
  }else if(is.numeric(threshold)){
    threshold <- matrix(threshold,1,n.series)
  }
  
  #if(n.series <2)
  # stop("Please provide at least 2-dimensional time series")
  
  # allocate memory
  ser.X <- vector(n.series, mode="list")     # data in 'x'
  ser.times <- vector(n.series, mode="list") # time index in 'x'
  ser.diffX <- vector(n.series, mode="list") # difference of data
  ser.rho <- vector(n.series, mode="list")
  
  for(i in 1:n.series){
    # set data and time index
    ser.X[[i]] <- as.numeric(data[[i]]) # we need to transform data into numeric to avoid problem with duplicated indexes below
    ser.times[[i]] <- as.numeric(time(data[[i]]))
    # set difference of the data with truncation
    ser.diffX[[i]] <- diff( ser.X[[i]] )
    
    if(is.numeric(threshold)){
      ser.rho[[i]] <- rep(threshold[i],length(ser.diffX[[i]]))
    }else{
      ser.rho[[i]] <- as.numeric(threshold[[i]])
    }
    # thresholding
    ser.diffX[[i]][ser.diffX[[i]]^2>ser.rho[[i]]] <- 0
  }
  
  # core part of cce
  
  cmat <- matrix(0, n.series, n.series)  # cov
  for(i in 1:n.series){
    for(j in i:n.series){ 
      #I <- rep(1,n.series)
      #Checking Starting Point
      #repeat{
      #  if(ser.times[[i]][I[i]] >= ser.times[[j]][I[j]+1]){
      #    I[j] <- I[j]+1	 
      #  }else if(ser.times[[i]][I[i]+1] <= ser.times[[j]][I[j]]){
      #    I[i] <- I[i]+1	 
      #  }else{
      #    break
      #  }
      #}
      
      
      #Main Component
      if(i!=j){
        #while((I[i]<length(ser.times[[i]])) && (I[j]<length(ser.times[[j]]))) {
        #  cmat[j,i] <- cmat[j,i] + (ser.diffX[[i]])[I[i]]*(ser.diffX[[j]])[I[j]]
        #  if(ser.times[[i]][I[i]+1]>ser.times[[j]][I[j]+1]){
        #    I[j] <- I[j] + 1
        #  }else if(ser.times[[i]][I[i]+1]<ser.times[[j]][I[j]+1]){
        #    I[i] <- I[i] +1
        #  }else{
        #    I[i] <- I[i]+1
        #    I[j] <- I[j]+1
        #  }
        #}
        cmat[j,i] <- .C("HayashiYoshida",as.integer(length(ser.times[[i]])),
                        as.integer(length(ser.times[[j]])),as.double(ser.times[[i]]),
                        as.double(ser.times[[j]]),as.double(ser.diffX[[i]]),
                        as.double(ser.diffX[[j]]),value=double(1))$value
      }else{
        cmat[i,j] <- sum(ser.diffX[[i]]^2)
      }
      cmat[i,j] <- cmat[j,i]  
    }
    
  }
  
  return(cmat)
}

#########################################################

# Pre-averaged truncated Hayashi-Yoshida estimator
## data: a list of zoos  theta: a postive number
## g: a real-valued function on [0,1] (a weight function)
## threshold: a numeric vector or a list of numeric vectors or zoos
## refreshing: a logical value (if TRUE we use the refreshed data)

PTHY <- function(data,theta,kn,g,threshold,refreshing=TRUE,
                 cwise=TRUE,eps=0.2){
  
  n.series <- length(data)
  
  #if(missing(theta)&&missing(kn))
  #  theta <- selectParam.pavg(data,a.theta=7585/1161216,
  #                            b.theta=151/20160,c.theta=1/24)
  if(missing(theta)) theta <- 0.15
  
  cmat <- matrix(0, n.series, n.series)  # cov
  
  if(refreshing){
    
    if(cwise){
      
      if(missing(kn)){# refreshing,cwise,missing kn
        
        theta <- matrix(theta,n.series,n.series) 
        theta <- (theta+t(theta))/2
        
        for(i in 1:n.series){
          for(j in i:n.series){ 
            if(i!=j){
              
              resample <- refresh_sampling.PHY(list(data[[i]],data[[j]]))
              dat <- resample$data
              ser.numX <- length(resample$refresh.times)-1
              
              # set data and time index
              ser.X <- lapply(dat,"as.numeric")
              ser.times <- lapply(lapply(dat,"time"),"as.numeric")
              
              # set difference of the data
              ser.diffX <- lapply(ser.X,"diff")
              
              # pre-averaging
              kn <- min(max(ceiling(theta[i,j]*sqrt(ser.numX)),2),ser.numX)
              weight <- sapply((1:(kn-1))/kn,g)
              psi.kn <- sum(weight)
              
              #ser.barX <- list(rollapplyr(ser.diffX[[1]],width=kn-1,FUN="%*%",weight),
              #                 rollapplyr(ser.diffX[[2]],width=kn-1,FUN="%*%",weight))
              ser.barX <- list(filter(ser.diffX[[1]],rev(weight),method="c",
                                      sides=1)[(kn-1):length(ser.diffX[[1]])],
                               filter(ser.diffX[[2]],rev(weight),method="c",
                                      sides=1)[(kn-1):length(ser.diffX[[2]])])
              
              ser.num.barX <- sapply(ser.barX,"length")
              
              # thresholding
              if(missing(threshold)){

                for(ii in 1:2){
                  
                  K <- ceiling(ser.num.barX[ii]^(3/4))
                  
                  obj0 <- abs(ser.barX[[ii]])
                  obj1 <- (pi/2)*obj0[1:(ser.num.barX[ii]-kn)]*obj0[-(1:kn)]
                  if(min(K,ser.num.barX[ii]-1)<2*kn){
                    #v.hat <- (median(obj0)/0.6745)^2
                    v.hat <- mean(obj1)
                  }else{
                    v.hat <- double(ser.num.barX[ii])
                    #v.hat[1:K] <- (median(obj0[1:K])/0.6745)^2
                    v.hat[-(1:K)] <- rollmean(obj1[1:(ser.num.barX[ii]-2*kn)],k=K-2*kn+1,align="left")
                    v.hat[1:K] <- v.hat[K+1]
                  }
                  #rho <- 2*log(ser.num.barX[ii])*
                  #  (median(abs(ser.barX[[ii]])/sqrt(v.hat))/0.6745)^2*v.hat
                  rho <- 2*log(ser.num.barX[ii])^(1+eps)*v.hat
                  ser.barX[[ii]][ser.barX[[ii]]^2>rho] <- 0
                }
              }else if(is.numeric(threshold)){
                threshold <- matrix(threshold,1,n.series)
                ser.barX[[1]][ser.barX[[1]]^2>threshold[i]] <- 0
                ser.barX[[2]][ser.barX[[2]]^2>threshold[j]] <- 0
              }else{
                ser.barX[[1]][ser.barX[[1]]^2>threshold[[i]][1:ser.num.barX[1]]] <- 0
                ser.barX[[2]][ser.barX[[2]]^2>threshold[[j]][1:ser.num.barX[2]]] <- 0
              }
              
              ser.num.barX <- ser.num.barX-1
              
              # core part of cce
              #start <- kn+1
              #end <- 1
              #for(k in 1:ser.num.barX[1]){
              #  while(!(ser.times[[1]][k]<ser.times[[2]][start])&&((start-kn)<ser.num.barX[2])){
              #    start <- start + 1
              #  }
              #  while((ser.times[[1]][k+kn]>ser.times[[2]][end+1])&&(end<ser.num.barX[2])){
              #    end <- end + 1
              #  }
              #  cmat[i,j] <- cmat[i,j] + ser.barX[[1]][k]*sum(ser.barX[[2]][(start-kn):end])
              #}
              cmat[i,j] <- .C("pHayashiYoshida",
                              as.integer(kn),
                              as.integer(ser.num.barX[1]),
                              as.integer(ser.num.barX[2]),
                              as.double(ser.times[[1]]),
                              as.double(ser.times[[2]]),
                              as.double(ser.barX[[1]]),
                              as.double(ser.barX[[2]]),
                              value=double(1))$value
              
              cmat[i,j] <- cmat[i,j]/(psi.kn^2)
              cmat[j,i] <- cmat[i,j]
              
            }else{
              
              diffX <- diff(as.numeric(data[[i]]))
              
              # pre-averaging
              kn <- min(max(ceiling(theta[i,i]*sqrt(length(diffX))),2),
                        length(data[[i]]))
              weight <- sapply((1:(kn-1))/kn,g)
              psi.kn <- sum(weight)
              
              #barX <- rollapplyr(diffX,width=kn-1,FUN="%*%",weight)
              barX <- filter(diffX,rev(weight),method="c",
                             sides=1)[(kn-1):length(diffX)]
              num.barX <- length(barX)
              
              # thresholding
              if(missing(threshold)){
                K <- ceiling(num.barX^(3/4))
                #K <- ceiling(kn^(3/2))
                
                obj0 <- abs(barX)
                obj1 <- (pi/2)*obj0[1:(num.barX-kn)]*obj0[-(1:kn)]
                if(min(K,num.barX-1)<2*kn){
                  #v.hat <- (median(obj0)/0.6745)^2
                  v.hat <- mean(obj1)
                }else{
                  v.hat <- double(num.barX)
                  #v.hat[1:K] <- (median(obj0[1:K])/0.6745)^2
                  v.hat[-(1:K)] <- rollmean(obj1[1:(num.barX-2*kn)],k=K-2*kn+1,align="left")
                  v.hat[1:K] <- v.hat[K+1]
                }
                #rho <- 2*log(num.barX)*
                #  (median(abs(barX)/sqrt(v.hat))/0.6745)^2*v.hat
                rho <- 2*log(num.barX)^(1+eps)*v.hat
                barX[barX^2>rho] <- 0
              }else if(is.numeric(threshold)){
                threshold <- matrix(threshold,1,n.series)
                barX[barX^2>threshold[i]] <- 0
              }else{
                barX[barX^2>threshold[[i]][1:num.barX]] <- 0
              }
              
              tmp <- barX[-num.barX]
              #cmat[i,j] <- tmp%*%rollapply(tmp,width=2*(kn-1)+1,FUN="sum",partial=TRUE)/(psi.kn)^2
              cmat[i,j] <- tmp%*%rollsum(c(double(kn-1),tmp,double(kn-1)),
                                         k=2*(kn-1)+1)/(psi.kn)^2
              
            }
          }
        }
      }else{# refreshing,cwise,not missing kn
        
        kn <- matrix(kn,n.series,n.series)
        
        for(i in 1:n.series){
          for(j in i:n.series){
            
            if(i!=j){
              
              resample <- refresh_sampling.PHY(list(data[[i]],data[[j]]))
              dat <- resample$data
              ser.numX <- length(resample$refresh.times)-1
              
              knij <- min(max(kn[i,j],2),ser.numX+1)
              weight <- sapply((1:(knij-1))/knij,g)
              psi.kn <- sum(weight)
              
              # set data and time index
              ser.X <- lapply(dat,"as.numeric")
              ser.times <- lapply(lapply(dat,"time"),"as.numeric")
              
              # set difference of the data
              ser.diffX <- lapply(ser.X,"diff")
              
              # pre-averaging
              #ser.barX <- list(rollapplyr(ser.diffX[[1]],width=knij-1,FUN="%*%",weight),
              #                 rollapplyr(ser.diffX[[2]],width=knij-1,FUN="%*%",weight))
              ser.barX <- list(filter(ser.diffX[[1]],rev(weight),method="c",
                                      sides=1)[(knij-1):length(ser.diffX[[1]])],
                               filter(ser.diffX[[2]],rev(weight),method="c",
                                      sides=1)[(knij-1):length(ser.diffX[[2]])])
              
              ser.num.barX <- sapply(ser.barX,"length")
              
              # thresholding
              if(missing(threshold)){
                for(ii in 1:2){
                  
                  K <- ceiling(ser.num.barX[ii]^(3/4))
                  
                  obj0 <- abs(ser.barX[[ii]])
                  obj1 <- (pi/2)*obj0[1:(ser.num.barX[ii]-knij)]*obj0[-(1:knij)]
                  if(min(K,ser.num.barX[ii]-1)<2*knij){
                    #v.hat <- (median(obj0)/0.6745)^2
                    v.hat <- mean(obj1)
                  }else{
                    v.hat <- double(ser.num.barX[ii])
                    #v.hat[1:K] <- (median(obj0[1:K])/0.6745)^2
                    v.hat[-(1:K)] <- rollmean(obj1[1:(ser.num.barX[ii]-2*knij)],k=K-2*knij+1,align="left")
                    v.hat[1:K] <- v.hat[K+1]
                  }
                  #rho <- 2*log(ser.num.barX[ii])*
                  #  (median(abs(ser.barX[[ii]])/sqrt(v.hat))/0.6745)^2*v.hat
                  rho <- 2*log(ser.num.barX[ii])^(1+eps)*v.hat
                  ser.barX[[ii]][ser.barX[[ii]]^2>rho] <- 0
                }
              }else if(is.numeric(threshold)){
                threshold <- matrix(threshold,1,n.series)
                ser.barX[[1]][ser.barX[[1]]^2>threshold[i]] <- 0
                ser.barX[[2]][ser.barX[[2]]^2>threshold[j]] <- 0
              }else{
                ser.barX[[1]][ser.barX[[1]]^2>threshold[[i]][1:ser.num.barX[1]]] <- 0
                ser.barX[[2]][ser.barX[[2]]^2>threshold[[j]][1:ser.num.barX[2]]] <- 0
              }
              
              ser.num.barX <- ser.num.barX-1
              
              # core part of cce
              #start <- knij+1
              #end <- 1
              #for(k in 1:ser.num.barX[1]){
              #  while(!(ser.times[[1]][k]<ser.times[[2]][start])&&((start-knij)<ser.num.barX[2])){
              #    start <- start + 1
              #  }
              #  while((ser.times[[1]][k+knij]>ser.times[[2]][end+1])&&(end<ser.num.barX[2])){
              #    end <- end + 1
              #  }
              #  cmat[i,j] <- cmat[i,j] + ser.barX[[1]][k]*sum(ser.barX[[2]][(start-knij):end])
              #}
              cmat[i,j] <- .C("pHayashiYoshida",
                              as.integer(knij),
                              as.integer(ser.num.barX[1]),
                              as.integer(ser.num.barX[2]),
                              as.double(ser.times[[1]]),
                              as.double(ser.times[[2]]),
                              as.double(ser.barX[[1]]),
                              as.double(ser.barX[[2]]),
                              value=double(1))$value
              
              cmat[i,j] <- cmat[i,j]/(psi.kn^2)
              cmat[j,i] <- cmat[i,j]
              
            }else{
              
              diffX <- diff(as.numeric(data[[i]]))
              # pre-averaging
              kni <- min(max(kni,2),length(data[[i]]))
              weight <- sapply((1:(kni-1))/kni,g)
              psi.kn <- sum(weight)
              
              #barX <- rollapplyr(diffX,width=kni-1,FUN="%*%",weight)
              barX <- filter(diffX,rev(weight),method="c",
                             sides=1)[(kni-1):length(diffX)]
              num.barX <- length(barX)
              
              # thresholding
              if(missing(threshold)){
                
                K <- ceiling(num.barX^(3/4))
                #K <- ceiling(kn^(3/2))
                
                obj0 <- abs(barX)
                obj1 <- (pi/2)*obj0[1:(num.barX-kni)]*obj0[-(1:kni)]
                if(min(K,num.barX-1)<2*kni){
                  #v.hat <- (median(obj0)/0.6745)^2
                  v.hat <- mean(obj1)
                }else{
                  v.hat <- double(num.barX)
                  #v.hat[1:K] <- (median(obj0[1:K])/0.6745)^2
                  v.hat[-(1:K)] <- rollmean(obj1[1:(num.barX-2*kni)],k=K-2*kni+1,align="left")
                  v.hat[1:K] <- v.hat[K+1]
                }
                #rho <- 2*log(num.barX)*
                #  (median(abs(barX)/sqrt(v.hat))/0.6745)^2*v.hat
                rho <- 2*log(num.barX)^(1+eps)*v.hat
                barX[barX^2>rho] <- 0
              }else if(is.numeric(threshold)){
                threshold <- matrix(threshold,1,n.series)
                barX[barX^2>threshold[i]] <- 0
              }else{
                barX[barX^2>threshold[[i]][1:num.barX]] <- 0
              }
              
              tmp <- barX[-num.barX]
              #cmat[i,j] <- tmp%*%rollapply(tmp,width=2*(kni-1)+1,FUN="sum",partial=TRUE)/(psi.kn)^2
              cmat[i,j] <- tmp%*%rollsum(c(double(kni-1),tmp,double(kni-1)),
                                         k=2*(kni-1)+1)/(psi.kn)^2
              
            }
          }
        }
      }
    }else{# refreshing,non-cwise
      
      # synchronize the data
      resample <- refresh_sampling.PHY(data)
      data <- resample$data
      ser.numX <- length(resample$refresh.times)-1
      
      # if missing kn, we choose it following Barndorff-Nielsen et al. (2011)
      if(missing(kn)){
        kn <- min(max(ceiling(mean(theta)*sqrt(ser.numX)),2),ser.numX+1)
      }
      kn <- kn[1]
      
      weight <- sapply((1:(kn-1))/kn,g)
      
      # allocate memory
      ser.X <- vector(n.series, mode="list")     # data in 'x'
      ser.times <- vector(n.series, mode="list") # time index in 'x'
      ser.diffX <- vector(n.series, mode="list") # difference of data
      ser.barX <- vector(n.series, mode="list")  # pre-averaged data
      ser.num.barX <- integer(n.series)         
      
      for(i in 1:n.series){
        # set data and time index
        ser.X[[i]] <- as.numeric(data[[i]]) # we need to transform data into numeric to avoid problem with duplicated indexes below
        ser.times[[i]] <- as.numeric(time(data[[i]]))
        # set difference of the data 
        ser.diffX[[i]] <- diff( ser.X[[i]] )
      }
      
      # thresholding
      if(missing(threshold)){
        #coef <- 2*coef*kn/sqrt(ser.numX)
        for(i in 1:n.series){
          
          #ser.num.barX[i] <- ser.numX[i]-kn+2
          #ser.barX[[i]] <- double(ser.num.barX[i])
          
          #for(j in 1:ser.num.barX[i]){
          #  ser.barX[[i]][j] <- sapply((1:(kn-1))/kn,g)%*%ser.diffX[[i]][j:(j+kn-2)]
          #}
          
          #ser.barX[[i]] <- rollapplyr(ser.diffX[[i]],width=kn-1,FUN="%*%",weight)
          ser.barX[[i]] <- filter(ser.diffX[[i]],rev(weight),method="c",
                                  sides=1)[(kn-1):length(ser.diffX[[i]])]
          ser.num.barX[i] <- length(ser.barX[[i]])
          
          K <- ceiling(ser.num.barX[i]^(3/4))
          
          obj0 <- abs(ser.barX[[i]])
          obj1 <- (pi/2)*obj0[1:(ser.num.barX[i]-kn)]*obj0[-(1:kn)]
          if(min(K,ser.num.barX[i]-1)<2*kn){
            #v.hat <- (median(obj0)/0.6745)^2
            v.hat <- mean(obj1)
          }else{
            v.hat <- double(ser.num.barX[i])
            #v.hat[1:K] <- (median(obj0[1:K])/0.6745)^2
            v.hat[-(1:K)] <- rollmean(obj1[1:(ser.num.barX[i]-2*kn)],
                                      k=K-2*kn+1,align="left")
            v.hat[1:K] <- v.hat[K+1]
          }
          #rho <- 2*log(ser.num.barX[ii])*
          #  (median(abs(ser.barX[[ii]])/sqrt(v.hat))/0.6745)^2*v.hat
          rho <- 2*log(ser.num.barX[i])^(1+eps)*v.hat
          ser.barX[[i]][ser.barX[[i]]^2>rho] <- 0
        }
      }else if(is.numeric(threshold)){
        threshold <- matrix(threshold,1,n.series)
        for(i in 1:n.series){
          #ser.num.barX[i] <- ser.numX[i]-kn+2
          #ser.barX[[i]] <- double(ser.num.barX[i])
          #for(j in 1:ser.num.barX[i]){
          #  tmp <- sapply((1:(kn-1))/kn,g)%*%ser.diffX[[i]][j:(j+kn-2)]
          #  if(tmp^2>threshold[i]){
          #    ser.barX[[i]][j] <- 0
          #  }else{
          #    ser.barX[[i]][j] <- tmp
          #  }
          #}
          #ser.barX[[i]] <- rollapplyr(ser.diffX[[i]],width=kn-1,FUN="%*%",weight)
          ser.barX[[i]] <- filter(ser.diffX[[i]],rev(weight),method="c",
                                  sides=1)[(kn-1):length(ser.diffX[[i]])]
          ser.num.barX[i] <- length(ser.barX[[i]])
          ser.barX[[i]][ser.barX[[i]]^2>threshold[i]] <- 0
        }
      }else{
        for(i in 1:n.series){
          #ser.num.barX[i] <- ser.numX[i]-kn+2
          #ser.barX[[i]] <- double(ser.num.barX[i])
          #for(j in 1:ser.num.barX[i]){
          #  tmp <- sapply((1:(kn-1))/kn,g)%*%ser.diffX[[i]][j:(j+kn-2)]
          #  if(tmp^2>threshold[[i]][j]){
          #    ser.barX[[i]][j] <- 0
          #  }else{
          #    ser.barX[[i]][j] <- tmp
          #  }
          #}
          #ser.barX[[i]] <- rollapplyr(ser.diffX[[i]],width=kn-1,FUN="%*%",weight)
          ser.barX[[i]] <- filter(ser.diffX[[i]],rev(weight),method="c",
                                  sides=1)[(kn-1):length(ser.diffX[[i]])]
          ser.num.barX[i] <- length(ser.barX[[i]])
          ser.barX[[i]][ser.barX[[i]]^2>threshold[[i]][1:ser.num.barX[i]]] <- 0
        }
      }
      
      ser.num.barX <- ser.num.barX-1
      
      # core part of cce
      
      for(i in 1:n.series){
        for(j in i:n.series){ 
          if(i!=j){
            #start <- kn+1
            #end <- 1
            #for(k in 1:ser.num.barX[i]){
            #  while(!(ser.times[[i]][k]<ser.times[[j]][start])&&((start-kn)<ser.num.barX[j])){
            #    start <- start + 1
            #  }
            #  while((ser.times[[i]][k+kn]>ser.times[[j]][end+1])&&(end<ser.num.barX[j])){
            #    end <- end + 1
            #  }
            #  cmat[i,j] <- cmat[i,j] + ser.barX[[i]][k]*sum(ser.barX[[j]][(start-kn):end])
            #}
            cmat[i,j] <- .C("pHayashiYoshida",
                            as.integer(kn),
                            as.integer(ser.num.barX[i]),
                            as.integer(ser.num.barX[j]),
                            as.double(ser.times[[i]]),
                            as.double(ser.times[[j]]),
                            as.double(ser.barX[[i]]),
                            as.double(ser.barX[[j]]),
                            value=double(1))$value
            cmat[j,i] <- cmat[i,j]
          }else{
            tmp <- ser.barX[[i]][1:ser.num.barX[i]]
            #cmat[i,j] <- tmp%*%rollapply(tmp,width=2*(kn-1)+1,FUN="sum",partial=TRUE)
            cmat[i,j] <- tmp%*%rollsum(c(double(kn-1),tmp,double(kn-1)),
                                       k=2*(kn-1)+1) 
          }
        }
      }
      
      psi.kn <- sum(weight)
      cmat <- cmat/(psi.kn^2)
      
    }
  }else{# non-refreshing
    
    if(cwise){# non-refreshing,cwise
      
      theta <- matrix(theta,n.series,n.series) 
      theta <- (theta+t(theta))/2
      
      if(missing(kn)){
        ntmp <- matrix(sapply(data,"length"),n.series,n.series)
        ntmp <- ntmp+t(ntmp)
        diag(ntmp) <- diag(ntmp)/2
        kn <- ceiling(theta*sqrt(ntmp))
        kn[kn<2] <- 2 # kn must be larger than 1
      }
      kn <- matrix(kn,n.series,n.series)
      
      for(i in 1:n.series){
        for(j in i:n.series){ 
          if(i!=j){
            
            dat <- list(data[[i]],data[[j]])
            
            # set data and time index
            ser.X <- lapply(dat,"as.numeric")
            ser.times <- lapply(lapply(dat,"time"),"as.numeric")
            # set difference of the data
            ser.diffX <- lapply(ser.X,"diff")
            
            # pre-averaging
            knij <- min(kn[i,j],sapply(ser.X,"length")) # kn must be less than the numbers of the observations
            weight <- sapply((1:(knij-1))/knij,g)
            psi.kn <- sum(weight)
            
            #ser.barX <- list(rollapplyr(ser.diffX[[1]],width=knij-1,FUN="%*%",weight),
            #                 rollapplyr(ser.diffX[[2]],width=knij-1,FUN="%*%",weight))
            ser.barX <- list(filter(ser.diffX[[1]],rev(weight),method="c",
                                    sides=1)[(knij-1):length(ser.diffX[[1]])],
                             filter(ser.diffX[[2]],rev(weight),method="c",
                                    sides=1)[(knij-1):length(ser.diffX[[2]])])
            ser.num.barX <- sapply(ser.barX,"length")
            
            # thresholding
            if(missing(threshold)){
              
              for(ii in 1:2){
                
                K <- ceiling(ser.num.barX[ii]^(3/4))
                
                obj0 <- abs(ser.barX[[ii]])
                obj1 <- (pi/2)*obj0[1:(ser.num.barX[ii]-knij)]*obj0[-(1:knij)]
                if(min(K,ser.num.barX[ii]-1)<2*knij){
                  #v.hat <- (median(obj0)/0.6745)^2
                  v.hat <- mean(obj1)
                }else{
                  v.hat <- double(ser.num.barX[ii])
                  #v.hat[1:K] <- (median(obj0[1:K])/0.6745)^2
                  v.hat[-(1:K)] <- rollmean(obj1[1:(ser.num.barX[ii]-2*knij)],k=K-2*knij+1,align="left")
                  v.hat[1:K] <- v.hat[K+1]
                }
                #rho <- 2*log(ser.num.barX[ii])*
                #  (median(abs(ser.barX[[ii]])/sqrt(v.hat))/0.6745)^2*v.hat
                rho <- 2*log(ser.num.barX[ii])^(1+eps)*v.hat
                ser.barX[[ii]][ser.barX[[ii]]^2>rho] <- 0
              }
            }else if(is.numeric(threshold)){
              threshold <- matrix(threshold,1,n.series)
              ser.barX[[1]][ser.barX[[1]]^2>threshold[i]] <- 0
              ser.barX[[2]][ser.barX[[2]]^2>threshold[j]] <- 0
            }else{
              ser.barX[[1]][ser.barX[[1]]^2>threshold[[i]][1:ser.num.barX[1]]] <- 0
              ser.barX[[2]][ser.barX[[2]]^2>threshold[[j]][1:ser.num.barX[2]]] <- 0
            }
            
            ser.num.barX <- ser.num.barX-1
            
            # core part of cce
            #start <- knij+1
            #end <- 1
            #for(k in 1:ser.num.barX[1]){
            #  while(!(ser.times[[1]][k]<ser.times[[2]][start])&&((start-knij)<ser.num.barX[2])){
            #    start <- start + 1
            #  }
            #  while((ser.times[[1]][k+knij]>ser.times[[2]][end+1])&&(end<ser.num.barX[2])){
            #    end <- end + 1
            #  }
            #  cmat[i,j] <- cmat[i,j] + ser.barX[[1]][k]*sum(ser.barX[[2]][(start-knij):end])
            #}
            cmat[i,j] <- .C("pHayashiYoshida",
                            as.integer(knij),
                            as.integer(ser.num.barX[1]),
                            as.integer(ser.num.barX[2]),
                            as.double(ser.times[[1]]),
                            as.double(ser.times[[2]]),
                            as.double(ser.barX[[1]]),
                            as.double(ser.barX[[2]]),
                            value=double(1))$value
            
            cmat[i,j] <- cmat[i,j]/(psi.kn^2)
            cmat[j,i] <- cmat[i,j]
            
          }else{
            
            diffX <- diff(as.numeric(data[[i]]))
            
            # pre-averaging
            kni <- min(kn[i,i],length(data[[i]])) # kn must be less than the number of the observations
            weight <- sapply((1:(kni-1))/kni,g)
            psi.kn <- sum(weight)
            
            #barX <- rollapplyr(diffX,width=kni-1,FUN="%*%",weight)
            barX <- filter(diffX,rev(weight),method="c",
                           sides=1)[(kni-1):length(diffX)]
            num.barX <- length(barX)
            
            # thrsholding
            if(missing(threshold)){
              
              K <- ceiling(num.barX^(3/4))
              #K <- ceiling(kn^(3/2))
              
              obj0 <- abs(barX)
              obj1 <- (pi/2)*obj0[1:(num.barX-kni)]*obj0[-(1:kni)]
              if(min(K,num.barX-1)<2*kni){
                #v.hat <- (median(obj0)/0.6745)^2
                v.hat <- mean(obj1)
              }else{
                v.hat <- double(num.barX)
                #v.hat[1:K] <- (median(obj0[1:K])/0.6745)^2
                v.hat[-(1:K)] <- rollmean(obj1[1:(num.barX-2*kni)],k=K-2*kni+1,align="left")
                v.hat[1:K] <- v.hat[K+1]
              }
              #rho <- 2*log(num.barX)*
              #  (median(abs(barX)/sqrt(v.hat))/0.6745)^2*v.hat
              rho <- 2*log(num.barX)^(1+eps)*v.hat
              barX[barX^2>rho] <- 0
            }else if(is.numeric(threshold)){
              threshold <- matrix(threshold,1,n.series)
              barX[barX^2>threshold[i]] <- 0
            }else{
              barX[barX^2>threshold[[i]][1:num.barX]] <- 0
            }
            
            tmp <- barX[-num.barX]
            #cmat[i,j] <- tmp%*%rollapply(tmp,width=2*(kni-1)+1,FUN="sum",partial=TRUE)/(psi.kn)^2
            cmat[i,j] <- tmp%*%rollsum(c(double(kni-1),tmp,double(kni-1)),
                                       k=2*(kni-1)+1)/(psi.kn)^2
            
          }
        }
      }
    }else{# non-refreshing, non-cwise
      
      # allocate memory
      ser.X <- vector(n.series, mode="list")     # data in 'x'
      ser.times <- vector(n.series, mode="list") # time index in 'x'
      ser.diffX <- vector(n.series, mode="list") # difference of data
      
      ser.numX <- integer(n.series)
      ser.barX <- vector(n.series, mode="list")
      
      for(i in 1:n.series){
        # set data and time index
        ser.X[[i]] <- as.numeric(data[[i]]) # we need to transform data into numeric to avoid problem with duplicated indexes below
        ser.times[[i]] <- as.numeric(time(data[[i]]))
        
        # set difference of the data 
        ser.diffX[[i]] <- diff( ser.X[[i]] )    
        ser.numX[i] <- length(ser.diffX[[i]])
      }
      
      # if missing kn, we select it following Barndorff-Nielsen et al.(2011)  
      if(missing(kn)){
        kn <- min(max(ceiling(mean(theta)*sqrt(sum(ser.numX))),2),
                  ser.numX+1)
      }
      kn <- kn[1]
      
      weight <- sapply((1:(kn-1))/kn,g)
      psi.kn <- sum(weight)
      
      ser.num.barX <- integer(n.series)
      
      # pre-averaging and thresholding
      if(missing(threshold)){
        
        for(i in 1:n.series){
          
          #ser.barX[[i]] <- rollapplyr(ser.diffX[[i]],width=kn-1,FUN="%*%",weight)
          ser.barX[[i]] <- filter(ser.diffX[[i]],rev(weight),method="c",
                                  sides=1)[(kn-1):length(ser.diffX[[i]])]
          ser.num.barX[i] <- length(ser.barX[[i]])
          
          K <- ceiling(ser.num.barX[i]^(3/4))
          
          obj0 <- abs(ser.barX[[i]])
          obj1 <- (pi/2)*obj0[1:(ser.num.barX[i]-kn)]*obj0[-(1:kn)]
          if(min(K,ser.num.barX[i]-1)<2*kn){
            #v.hat <- (median(obj0)/0.6745)^2
            v.hat <- mean(obj1)
          }else{
            v.hat <- double(ser.num.barX[i])
            #v.hat[1:K] <- (median(obj0[1:K])/0.6745)^2
            v.hat[-(1:K)] <- rollmean(obj1[1:(ser.num.barX[i]-2*kn)],k=K-2*kn+1,align="left")
            v.hat[1:K] <- v.hat[K+1]
          }
          #rho <- 2*log(ser.num.barX[ii])*
          #  (median(abs(ser.barX[[ii]])/sqrt(v.hat))/0.6745)^2*v.hat
          rho <- 2*log(ser.num.barX[i])^(1+eps)*v.hat
          ser.barX[[i]][ser.barX[[i]]^2>rho] <- 0
        }
      }else if(is.numeric(threshold)){
        threshold <- matrix(threshold,1,n.series)
        for(i in 1:n.series){
          ser.barX[[i]] <- rollapplyr(ser.diffX[[i]],width=kn-1,FUN="%*%",weight)
          ser.num.barX[i] <- length(ser.barX[[i]])
          ser.barX[[i]][ser.barX[[i]]^2>threshold[i]] <- 0
        }
      }else{
        for(i in 1:n.series){
          ser.barX[[i]] <- rollapplyr(ser.diffX[[i]],width=kn-1,FUN="%*%",weight)
          ser.num.barX[i] <- length(ser.barX[[i]])
          ser.barX[[i]][ser.barX[[i]]^2>threshold[[i]][1:ser.num.barX[i]]] <- 0
        }
      }
      
      ser.num.barX <- ser.num.barX-1
      
      # core part of cce
      
      cmat <- matrix(0, n.series, n.series)  # cov
      for(i in 1:n.series){
        for(j in i:n.series){ 
          if(i!=j){
            #start <- kn+1
            #end <- 1
            #for(k in 1:ser.num.barX[i]){
            #  while(!(ser.times[[i]][k]<ser.times[[j]][start])&&((start-kn)<ser.num.barX[j])){
            #    start <- start + 1
            #  }
            #  while((ser.times[[i]][k+kn]>ser.times[[j]][end+1])&&(end<ser.num.barX[j])){
            #    end <- end + 1
            #  }
            #  cmat[i,j] <- cmat[i,j] + ser.barX[[i]][k]*sum(ser.barX[[j]][(start-kn):end])
            #}
            cmat[i,j] <- .C("pHayashiYoshida",
                            as.integer(kn),
                            as.integer(ser.num.barX[i]),
                            as.integer(ser.num.barX[j]),
                            as.double(ser.times[[i]]),
                            as.double(ser.times[[j]]),
                            as.double(ser.barX[[i]]),
                            as.double(ser.barX[[j]]),
                            value=double(1))$value
            cmat[j,i] <- cmat[i,j]
          }else{
            tmp <- ser.barX[[i]][1:ser.num.barX[i]]
            #cmat[i,j] <- tmp%*%rollapply(tmp,width=2*(kn-1)+1,FUN="sum",partial=TRUE) 
            cmat[i,j] <- tmp%*%rollsum(c(double(kn-1),tmp,double(kn-1)),
                                       k=2*(kn-1)+1)
          }
        }
      }
      
      cmat <- cmat/(psi.kn^2)
      
    }
  }
  
  return(cmat)
}


###################################################################

# subsampled realized covariance

SRC <- function(data,frequency=300,avg=TRUE,utime){
  
  d.size <- length(data)
  
  if(missing(utime)) utime <- ifelse(is.numeric(time(data[[1]])),23400,1)
  
  # allocate memory
  ser.X <- vector(d.size, mode="list")     # data in 'x'
  ser.times <- vector(d.size, mode="list") # time index in 'x'
  ser.numX <- double(d.size)
  
  for(d in 1:d.size){
    # set data and time index
    ser.X[[d]] <- as.numeric(data[[d]]) # we need to transform data into numeric to avoid problem with duplicated indexes below
    ser.times[[d]] <- as.numeric(time(data[[d]]))*utime
    ser.numX[d] <- length(ser.X[[d]])
  }
  
  Init <- max(sapply(ser.times,FUN="head",n=1))
  Terminal <- max(sapply(ser.times,FUN="tail",n=1))
  
  grid <- seq(Init,Terminal,by=frequency)
  n.sparse <- length(grid)
  #I <- matrix(1,d.size,n.sparse)
  
  K <- floor(Terminal-grid[n.sparse]) + 1
  
  sdiff1 <- array(0,dim=c(d.size,n.sparse-1,K))
  sdiff2 <- array(0,dim=c(d.size,n.sparse-2,frequency-K))
  
  for(d in 1:d.size){
    
    subsample <- matrix(.C("ctsubsampling",as.double(ser.X[[d]]),
                           as.double(ser.times[[d]]),
                           as.integer(frequency),as.integer(n.sparse),
                           as.integer(ser.numX[d]),as.double(grid),
                           result=double(frequency*n.sparse))$result,
                        n.sparse,frequency)
    
    sdiff1[d,,] <- diff(subsample[,1:K])
    sdiff2[d,,] <- diff(subsample[-n.sparse,-(1:K)])
    
  }
  
  if(avg){
    
    #cmat <- matrix(0,d.size,d.size)
    
    #for(t in 1:frequency){
    #  sdiff <- matrix(0,d.size,n.sparse-1)
    #  for(d in 1:d.size){
    #    for(i in 1:n.sparse){
    #      while((ser.times[[d]][I[d,i]+1]<=grid[i])&&(I[d,i]<ser.numX[d])){
    #        I[d,i] <- I[d,i]+1
    #      }
    #    }
    #    sdiff[d,] <- diff(ser.X[[d]][I[d,]])
    #  }
    #  cmat <- cmat + sdiff%*%t(sdiff)
      
    #  grid <- grid+rep(1,n.sparse)
    #  if(grid[n.sparse]>Terminal){
    #    grid <- grid[-n.sparse]
        #I <- I[,-n.sparse]
    #    n.sparse <- n.sparse-1
    #    I <- matrix(I[,-n.sparse],d.size,n.sparse)
    #  }
    #}
    
    #cmat <- cmat/frequency
    cmat <- matrix(rowMeans(cbind(apply(sdiff1,3,FUN=function(x) x %*% t(x)),
                                  apply(sdiff2,3,FUN=function(x) x %*% t(x)))),
                   d.size,d.size)
    
  }else{
    
    #sdiff <- matrix(0,d.size,n.sparse-1)
    #for(d in 1:d.size){
    #  for(i in 1:n.sparse){
    #    while((ser.times[[d]][I[d,i]+1]<=grid[i])&&(I[d,i]<ser.numX[d])){
    #      I[d,i] <- I[d,i]+1
    #    }
    #  }
    #  sdiff[d,] <- diff(ser.X[[d]][I[d,]])
    #}
    
    #cmat <- sdiff%*%t(sdiff)
    if(K>0){
      cmat <- sdiff1[,,1]%*%t(sdiff1[,,1])
    }else{
      cmat <- sdiff2[,,1]%*%t(sdiff2[,,1])
    }
  }
  
  return(cmat)
}


##############################################################

# subsampled realized bipower covariation

BPC <- function(sdata){
  
  d.size <- nrow(sdata)
  cmat <- matrix(0,d.size,d.size)
  
  for(i in 1:d.size){
    for(j in i:d.size){
      if(i!=j){
        cmat[i,j] <- (BPV(sdata[i,]+sdata[j,])-BPV(sdata[i,]-sdata[j,]))/4
        cmat[j,i] <- cmat[i,j]
      }else{
        cmat[i,i] <- BPV(sdata[i,])
      }
    }
  }
  
  return(cmat)
}

SBPC <- function(data,frequency=300,avg=TRUE,utime){
  
  d.size <- length(data)
  
  if(missing(utime)) utime <- ifelse(is.numeric(time(data[[1]])),23400,1)
  
  # allocate memory
  ser.X <- vector(d.size, mode="list")     # data in 'x'
  ser.times <- vector(d.size, mode="list") # time index in 'x'
  ser.numX <- double(d.size)
  
  for(d in 1:d.size){
    # set data and time index
    ser.X[[d]] <- as.numeric(data[[d]]) # we need to transform data into numeric to avoid problem with duplicated indexes below
    ser.times[[d]] <- as.numeric(time(data[[d]]))*utime
    ser.numX[d] <- length(ser.X[[d]])
  }
  
  Init <- max(sapply(ser.times,FUN="head",n=1))
  Terminal <- max(sapply(ser.times,FUN="tail",n=1))
  
  grid <- seq(Init,Terminal,by=frequency)
  n.sparse <- length(grid)
  #I <- matrix(1,d.size,n.sparse)
  
  K <- floor(Terminal-grid[n.sparse]) + 1
  
  sdata1 <- array(0,dim=c(d.size,n.sparse,K))
  sdata2 <- array(0,dim=c(d.size,n.sparse-1,frequency-K))
  
  for(d in 1:d.size){
    
    subsample <- matrix(.C("ctsubsampling",as.double(ser.X[[d]]),
                           as.double(ser.times[[d]]),
                           as.integer(frequency),as.integer(n.sparse),
                           as.integer(ser.numX[d]),as.double(grid),
                           result=double(frequency*n.sparse))$result,
                        n.sparse,frequency)
    
    sdata1[d,,] <- subsample[,1:K]
    sdata2[d,,] <- subsample[-n.sparse,-(1:K)]
    
  }
  
  if(avg){
    
    cmat <- matrix(0,d.size,d.size)
    
    #for(t in 1:frequency){
    #  sdata <- matrix(0,d.size,n.sparse)
    #  for(d in 1:d.size){
    #    for(i in 1:n.sparse){
    #      while((ser.times[[d]][I[d,i]+1]<=grid[i])&&(I[d,i]<ser.numX[d])){
    #        I[d,i] <- I[d,i]+1
    #      }
    #    }
    #    sdata[d,] <- ser.X[[d]][I[d,]]
    #  }
    #  cmat <- cmat + BPC(sdata)
      
    #  grid <- grid+rep(1,n.sparse)
    #  if(grid[n.sparse]>Terminal){
    #    grid <- grid[-n.sparse]
        #I <- I[,-n.sparse]
    #    n.sparse <- n.sparse-1
    #    I <- matrix(I[,-n.sparse],d.size,n.sparse)
    #  }
    #}
    
    #cmat <- cmat/frequency
    cmat <- matrix(rowMeans(cbind(apply(sdata1,3,FUN=BPC),
                                  apply(sdata2,3,FUN=BPC))),
                   d.size,d.size)
    
  }else{
    
    #sdata <- matrix(0,d.size,n.sparse)
    #for(d in 1:d.size){
    #  for(i in 1:n.sparse){
    #    while((ser.times[[d]][I[d,i]+1]<=grid[i])&&(I[d,i]<ser.numX[d])){
    #      I[d,i] <- I[d,i]+1
    #    }
    #  }
    #  sdata[d,] <- ser.X[[d]][I[d,]]
    #}
    
    #cmat <- BPC(sdata)
    if(K>0){
      cmat <- BPC(sdata1[,,1])
    }else{
      cmat <- BPC(sdata2[,,1])
    }
  }
  
  return(cmat)
}

#
# CumulativeCovarianceEstimator
#

# returns a matrix of var-cov estimates

setGeneric("cce",
           function(x,method="HY",theta,kn,g=function(x)min(x,1-x),
                    refreshing=TRUE,cwise=TRUE,
                    delta=0,adj=TRUE,K,c.two,J=1,c.multi,
                    kernel,H,c.RK,eta=3/5,m=2,ftregion=0,
                    vol.init=NA,covol.init=NA,nvar.init=NA,ncov.init=NA,
                    mn,alpha=0.4,frequency=300,avg=TRUE,
                    threshold,utime,psd=FALSE)
             standardGeneric("cce"))

setMethod("cce",signature(x="yuima"),
          function(x,method="HY",theta,kn,g=function(x)min(x,1-x),
                   refreshing=TRUE,cwise=TRUE,
                   delta=0,adj=TRUE,K,c.two,J=1,c.multi,
                   kernel,H,c.RK,eta=3/5,m=2,ftregion=0,
                   vol.init=NA,covol.init=NA,
                   nvar.init=NA,ncov.init=NA,mn,alpha=0.4,
                   frequency=300,avg=TRUE,threshold,utime,psd=FALSE)
            cce(x@data,method=method,theta=theta,kn=kn,g=g,refreshing=refreshing,
                cwise=cwise,delta=delta,adj=adj,K=K,c.two=c.two,J=J,
                c.multi=c.multi,kernel=kernel,H=H,c.RK=c.RK,eta=eta,m=m,
                ftregion=ftregion,vol.init=vol.init,
                covol.init=covol.init,nvar.init=nvar.init,
                ncov.init=ncov.init,mn=mn,alpha=alpha,
                frequency=frequency,avg=avg,threshold=threshold,
                utime=utime,psd=psd))

setMethod("cce",signature(x="yuima.data"),
          function(x,method="HY",theta,kn,g=function(x)min(x,1-x),
                   refreshing=TRUE,cwise=TRUE,
                   delta=0,adj=TRUE,K,c.two,J=1,c.multi,
                   kernel,H,c.RK,eta=3/5,m=2,ftregion=0,
                   vol.init=NA,covol.init=NA,
                   nvar.init=NA,ncov.init=NA,mn,alpha=0.4,
                   frequency=300,avg=TRUE,threshold,utime,psd=FALSE){
            
data <- get.zoo.data(x)
d.size <- length(data)

for(i in 1:d.size){
  
  # NA data must be skipped
  idt <- which(is.na(data[[i]]))
  if(length(idt>0)){
    data[[i]] <- data[[i]][-idt]
  }
  if(length(data[[i]])<2) {
    stop("length of data (w/o NA) must be more than 1")
  }
  
}

cmat <- NULL

switch(method,
       "HY"="<-"(cmat,HY(data)),
       "PHY"="<-"(cmat,PHY(data,theta,kn,g,refreshing,cwise)),
       "MRC"="<-"(cmat,MRC(data,theta,kn,g,delta,adj)),
       "TSCV"="<-"(cmat,TSCV(data,K,c.two,J,adj,utime)),  
       "GME"="<-"(cmat,GME(data,c.multi,utime)),
       "RK"="<-"(cmat,RK(data,kernel,H,c.RK,eta,m,ftregion,utime)),
       "QMLE"="<-"(cmat,cce.qmle(data,vol.init,covol.init,
                                 nvar.init,ncov.init)),
       "SIML"="<-"(cmat,SIML(data,mn,alpha)),
       "THY"="<-"(cmat,THY(data,threshold)),
       "PTHY"="<-"(cmat,PTHY(data,theta,kn,g,threshold,
                             refreshing,cwise)),
       "SRC"="<-"(cmat,SRC(data,frequency,avg,utime)),
       "SBPC"="<-"(cmat,SBPC(data,frequency,avg,utime)))

if(is.null(cmat))
  stop("method is not available")

if(psd){
  tmp <- svd(cmat%*%cmat)
  cmat <- tmp$u%*%diag(sqrt(tmp$d))%*%t(tmp$v)
}

if(d.size>1){
  if(all(diag(cmat)>0)){
    sdmat <- diag(sqrt(diag(cmat)))
    cormat <- solve(sdmat) %*% cmat %*% solve(sdmat)
  }else{
    cormat <- NA
  }
}else{
  cormat <- as.matrix(1)
}
rownames(cmat) <- names(data)
colnames(cmat) <- names(data)
rownames(cormat) <- names(data)
colnames(cormat) <- names(data)
return(list(covmat=cmat,cormat=cormat))
})
