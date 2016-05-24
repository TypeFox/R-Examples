#lead-lag estimation

#x:data

#setGeneric( "llag", function(x,verbose=FALSE) standardGeneric("llag") )
#setMethod( "llag", "yuima", function(x,verbose=FALSE) llag(x@data,verbose=verbose ))
#setMethod( "llag", "yuima.data", function(x,verbose=FALSE) {

#setGeneric( "llag",
#		function(x,from=FALSE,to=FALSE,division=FALSE,verbose=FALSE)
#		standardGeneric("llag") )
#setMethod( "llag", "yuima",
#		function(x,from=FALSE,to=FALSE,division=FALSE,verbose=FALSE)
#		llag(x@data,from=FALSE,to=FALSE,division=FALSE,verbose=verbose ))
#setMethod( "llag", "yuima.data", function(x,from=FALSE,to=FALSE,division=FALSE,verbose=FALSE) {


## function to make the grid if it is missing
make.grid <- function(d, d.size, x, from, to, division){
  
  out <- vector(d.size, mode = "list")
  
  if(length(from) != d.size){
    from <- c(from,rep(-Inf,d.size - length(from)))
  }
  
  if(length(to) != d.size){
    to <- c(to,rep(Inf,d.size - length(to)))
  }
  
  if(length(division) == 1){
    division <- rep(division,d.size)
  }
  
  for(i in 1:(d-1)){
    for(j in (i+1):d){
      
      time1 <- as.numeric(time(x[[i]]))
      time2 <- as.numeric(time(x[[j]]))
      
      # calculate the maximum of correlation by substituting theta to lagcce
      
      #n:=2*length(data)
      
      num <- d*(i-1) - (i-1)*i/2 + (j-i)
      
      if(division[num]==FALSE){
        n <- round(2*max(length(time1),length(time2)))+1
      }else{
        n <- division[num]
      }
      
      # maximum value of lagcce
      
      tmptheta <- time2[1]-time1[1] # time lag (sec)
      
      num1 <- time1[length(time1)]-time1[1] # elapsed time for series 1
      num2 <- time2[length(time2)]-time2[1] # elapsed time for series 2
      
      # modified
      
      if(is.numeric(from[num])==TRUE && is.numeric(to[num])==TRUE){
        num2 <- min(-from[num],num2+tmptheta)
        num1 <- min(to[num],num1-tmptheta)
        tmptheta <- 0
        
        if(-num2 >= num1){
          print("llag:invalid range")
          return(NULL)
        }
      }else if(is.numeric(from[num])==TRUE){
        num2 <- min(-from[num],num2+tmptheta)
        num1 <- num1-tmptheta
        tmptheta <- 0
        
        if(-num2 >= num1){
          print("llag:invalid range")
          return(NULL)
        }
      }else if(is.numeric(to[num])==TRUE){
        num2 <- num2+tmptheta
        num1 <- min(to[num],num1-tmptheta)
        tmptheta <- 0
        
        if(-num2 >= num1){
          print("llag:invalid range")
          return(NULL)
        }
      }
      
      out[[num]] <- seq(-num2-tmptheta,num1-tmptheta,length=n)[2:(n-1)]
      
    }
  }
  
  return(out)
}


## function to compute asymptotic variances
llag.avar <- function(x, grid, bw, alpha, fisher, ser.diffX, ser.times, vol, cormat, ccor, idx, G, d, d.size){
  
  # treatment of the bandwidth
  if(missing(bw)){
    
    bw <- matrix(0, d, d)
    
    for(i in 1:(d - 1)){
      
      for(j in (i + 1):d){
        
        Init <- min(ser.times[[i]][1], ser.times[[j]][1])
        Term <- max(tail(ser.times[[i]], n = 1), tail(ser.times[[j]], n = 1))
        bw[i, j] <- (Term - Init) * min(length(ser.times[[i]]), length(ser.times[[j]]))^(-0.45)
        bw[j, i] <- bw[i, j]
        
      }
      
    }
    
  }
  
  bw <- matrix(bw, d, d)
  
  p <- diag(d) # p-values
  avar <- vector(d.size,mode="list") # asymptotic variances
  CI <- vector(d.size,mode="list") # confidence intervals
  
  names(avar) <- names(ccor)
  names(CI) <- names(ccor)
  
  vv <- (2/3) * sapply(ser.diffX, FUN = function(x) sum(x^4)) # AVAR for RV
  
  for(i in 1:(d-1)){
    for(j in (i+1):d){
      
      time1 <- ser.times[[i]]
      time2 <- ser.times[[j]] 
      
      num <- d*(i-1) - (i-1)*i/2 + (j-i)
      
      ## computing conficence intervals
      N.max <- max(length(time1), length(time2))
      tmp <- rev(as.numeric(ccor[[num]]))
      
      d1 <- -0.5 * tmp * sqrt(vol[j]/vol[i])
      d2 <- -0.5 * tmp * sqrt(vol[i]/vol[j])
      
      avar.tmp <- .C("hycrossavar",
                     as.double(G[[num]]),
                     as.double(time1),
                     as.double(time2),
                     as.integer(length(G[[num]])),
                     as.integer(length(time1)),
                     as.integer(length(time2)),
                     as.double(x[[i]]),
                     as.double(x[[j]]),
                     as.integer(N.max),
                     as.double(bw[i, j]),
                     as.double(vv[i]),
                     as.double(vv[j]),
                     as.double(d1^2),
                     as.double(d2^2),
                     as.double(2 * d1),
                     as.double(2 * d2),
                     as.double(2 * d1 * d2),
                     covar = double(length(G[[num]])),
                     corr = double(length(G[[num]])))$corr / (vol[i] * vol[j])
      
      #avar[[num]][avar[[num]] <= 0] <- NA
      
      if(fisher == TRUE){
        z <- atanh(tmp) # the Fisher z transformation of the estimated correlation
        se.z <- sqrt(avar.tmp)/(1 - tmp^2)
        p[i,j] <- pchisq((atanh(cormat[i,j])/se.z[idx[num]])^2, df=1, lower.tail=FALSE)
        CI[[num]] <- zoo(tanh(qnorm(1 - alpha/2) * se.z), -grid[[num]])
      }else{
        p[i,j] <- pchisq(cormat[i,j]^2/avar.tmp[idx[num]], df=1, lower.tail=FALSE)
        c.alpha <- sqrt(qchisq(alpha, df=1, lower.tail = FALSE) * avar.tmp)
        CI[[num]] <- zoo(c.alpha, -grid[[num]])
      }
      
      p[j,i] <- p[i,j]
      avar[[num]] <- zoo(avar.tmp, -grid[[num]])
      
    }
  }
  
  return(list(p = p, avar = avar, CI = CI))
}


## main body
setGeneric( "llag", function(x, from = -Inf, to = Inf, division = FALSE, 
                             verbose = (ci || ccor), grid, psd = TRUE, plot = ci,
                             ccor = ci, ci = FALSE, alpha = 0.01, fisher = TRUE, bw) standardGeneric("llag") )

## yuima-method
setMethod("llag", "yuima", function(x, from, to, division, verbose, grid, psd, plot, 
                                         ccor, ci, alpha, fisher, bw)
  llag(x@data, from, to, division, verbose, grid, psd, plot, ccor, ci, alpha, fisher, bw))

## yuima.data-method
setMethod("llag", "yuima.data", function(x, from, to, division, verbose, grid, psd, plot, 
                                         ccor, ci, alpha, fisher, bw)
  llag(x@zoo.data, from, to, division, verbose, grid, psd, plot, ccor, ci, alpha, fisher, bw))

## list-method
setMethod("llag", "list", function(x, from, to, division, verbose, grid, psd, plot, 
                                   ccor, ci, alpha, fisher, bw) {
  
  d <- length(x)
  
  # allocate memory
  ser.times <- vector(d, mode="list") # time index in 'x'
  ser.diffX <- vector(d, mode="list") # difference of data
  vol <- double(d)
  
  # Set the tolerance to avoid numerical erros in comparison
  tol <- 1e-6
  
  for(i in 1:d){
    
    # NA data must be skipped
    idt <- which(is.na(x[[i]]))
    if(length(idt>0)){
      x[[i]] <- x[[i]][-idt]
    }
    if(length(x[[i]])<2) {
      stop("length of data (w/o NA) must be more than 1")
    }
    
    # set data and time index
    ser.times[[i]] <- as.numeric(time(x[[i]]))/tol
    # set difference of the data 
    ser.diffX[[i]] <- diff( as.numeric(x[[i]]) )
    vol[i] <- sum(ser.diffX[[i]]^2)
  }
  
  theta <- matrix(0,d,d)
  
  d.size <- d*(d-1)/2
  crosscor <- vector(d.size,mode="list")
  idx <- integer(d.size)
  
  # treatment of the grid
  if(missing(grid)) 
    grid <- make.grid(d, d.size, x, from, to, division)
  
  if(is.list(grid)){
    G <- relist(unlist(grid)/tol, grid)
  }else{
    if(is.numeric(grid)){
      G <- data.frame(matrix(grid/tol,length(grid),d.size))
      grid <- data.frame(matrix(grid,length(grid),d.size))
    }else{
      print("llag:invalid grid")
      return(NULL)
    }
  }
  
  # core part
  if(psd){ # positive semidefinite correction is implemented
    
    cormat <- diag(d)
    LLR <- diag(d)
    
    for(i in 1:(d-1)){
      for(j in (i+1):d){
        
        time1 <- ser.times[[i]]
        time2 <- ser.times[[j]] 
        
        num <- d*(i-1) - (i-1)*i/2 + (j-i)
        
        names(crosscor)[num] <- paste("(",i,",",j,")", sep = "")
        
        tmp <- .C("HYcrosscorr",
                  as.integer(length(G[[num]])),
                  as.integer(length(time1)),
                  as.integer(length(time2)),
                  as.double(G[[num]]),
                  as.double(time1),
                  as.double(time2),
                  double(length(time2)),
                  as.double(ser.diffX[[i]]),
                  as.double(ser.diffX[[j]]),
                  as.double(vol[i]),
                  as.double(vol[j]),
                  value=double(length(G[[num]])))$value
        
        idx[num] <- which.max(abs(tmp))
        mlag <- -grid[[num]][idx[num]] # make the first timing of max or min
        cor <- tmp[idx[num]]
        
        theta[i,j] <- mlag
        cormat[i,j] <- cor
        theta[j,i] <- -mlag
        cormat[j,i] <- cormat[i,j]
        
        LLR[i,j] <- sum(tmp[grid[[num]] < 0]^2)/sum(tmp[grid[[num]] > 0]^2)
        LLR[j,i] <- 1/LLR[i,j]
        
        crosscor[[num]] <- zoo(tmp,-grid[[num]])
      }
    }
    
    covmat <- diag(sqrt(vol))%*%cormat%*%diag(sqrt(vol))
    
  }else{# non-psd
    
    covmat <- diag(vol)
    LLR <- diag(d)
    
    for(i in 1:(d-1)){
      for(j in (i+1):d){
        
        time1 <- ser.times[[i]]
        time2 <- ser.times[[j]] 
        
        num <- d*(i-1) - (i-1)*i/2 + (j-i)
        
        names(crosscor)[num] <- paste("(",i,",",j,")", sep = "")
        
        tmp <- .C("HYcrosscov",
                  as.integer(length(G[[num]])),
                  as.integer(length(time1)),
                  as.integer(length(time2)),
                  as.double(G[[num]]),
                  as.double(time1),
                  as.double(time2),
                  double(length(time2)),
                  as.double(ser.diffX[[i]]),
                  as.double(ser.diffX[[j]]),
                  value=double(length(G[[num]])))$value
        
        idx[num] <- which.max(abs(tmp))
        mlag <- -grid[[num]][idx[num]] # make the first timing of max or min
        cov <- tmp[idx[num]]
        
        theta[i,j] <- mlag
        covmat[i,j] <- cov
        theta[j,i] <- -mlag
        covmat[j,i] <- covmat[i,j]
        
        LLR[i,j] <- sum(tmp[grid[[num]] < 0]^2)/sum(tmp[grid[[num]] > 0]^2)
        LLR[j,i] <- 1/LLR[i,j]
        
        crosscor[[num]] <- zoo(tmp,-grid[[num]])/sqrt(vol[i]*vol[j])
      }
    }
    
    cormat <- diag(1/sqrt(diag(covmat)))%*%covmat%*%diag(1/sqrt(diag(covmat)))
    
  }
  
  if(ci){ # computing confidence intervals
    
    out <- llag.avar(x = x, grid = grid, bw = bw, alpha = alpha, fisher = fisher,
                     ser.diffX = ser.diffX, ser.times = ser.times, 
                     vol = vol, cormat = cormat, ccor = crosscor, idx = idx, 
                     G = G, d = d, d.size = d.size)
    
    p <- out$p
    CI <- out$CI
    avar <- out$avar
    
    colnames(theta) <- names(x)
    rownames(theta) <- names(x)
    colnames(p) <- names(x)
    rownames(p) <- names(x)
    
    if(plot){
      for(i in 1:(d-1)){
        for(j in (i+1):d){
          
          num <- d*(i-1) - (i-1)*i/2 + (j-i)
          y.max <- max(abs(as.numeric(crosscor[[num]])),as.numeric(CI[[num]]))
          
          plot(crosscor[[num]],
               main=paste(i,"vs",j,"(positive",expression(theta),"means",i,"leads",j,")"),
               xlab=expression(theta),ylab=expression(U(theta)),
               ylim=c(-y.max,y.max))
          
          lines(CI[[num]],lty=2,col="blue")
          lines(-CI[[num]],lty=2,col="blue")
        }
      }
    }
    
    if(verbose==TRUE){
      
      colnames(covmat) <- names(x)
      rownames(covmat) <- names(x)
      colnames(cormat) <- names(x)
      rownames(cormat) <- names(x)
      colnames(LLR) <- names(x)
      rownames(LLR) <- names(x)
      
      if(ccor){
        result <- list(lagcce=theta,p.values=p,covmat=covmat,cormat=cormat,LLR=LLR,ccor=crosscor,avar=avar)
      }else{
        result <- list(lagcce=theta,p.values=p,covmat=covmat,cormat=cormat,LLR=LLR)
      }
      
    }else{
      return(theta)
    }
    
  }else{
    
    colnames(theta) <- names(x)
    rownames(theta) <- names(x)
    
    if(plot){
      for(i in 1:(d-1)){
        for(j in (i+1):d){
          num <- d*(i-1) - (i-1)*i/2 + (j-i)
          plot(crosscor[[num]],
               main=paste(i,"vs",j,"(positive",expression(theta),"means",i,"leads",j,")"),
               xlab=expression(theta),ylab=expression(U(theta)))
        }
      }
    }
    
    if(verbose==TRUE){
      
      colnames(covmat) <- names(x)
      rownames(covmat) <- names(x)
      colnames(cormat) <- names(x)
      rownames(cormat) <- names(x)
      colnames(LLR) <- names(x)
      rownames(LLR) <- names(x)
      
      if(ccor){
        result <- list(lagcce=theta,covmat=covmat,cormat=cormat,LLR=LLR,ccor=crosscor)
      }else{
        result <- list(lagcce=theta,covmat=covmat,cormat=cormat,LLR=LLR)
      }
    }else{
      return(theta)
    }
    
  }
  
  class(result) <- "yuima.llag"
  
  return(result)
})

# print method for yuima.llag-class
print.yuima.llag <- function(x, ...){
  
  if(is.null(x$p.values)){
    
    cat("Estimated lead-lag parameters\n")
    print(x$lagcce, ...)
    cat("Corresponding covariance matrix\n")
    print(x$covmat, ...)
    cat("Corresponding correlation matrix\n")
    print(x$cormat, ...)
    cat("Lead-lag ratio\n")
    print(x$LLR, ...)
    
  }else{
    
    cat("Estimated lead-lag parameters\n")
    print(x$lagcce, ...)
    cat("Corresponding p-values\n")
    print(x$p.values, ...)
    cat("Corresponding covariance matrix\n")
    print(x$covmat, ...)
    cat("Corresponding correlation matrix\n")
    print(x$cormat, ...)
    cat("Lead-lag ratio\n")
    print(x$LLR, ...)
    
  }
  
}

# plot method for yuima.llag-class
plot.yuima.llag <- function(x, alpha = 0.01, fisher = TRUE, ...){
  
  if(is.null(x$ccor)){
    warning("cross-correlation functions were not returned by llag. Set verbose = TRUE and ccor = TRUE to return them.",
            call. = FALSE)
    return(NULL)
  }else{
    
    d <- nrow(x$LLR)
    
    if(is.null(x$avar)){
      for(i in 1:(d-1)){
        for(j in (i+1):d){
          
          num <- d*(i-1) - (i-1)*i/2 + (j-i)
          
          plot(x$ccor[[num]],
               main=paste(i,"vs",j,"(positive",expression(theta),"means",i,"leads",j,")"),
               xlab=expression(theta),ylab=expression(U(theta)))
          
        }
      }
      
    }else{
      for(i in 1:(d-1)){
        for(j in (i+1):d){
          
          num <- d*(i-1) - (i-1)*i/2 + (j-i)
          G <- time(x$ccor[[num]])
          
          if(fisher == TRUE){
            z <- atanh(x$ccor[[num]]) # the Fisher z transformation of the estimated correlation
            se.z <- sqrt(x$avar[[num]])/(1 - x$ccor[[num]]^2)
            CI <- zoo(tanh(qnorm(1 - alpha/2) * se.z), G)
          }else{
            c.alpha <- sqrt(qchisq(alpha, df=1, lower.tail = FALSE) * x$avar[[num]])
            CI <- zoo(c.alpha, G)
          }
          
          y.max <- max(abs(as.numeric(x$ccor[[num]])),as.numeric(CI))
          
          plot(x$ccor[[num]],
               main=paste(i,"vs",j,"(positive",expression(theta),"means",i,"leads",j,")"),
               xlab=expression(theta),ylab=expression(U(theta)),
               ylim=c(-y.max,y.max))
          
          lines(CI,lty=2,col="blue")
          lines(-CI,lty=2,col="blue")
          
        }
      }
    }
    
  }
  
}


## Old version until Oct. 10, 2015
if(0){
setGeneric( "llag",
		function(x,from=-Inf,to=Inf,division=FALSE,verbose=FALSE,grid,psd=TRUE,
             plot=FALSE,ccor=FALSE)
		standardGeneric("llag") )
setMethod( "llag", "yuima",
		function(x,from=-Inf,to=Inf,division=FALSE,verbose=FALSE,grid,psd=TRUE,
             plot=FALSE,ccor=FALSE)
		llag(x@data,from=from,to=to,division=division,verbose=verbose,grid=grid,
         psd=psd,plot=plot))
setMethod( "llag", "yuima.data", function(x,from=-Inf,to=Inf,division=FALSE,
                                          verbose=FALSE,grid,psd=TRUE,
                                          plot=FALSE,ccor=FALSE) {
  
  if((is(x)=="yuima")||(is(x)=="yuima.data")){
    zdata <- get.zoo.data(x)
  }else{
    print("llag:invalid argument")
    return(NULL)
  }
  
  d <- length(zdata)
  
  # allocate memory
  ser.times <- vector(d, mode="list") # time index in 'x'
  ser.diffX <- vector(d, mode="list") # difference of data
  vol <- double(d)
  
  # Set the tolerance to avoid numerical erros in comparison
  tol <- 1e-6
  
  for(i in 1:d){
    
    # NA data must be skipped
    idt <- which(is.na(zdata[[i]]))
    if(length(idt>0)){
      zdata[[i]] <- zdata[[i]][-idt]
    }
    if(length(zdata[[i]])<2) {
      stop("length of data (w/o NA) must be more than 1")
    }
    
    # set data and time index
    ser.times[[i]] <- as.numeric(time(zdata[[i]]))
    # set difference of the data 
    ser.diffX[[i]] <- diff( as.numeric(zdata[[i]]) )
    vol[i] <- sum(ser.diffX[[i]]^2)
  }
  
  theta <- matrix(0,d,d)
  
  d.size <- d*(d-1)/2
  crosscor <- vector(d.size,mode="list")
  
  if(psd){
    
    cormat <- diag(d)
    
    if(missing(grid)){
      
      if(length(from) != d.size){
        from <- c(from,rep(-Inf,d.size - length(from)))
      }
      
      if(length(to) != d.size){
        to <- c(to,rep(Inf,d.size - length(to)))
      }
      
      if(length(division) == 1){
        division <- rep(division,d.size)
      }
      
      for(i in 1:(d-1)){
        for(j in (i+1):d){
          
          time1 <- ser.times[[i]]
          time2 <- ser.times[[j]] 
          
          # calculate the maximum of correlation by substituting theta to lagcce
          
          #n:=2*length(data)
          
          num <- d*(i-1) - (i-1)*i/2 + (j-i)
          
          if(division[num]==FALSE){
            n <- round(2*max(length(time1),length(time2)))+1
          }else{
            n <- division[num]
          }
          
          # maximum value of lagcce
          
          tmptheta <- time2[1]-time1[1] # time lag (sec)
          
          num1 <- time1[length(time1)]-time1[1] # elapsed time for series 1
          num2 <- time2[length(time2)]-time2[1] # elapsed time for series 2
          
          # modified
          
          if(is.numeric(from[num])==TRUE && is.numeric(to[num])==TRUE){
            num2 <- min(-from[num],num2+tmptheta)
            num1 <- min(to[num],num1-tmptheta)
            tmptheta <- 0
            
            if(-num2 >= num1){
              print("llag:invalid range")
              return(NULL)
            }
          }else if(is.numeric(from[num])==TRUE){
            num2 <- min(-from[num],num2+tmptheta)
            num1 <- num1-tmptheta
            tmptheta <- 0
            
            if(-num2 >= num1){
              print("llag:invalid range")
              return(NULL)
            }
          }else if(is.numeric(to[num])==TRUE){
            num2 <- num2+tmptheta
            num1 <- min(to[num],num1-tmptheta)
            tmptheta <- 0
            
            if(-num2 >= num1){
              print("llag:invalid range")
              return(NULL)
            }
          }
          
          y <- seq(-num2-tmptheta,num1-tmptheta,length=n)[2:(n-1)]
          
          tmp <- .C("HYcrosscorr",
                    as.integer(n-2),
                    as.integer(length(time1)),
                    as.integer(length(time2)),
                    as.double(y/tol),
                    as.double(time1/tol),
                    as.double(time2/tol),
                    double(length(time2)),
                    as.double(ser.diffX[[i]]),
                    as.double(ser.diffX[[j]]),
                    as.double(vol[i]),
                    as.double(vol[j]),
                    value=double(n-2))$value
          
          idx <- which.max(abs(tmp))
          mlag <- -y[idx] # make the first timing of max or min
          corr <- tmp[idx]
          
          theta[i,j] <- mlag
          cormat[i,j] <- corr
          theta[j,i] <- -mlag
          cormat[j,i] <- cormat[i,j]
          
          crosscor[[num]] <- zoo(tmp,-y)
        }
      }
      
    }else{
      
      if(!is.list(grid)){
        if(is.numeric(grid)){
          grid <- data.frame(matrix(grid,length(grid),d.size))
        }else{
          print("llag:invalid grid")
          return(NULL)
        }
      }
      
      for(i in 1:(d-1)){
        for(j in (i+1):d){
          
          time1 <- ser.times[[i]]
          time2 <- ser.times[[j]] 
          
          num <- d*(i-1) - (i-1)*i/2 + (j-i)
          
          tmp <- .C("HYcrosscorr",
                    as.integer(length(grid[[num]])),
                    as.integer(length(time1)),
                    as.integer(length(time2)),
                    as.double(grid[[num]]/tol),
                    as.double(time1/tol),
                    as.double(time2/tol),
                    double(length(time2)),
                    as.double(ser.diffX[[i]]),
                    as.double(ser.diffX[[j]]),
                    as.double(vol[i]),
                    as.double(vol[j]),
                    value=double(length(grid[[num]])))$value
          
          idx <- which.max(abs(tmp))
          mlag <- -grid[[num]][idx] # make the first timing of max or min
          cor <- tmp[idx]
          
          theta[i,j] <- mlag
          cormat[i,j] <- cor
          theta[j,i] <- -mlag
          cormat[j,i] <- cormat[i,j]
          
          crosscor[[num]] <- zoo(tmp,-grid[[num]])
        }
      }
    }
    
    covmat <- diag(sqrt(vol))%*%cormat%*%diag(sqrt(vol))
    
  }else{# non-psd
    
    covmat <- diag(vol)
    
    if(missing(grid)){
      
      if(length(from) != d.size){
        from <- c(from,rep(-Inf,d.size - length(from)))
      }
      
      if(length(to) != d.size){
        to <- c(to,rep(Inf,d.size - length(to)))
      }
      
      if(length(division) == 1){
        division <- rep(division,d.size)
      }
      
      for(i in 1:(d-1)){
        for(j in (i+1):d){
          
          time1 <- ser.times[[i]]
          time2 <- ser.times[[j]] 
          
          # calculate the maximum of correlation by substituting theta to lagcce
          
          #n:=2*length(data)
          
          num <- d*(i-1) - (i-1)*i/2 + (j-i)
          
          if(division[num]==FALSE){
            n <- round(2*max(length(time1),length(time2)))+1
          }else{
            n <- division[num]
          }
          
          # maximum value of lagcce
          
          tmptheta <- time2[1]-time1[1] # time lag (sec)
          
          num1 <- time1[length(time1)]-time1[1] # elapsed time for series 1
          num2 <- time2[length(time2)]-time2[1] # elapsed time for series 2
          
          # modified
          
          if(is.numeric(from[num])==TRUE && is.numeric(to[num])==TRUE){
            num2 <- min(-from[num],num2+tmptheta)
            num1 <- min(to[num],num1-tmptheta)
            tmptheta <- 0
            
            if(-num2 >= num1){
              print("llag:invalid range")
              return(NULL)
            }
          }else if(is.numeric(from[num])==TRUE){
            num2 <- min(-from[num],num2+tmptheta)
            num1 <- num1-tmptheta
            tmptheta <- 0
            
            if(-num2 >= num1){
              print("llag:invalid range")
              return(NULL)
            }
          }else if(is.numeric(to[num])==TRUE){
            num2 <- num2+tmptheta
            num1 <- min(to[num],num1-tmptheta)
            tmptheta <- 0
            
            if(-num2 >= num1){
              print("llag:invalid range")
              return(NULL)
            }
          }
          
          y <- seq(-num2-tmptheta,num1-tmptheta,length=n)[2:(n-1)]
          
          tmp <- .C("HYcrosscov",
                    as.integer(n-2),
                    as.integer(length(time1)),
                    as.integer(length(time2)),
                    as.double(y/tol),
                    as.double(time1/tol),
                    as.double(time2/tol),
                    double(length(time2)),
                    as.double(ser.diffX[[i]]),
                    as.double(ser.diffX[[j]]),
                    value=double(n-2))$value
          
          idx <- which.max(abs(tmp))
          mlag <- -y[idx] # make the first timing of max or min
          cov <- tmp[idx]
          
          theta[i,j] <- mlag
          covmat[i,j] <- cov
          theta[j,i] <- -mlag
          covmat[j,i] <- covmat[i,j]
          
          crosscor[[num]] <- zoo(tmp,-y)/sqrt(vol[i]*vol[j])
        }
      }
      
    }else{
      
      if(!is.list(grid)){
        if(is.numeric(grid)){
          grid <- data.frame(matrix(grid,length(grid),d.size))
        }else{
          print("llag:invalid grid")
          return(NULL)
        }
      }
      
      for(i in 1:(d-1)){
        for(j in (i+1):d){
          
          time1 <- ser.times[[i]]
          time2 <- ser.times[[j]] 
          
          num <- d*(i-1) - (i-1)*i/2 + (j-i)
          
          tmp <- .C("HYcrosscov",
                    as.integer(length(grid[[num]])),
                    as.integer(length(time1)),
                    as.integer(length(time2)),
                    as.double(grid[[num]]/tol),
                    as.double(time1/tol),
                    as.double(time2/tol),
                    double(length(time2)),
                    as.double(ser.diffX[[i]]),
                    as.double(ser.diffX[[j]]),
                    value=double(length(grid[[num]])))$value
          
          idx <- which.max(abs(tmp))
          mlag <- -grid[[num]][idx] # make the first timing of max or min
          cov <- tmp[idx]
          
          theta[i,j] <- mlag
          covmat[i,j] <- cov
          theta[j,i] <- -mlag
          covmat[j,i] <- covmat[i,j]
          
          crosscor[[num]] <- zoo(tmp,-grid[[num]])/sqrt(vol[i]*vol[j])
        }
      }
    }
    
    cormat <- diag(1/sqrt(diag(covmat)))%*%covmat%*%diag(1/sqrt(diag(covmat)))
  }
  
  colnames(theta) <- names(zdata)
  rownames(theta) <- names(zdata)
  
  if(plot){
    for(i in 1:(d-1)){
      for(j in (i+1):d){
        num <- d*(i-1) - (i-1)*i/2 + (j-i)
        plot(crosscor[[num]],
             main=paste(i,"vs",j,"(positive",expression(theta),"means",i,"leads",j,")"),
             xlab=expression(theta),ylab=expression(U(theta)))
      }
    }
  }
  
  if(verbose==TRUE){
    colnames(covmat) <- names(zdata)
    rownames(covmat) <- names(zdata)
    colnames(cormat) <- names(zdata)
    rownames(cormat) <- names(zdata)
    if(ccor){
      return(list(lagcce=theta,covmat=covmat,cormat=cormat,ccor=crosscor))
    }else{
      return(list(lagcce=theta,covmat=covmat,cormat=cormat))
    }
  }else{
    return(theta)
  }
})
}

## Old version
if(0){
setMethod( "llag", "yuima.data", function(x,from=-Inf,to=Inf,division=FALSE,verbose=FALSE) {
  
  
  if(!is(x)=="yuima.data"){
    if(is(x)=="yuima"){
      dat <- x@data
    }else{
      print("llag:invalid argument")
      return(NULL)
    }
  }else{
    dat <- x
  }
  
  d <- length(dat@zoo.data)
  
  lagccep <- function(datp,theta){
    time(datp[[2]]) <- time(datp[[2]])+theta
    return(cce(setData(datp))$covmat[1,2])
  }
  
  lagcce <- function(datzoo,theta){
    d <- dim(theta)[1]
    lcmat <- cce(setData(datzoo))$covmat
    
    for(i in 1:(d-1)){
      for(j in (i+1):d){
        datp <- datzoo[c(i,j)]
        lcmat[i,j] <- lagccep(datp,theta[i,j])
        lcmat[j,i] <- lcmat[i,j]
      }
    }		
    return(lcmat)	
  }
  
  d.size <- d*(d-1)/2
  
  if(length(from) != d.size){
    from <- c(from,rep(-Inf,d.size - length(from)))
  }
  
  if(length(to) != d.size){
    to <- c(to,rep(Inf,d.size - length(to)))
  }
  
  if(length(division) != d.size){
    division <- c(division,rep(FALSE,d.size - length(division)))
  }
  
  find_lag <- function(i,j){
    datp <- dat@zoo.data[c(i,j)]
    time1 <- time(datp[[1]])
    time2 <- time(datp[[2]])  
    
    # calculate the maximum of correlation by substituting theta to lagcce
    
    #n:=2*length(data)
    
    num <- d*(i-1) - (i-1)*i/2 + (j-i)
    
    if(division[num]==FALSE){
      n <- round(2*max(length(datp[[1]]),length(datp[[2]])))+1
    }else{
      n <- division[num]
    }
    
    # maximum value of lagcce
    
    tmptheta <- as.numeric(time2[1])-as.numeric(time1[1]) # time lag (sec)
    
    num1 <- as.numeric(time1[length(time1)])-as.numeric(time1[1]) # elapsed time for series 1
    num2 <- as.numeric(time2[length(time2)])-as.numeric(time2[1]) # elapsed time for series 2
    
    # modified
    
    if(is.numeric(from[num])==TRUE && is.numeric(to[num])==TRUE){
      num2 <- min(-from[num],num2+tmptheta)
      num1 <- min(to[num],num1-tmptheta)
      tmptheta <- 0
      
      if(-num2 >= num1){
        print("llag:invalid range")
        return(NULL)
      }
    }else if(is.numeric(from[num])==TRUE){
      num2 <- min(-from[num],num2+tmptheta)
      num1 <- num1-tmptheta
      tmptheta <- 0
      
      if(-num2 >= num1){
        print("llag:invalid range")
        return(NULL)
      }
    }else if(is.numeric(to[num])==TRUE){
      num2 <- num2+tmptheta
      num1 <- min(to[num],num1-tmptheta)
      tmptheta <- 0
      
      if(-num2 >= num1){
        print("llag:invalid range")
        return(NULL)
      }
    }
    
    y <- seq(-num2-tmptheta,num1-tmptheta,length=n)
    tmp <- double(n)
    
    for(i.tmp in 2:(n-1)){
      tmp[i.tmp] <- lagccep(datp,y[i.tmp])
    }
    
    mat <- cbind(y[2:(n-1)],tmp[2:(n-1)])
    
    #idx <- abs(mat[,2])==max(abs(mat[,2]))
    #mlag <- mat[,1][idx][1] # make the first timing of max or min
    #cov <- mat[,2][idx][1]
    idx <- which.max(abs(mat[,2]))
    mlag <- mat[,1][idx] # make the first timing of max or min
    cov <- mat[,2][idx]
    return(list(lag=-mlag,cov=cov))
  }
  
  theta <- matrix(numeric(d^2),ncol=d)
  covmat <- cce(dat)$covmat
  
  for(i in 1:(d-1)){
    for(j in (i+1):d){
      fl <- find_lag(i,j)
      theta[i,j] <- fl$lag
      covmat[i,j] <- fl$cov
      theta[j,i] <- -theta[i,j]
      covmat[j,i] <- covmat[i,j]
    }
  }
  
  
  #  covmat <- lagcce(dat@zoo.data,theta)
  cormat <- diag(1/sqrt(diag(covmat)))%*%covmat%*%diag(1/sqrt(diag(covmat)))
  if(verbose==TRUE){
    return(list(lagcce=theta,covmat=covmat,cormat=cormat))
  }else{
    return(theta)
  }
})
}