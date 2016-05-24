
gjamSimData <-
function( n=1000,S=10,q=5,nmiss=0,typeNames, effort = NULL ){
  
  # nmiss  = number of missing values in x
 # effort = used for DA data, format: 
  #          list(columns = 1:S, values = rep(1,n))
  
  pg <- .95
  
  if(length(typeNames) == 1)typeNames <- rep(typeNames,S)
  Q  <- q
  S1 <- S
  
  cuts <- numeric(0)
  
  tmp <- .gjamGetTypes(typeNames)
  typeCols <- tmp$typeCols
  typeFull <- tmp$typeFull
  typeCode <- tmp$TYPES[typeCols]
  allTypes <- sort(unique(typeCols))
  
  x <- matrix( rnorm(n*q,.1), n,q)  
  x[,1] <- 1
  
  beta <- matrix(0,Q,S)
  ss   <- diag(.01,S)    
  
  notOther <- c(1:S)
  
  for(k in allTypes){
    
    wk <- which(typeCols == k)
    nk <- length(wk)
    
    if(typeFull[wk[1]] == 'presenceAbsence'){
      diag(ss)[wk] <- 1
      beta[,wk]    <- runif(q*nk,-1.5,1.5)
    }
    if(typeFull[wk[1]] == 'continuous'){
      diag(ss)[wk] <- .4
      beta[,wk]    <- runif(q*nk,0,1.5)
    }
    if(typeFull[wk[1]] == 'discrete'){
      diag(ss)[wk] <- 1
      beta[,wk]    <- runif(q*nk,-.1,2)
    }
    if(typeFull[wk[1]] %in% c('fracComp','countComp')){
      
      if(length(wk) < 2)stop('composition data must have at least 2 columns')
      
      ym <- matrix(runif(n*nk,.4/nk,1/nk),n,nk)
      bb <- solve(crossprod(x[,-1]))%*%t(x[,-1])%*%ym
      beta[,wk]    <- rbind(runif(nk,0,1/nk),bb)
      diag(ss)[wk] <- (1/S/q)^2
    }
    if(typeFull[wk[1]] == 'ordinal'){
      diag(ss)[wk] <- 1
      beta[,wk]    <- runif(q*nk,-.4,2)
    }
  }

  xnames <- paste('x',1:q,sep='')
  snames <- cnames <- paste('S',1:S,sep='')   
  dnames <- paste(snames,typeCode,sep='~')
  
  sigma <- cov( .rMVN(S+10,0,ss) )       
  w     <- x%*%beta + .rMVN(n,0,sigma) 
  colnames(w) <- snames
  
  y  <- w
  y[y < 0] <- 0
  z  <- w*0
  z[w <= 0]   <- 1
  z[w > 0]    <- 2
  plo <- phi  <- w*0
  plo[z == 1] <- -Inf
  phi[z == 2] <- Inf
  
  for(k in allTypes){
    
    wk <- which(typeCols == k)
    nk <- length(wk)
    
    if( typeFull[wk[1]] == 'presenceAbsence' )y[,wk] <- z[,wk] - 1       

    if( typeFull[wk[1]] == 'discrete' ){
      
      if(!is.null(effort)){
        we     <- wk[wk %in% effort$columns]
        y[,we] <- round( w[,we]*effort$values,0 )
        wide   <- .5/effort$values
      } else {
        wide   <- .5
        y[,wk] <- round(w[,wk],0)
      }
      y[,wk][y[,wk] < 0] <- 0
      
      hi <- w[,wk] + wide
      lo <- w[,wk] - wide
      
      lo[w[,wk] < 0] <- -Inf
      plo[,wk]   <- lo
      phi[,wk]   <- hi
      hi[hi < 0] <- 0
      z[,wk]     <- y[,wk] + 1
    }
    
    if( typeFull[wk[1]] %in% c('fracComp','countComp' )){
      
      snames[[wk[nk]]] <- 'other'
      noto <- 1:(nk - 1)
      
      ww     <- w[,wk]
      w0     <- which(ww < 0)
      w1     <- which(ww >= 0)
      ww[w0] <- 0  
      
      yk  <- .gjamCompW2Y(ww,pg,notOther=noto)$ww

      lo <- plo[,wk]
      lo[lo == -Inf] <- -1
      hi <- phi[,wk]
      hi[hi == Inf] <- 1
      plo[,wk] <- lo
      phi[,wk] <- hi
      
      beta[,wk[nk]] <- 0 
      
      if(typeFull[wk[1]] == 'countComp'){
        
        ee <- rpois(n,rgamma(n,100))
        yy <- sweep(yk,1,ee,'*')
        
        ww <- ceiling(yy)
        ww[ww < 0] <- 0
        
        lo <- (ww - 1)/ee 
        hi <- ww /ee
        
        ll <- lo
        ll[ll < 0] <- 0
        ll <- matrix(rowSums(ll),n,nk)
        
        hi[ww == 0] <- 1 - ll[ww == 0]
        lo[lo < 0]  <- -1
        
        plo[,wk] <- lo
        phi[,wk] <- hi
        
        y[,wk] <- ww
        z[,wk] <- ww + 1
      }
    }
    
    if( typeFull[wk[1]] == 'ordinal' ){
      
      yy   <- w[,wk]
      ncut <- 8
      maxw <- floor(max(yy))
      
      cuts  <- t( matrix( c(-Inf, seq(0,(maxw-1),length=(ncut-2)) ,Inf),ncut,nk) )
      rownames(cuts) <- snames[wk]
      
      for(j in 1:nk){
        z[,wk[j]]   <- findInterval(yy[,j],cuts[j,])
        plo[,wk[j]] <- cuts[j,z[,wk[j]]]
        phi[,wk[j]] <- cuts[j,1 + z[,wk[j]]]
      }
      
      y[,wk] <- z[,wk] - 1
      
      cuts <- .gjamTheta2cuts(cuts,sigma[wk,wk])   
    }
  }
  
  if(nmiss > 0){
    x[ sample(length(x),nmiss) ] <- NA
    x[,1] <- 1
    wmiss <- which(is.na(x),arr.ind=T)
    nmiss <- nrow(wmiss)
  }
  
  xnames[1]      <- 'intercept'
  colnames(y)    <- snames
  colnames(beta) <- rownames(sigma) <- colnames(sigma) <- snames
  colnames(x)    <- rownames(beta) <- xnames
  
  form <- as.formula( paste('~ ',paste(colnames(x)[-1],collapse='+' )) )
  
  list(formula = form, xdata = data.frame(x), y = y, w = w,  
       typeNames = typeNames, effort = effort,
       trueValues = list(beta = beta, sigma = sigma, 
                         corSpec = .cov2Cor(sigma), cuts = cuts))
}
