
computeRealVolumeConstant <- function(model,integration.points,integration.weights=NULL,T){
	#a modifier pour pouvoir travailler sur 2n pts d'integration
  
  if(is.null(integration.weights)){
      case <- 1
  }else{
    if(nrow(integration.points) == length(integration.weights)){
      case <- 1
    }else{
      case <- 2
    }
  }
  
  if(case==1){
  	obj <- predict_nobias_km(object=model,newdata=integration.points,type="UK",se.compute=TRUE,cov.compute=TRUE)
  	sigma.M <- obj$cov
  	intpoints.oldmean <- obj$mean
  	intpoints.oldsd <- obj$sd 
  	
  	sigma.M.tab <- as.numeric(sigma.M)
  	M <- nrow(integration.points)
  	indices <- expand.grid(c(1:M),c(1:M))
  	col1 <- indices[,1]
  	col2 <- indices[,2]
  	
  	the.Phi <- pnorm((T - intpoints.oldmean)/intpoints.oldsd)
  	Phi.a <- the.Phi[col1]
  	Phi.b <- the.Phi[col2]
  	
  	a <- (T - intpoints.oldmean[col1])/intpoints.oldsd[col1]
  	b <- (T - intpoints.oldmean[col2])/intpoints.oldsd[col2]
	
	a[a==Inf]<- 1000 ;a[a== -Inf] <- 1000;b[b==Inf]<-0;b[b== -Inf]<-0
	a[is.nan(a)] <- 1000; b[is.nan(b)]<-0
	
  	correl <- pmin(1,sigma.M.tab /(intpoints.oldsd[col1] * intpoints.oldsd[col2]))
	correl <- pmin(correl,1);correl <- pmax(correl,-1)
	
  	Phi.ab <- pbivnorm(a,b,correl)
  	constant.tab <- 1 - Phi.a - Phi.b + Phi.ab
  	
  	if (is.null(integration.weights)) {
  		constant.result <- mean(constant.tab)
  	}else{ 
  		new.weights <- integration.weights[col1]*integration.weights[col2]
  		constant.result <- sum(constant.tab*new.weights)
  	}
  	return(constant.result)
  }else{
    #sample generated with jn importance sampling
    #implementation brutale
    n <- nrow(integration.points)/2
    xtab <- integration.points[c(1:n),]
    ytab <- integration.points[c((n+1):(2*n)),]
    
    if(model@d == 1) xtab <- matrix(xtab,ncol=1)
    if(model@d == 1) ytab <- matrix(ytab,ncol=1)
    
    batch <- 100
    n.batch <- ceiling(n / batch)
    
    kn <- NULL
    mn1 <- NULL
    sn1 <- NULL
    mn2 <- NULL
    sn2 <- NULL
    
    first <- 1
    last <- min(n,batch)
    size <- 1 + last - first
    indices <- c(first:last)
    
    for (i in 1:n.batch){
      #build the matrix of integration points where we compute the prediction
      tmp1 <- xtab[indices,]
      tmp2 <- ytab[indices,]
      if(model@d == 1) tmp1 <- matrix(tmp1,ncol=1)
      if(model@d == 1) tmp2 <- matrix(tmp2,ncol=1)
      
      mat <- rbind(tmp1,tmp2)
      obj <- predict_nobias_km(object=model,newdata=mat,type="UK",se.compute=TRUE,cov.compute=TRUE)
      
      mn1 <- c(mn1,obj$mean[1:size])
      sn1 <- c(sn1,obj$sd[1:size])
      
      mn2 <- c(mn2,obj$mean[(size+1):(2*size)])
      sn2 <- c(sn2,obj$sd[(size+1):(2*size)])
      
      for (j in 1:size){#extract kn
        onecov <- obj$cov[j,size+j]
        kn <- c(kn,onecov)
      }
      first <- last+1
      last <- min(n,last+batch)
      size <- 1 + last - first
      indices <- c(first:last)
    }
    #we have mn1,mn2, sn1,sn2 and kn
    pn1 <- pnorm((mn1 - T)/sn1);pn2 <- pnorm((mn2 - T)/sn2)
    a <- (T - mn1)/sn1; b <- (T - mn2)/sn2; 
	
	a[a==Inf]<- 1000 ;a[a== -Inf] <- 1000;b[b==Inf]<-0;b[b== -Inf]<-0
	a[is.nan(a)] <- 1000; b[is.nan(b)]<-0
	
    rho <- pmin(1,kn /(sn1*sn2));rho <- pmax(rho,-1)
	
  	pn.12 <- pbivnorm(a,b,rho)
    
    constant.tab <- pn1 + pn2 + pn.12 - 1
    res <- sum(constant.tab*integration.weights)
    return(res)
    
  }
}

