
moments.differentials<-function(z,W,n.boot=1000,standardized=FALSE){
	
	Pmatrix<-function(z){res<-NULL;if(is.null(dim(z))){
		                   res<-var(z)}else{res<-cov(z)};res;}

	weighted.var<-function(x, w, na.rm = FALSE) {
	    if (na.rm) {
	        w <- w[i <- !is.na(x)]
	        x <- x[i]
	    }
	    sum.w <- sum(w)
	    sum.w2 <- sum(w^2)
	    mean.w <- sum(x * w) / sum(w)
    (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, 
	                                       na.rm = na.rm)
	}

	weighted.covar<-function(x1,x2, w) {
	    w<-w/sum(w)
	    sum(w*(x1-weighted.mean(x1,w))*(x2-weighted.mean(x2,w))) / 
                                                 (1-sum(w^2))
	}

	  traits<-1
	  if(is.null(dim(z)[2])==FALSE) traits<-dim(z)[2]

	  # put it all in a function to ease bootstrapping
	  coef.moments<-function(z,W){	
		z<-as.data.frame(z)
		if(standardized){
		  z<-(z-matrix(sapply(as.data.frame(z),function(x){mean(x)}),
					dim(z)[1],dim(z)[2],byrow=TRUE))/
					matrix(sapply(as.data.frame(z),function(x){
						sd(x)}),dim(z)[1],dim(z)[2],byrow=TRUE)
		}

  	
  		res<-as.data.frame(matrix(NA,2*(2*traits+(traits^2-traits)/2),3))
		names(res)<-c("Coefficient","SE","P-value")
	
		Cmatrix<-matrix(NA,traits,traits)
	
		#differentials
		for(x in 1:traits){
			res[x,1]<-weighted.mean(z[,x],W)-mean(z[,x])
			rownames(res)[x]<-paste("S",x)
			rownames(res)[x+2*traits+(traits^2-traits)/2]<-paste("B",x)
			
			res[traits+x,1]<-weighted.var(z[,x],W)-var(z[,x])+res[x,1]^2
			rownames(res)[traits+x]<-paste("C",x)
			rownames(res)[traits+x+2*traits+
		                       (traits^2-traits)/2]<-paste("G",x)
			Cmatrix[x,x]<-res[traits+x,1]
		}
		if(traits>1){
	 	 index<-2*traits
	 	 for(x in 1:(traits-1)){
	 		for(y in (x+1):traits){
				index<-index+1
				res[index,1]<-weighted.covar(z[,x],z[,y],W)-
			                   cov(z[,x],z[,y])+res[x,1]*res[y,1]
				rownames(res)[index]<-paste("C ",x,",",y,sep="")
				Cmatrix[x,y]<-res[index,1]
				Cmatrix[y,x]<-res[index,1]
				rownames(res)[index+2*traits+(traits^2-traits)/2]<-
			                   paste("G ",x,",",y,sep="")
			}
		  }
		}
	
		#gradients
		res[(2*traits+(traits^2-traits)/2+1):(2*traits+(traits^2-
	                    traits)/2+traits),1]<-
						solve(Pmatrix(z))%*%res[1:traits,1]
		GammaMatrix<-solve(Pmatrix(z))%*%Cmatrix%*%solve(Pmatrix(z))
		if(traits>1)	res[(2*traits+(traits^2-traits)/2+traits+1):
	                    (4*traits+(traits^2-traits)),1]<-
	                    GammaMatrix[lower.tri(GammaMatrix)]
	
	
		res[(2*traits+(traits^2-traits)/2+traits+1):
	       (2*traits+(traits^2-traits)/2+2*traits),1]<-
	       diag(GammaMatrix)
	
		if(traits>1)	res[(2*traits+(traits^2-traits)/2+2*traits+1):
	       (2*(2*traits+(traits^2-traits)/2)),1]<-
	       GammaMatrix[lower.tri(GammaMatrix)]
		res
  	}
  
  	# a function to call coef.moments() in the way required by boot()
  	call.coef.moments<-function(data,original){
  		b<-data[original,]
  		t(coef.moments(z=as.data.frame(b[,1:traits]),
                         W=b[,traits+1])[,1])
  	}

  	# do the bootstrapping  
  	data.frame.for.bootstrapping<-as.data.frame(cbind(z,W))
  	boot.res<-boot(data.frame.for.bootstrapping, call.coef.moments,n.boot)$t
  
  	results<-coef.moments(z,W)
  	results[,2]<-sapply(as.data.frame(boot.res),function(x){sd(x)})
  	results[,3]<-sapply(as.data.frame(boot.res),function(x){
  	                   2*min(sum((x>0)+0),sum((x<0)+0))/n.boot })

  	results[1:(dim(results)[1]/2),]
}
