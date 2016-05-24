#wrapper function
fit.CLMM.2 <- function(data.y1, data.x1, data.z1, data.y2, data.x2, data.z2, n.clst, n.run=1){
	na.flag1<-data.y1;
	na.flag1[!is.na(na.flag1)]=FALSE;
	na.flag1[is.na(na.flag1)]=TRUE;
	na.flag2<-data.y2;
	na.flag2[!is.na(na.flag2)]=FALSE;
	na.flag2[is.na(na.flag2)]=TRUE;
	if((TRUE%in%na.flag1)|(TRUE%in%na.flag2)){
		time.NA1=apply(data.y1,c(1,2),function(x){if(FALSE%in%is.na(x)){y=0}else{y=1};return(y)});
		sample.NA1=apply(data.y1,c(1,3),function(x){if(FALSE%in%is.na(x)){y=0}else{y=1};return(y)});
		gene.NA1=apply(data.y1,c(2,3),function(x){if(FALSE%in%is.na(x)){y=0}else{y=1};return(y)});
		time.NA2=apply(data.y1,c(1,2),function(x){if(FALSE%in%is.na(x)){y=0}else{y=1};return(y)});
		sample.NA2=apply(data.y1,c(1,3),function(x){if(FALSE%in%is.na(x)){y=0}else{y=1};return(y)});
		gene.NA2=apply(data.y1,c(2,3),function(x){if(FALSE%in%is.na(x)){y=0}else{y=1};return(y)});
		if((1%in%time.NA1)|(1%in%sample.NA1)|(1%in%gene.NA1)|(1%in%time.NA2)|(1%in%sample.NA2)|(1%in%gene.NA2)){
			stop("All NAs are found in a pair of observations. Please recheck your input.");
		}else{
			if(TRUE%in%na.flag1){
				index.NA1=which(is.na(data.y1),arr.ind=TRUE);
				for(i in 1:nrow(index.NA1)){
					temp=data.y1[index.NA1[i,1],,index.NA1[i,3]];
					data.y1[index.NA1[i,1],index.NA1[i,2],index.NA1[i,3]]=mean(temp[!is.na(temp)]);
				};
			};
			if(TRUE%in%na.flag2){
				index.NA2=which(is.na(data.y2),arr.ind=TRUE);
				for(i in 1:nrow(index.NA2)){
					temp=data.y2[index.NA2[i,1],,index.NA2[i,3]];
					data.y2[index.NA2[i,1],index.NA2[i,2],index.NA2[i,3]]=mean(temp[!is.na(temp)]);
				};
			};
			fit.CLMM.simple.2data.NA(data.y1=data.y1,data.x1=data.x1,data.z1=data.x1,
                                     data.y2=data.y2,data.x2=data.x2,data.z2=data.x2,
					                 na.flag1=na.flag1,na.flag2=na.flag2,n.clst=n.clst,n.run=n.run);
		}
	}else{
		fit.CLMM.simple.2data(data.y1=data.y1,data.x1=data.x1,data.z1=data.x1,
                              data.y2=data.y2,data.x2=data.x2,data.z2=data.x2,
					          n.clst=n.clst,n.run=n.run);
	}
}
#############################################################################
#####							      				#####
###									  			  ###
#	      fit a CLMM model with ONE common clustering	  	    	    #
#	       	  and SAMPLE-specific covariates x_i	  	  	    #
#								  	    			    #
#								  	    			    #
#	    UPDATE: 1. compute "u.hat" slightly differently	  	          #
#		    2. D can ONLY be "Diagonal"	  	    		    	    #
#	    	    3. ALLOW FOR CLUSTER-SPECIFIC VARIANCE		    	    #
###									  			  ###
#####							      	        		#####
#############################################################################

##### INPUT
###
##	data.y - matrix of observations,
##		  data.y[j, i, l] for sample i, gene j, and measurement l;
###
## 	data.x - design matrix for fixed effects (common for all genes),
##		  data.x[i, l, p] for sample i, measurement l, and covariate p,
###
##	data.z - design matrix for random effects (common for all genes),	
##		  data.z[i, l, q] for sample i, measurement l, and covariate p,
###
##	n.clst - number of clusters for the beta associated with data.x;
###
##	var.type - type of the variance matrix for random effects; 
##		    either "D" for "Diagonal" or "G" for "General"
###
##	n.run - fit the CLMM "n.run" times, each with a different staring value
###
##### OUTPUT (a list of)
###
##	theta.hat - regression parameters estimated via EM algorithm;
###
##	data.u - clustering associated with data.x
##
###

fit.CLMM.simple.2data 	<- function(data.y1, data.x1, data.z1, data.y2, data.x2, data.z2, n.clst, n.run=1){
   ### this will be repeatedly used by M-step
   data.x.x.sum 	<- compute.x.z.sum.CLMM.simple.2(data.x1, data.x1) + compute.x.z.sum.CLMM.simple.2(data.x2, data.x2)	

   ### try different starting values
   llh 			<- -9999999999
   for(s in 1:n.run){
	### get "start values"
	theta.hat 		<-  fit.CLMM.simple.start.2(data.x1, data.y1, data.z1, data.x2, data.y2, data.z2, n.clst, start=s)

	### iterate btw E- and M- steps
	est.hat.new 	<- fit.CLMM.simple.EM.2(data.x1, data.y1, data.z1,data.x2, data.y2, data.z2, data.x.x.sum, theta.hat)

	if(est.hat.new$theta.hat$llh > llh){
	   est.hat 		<- est.hat.new
	   llh 		<- est.hat.new$theta.hat$llh
	}
   }
   return(est.hat)
}


### find a starting value for "zeta" in CLMM
library(cluster)
fit.CLMM.simple.start.2 <- function(data.x1, data.y1, data.z1, data.x2, data.y2, data.z2, n.clst, start=1){
   # number of genes and covariates
   J 	<- dim(data.y1)[1]
   m1 <- dim(data.x1)[1]
   m2 <- dim(data.x2)[1]
   P 	<- dim(data.x1)[3]
   Q 	<- dim(data.z1)[3]

   # get beta.hat for each gene-specific model
   beta.hat 		<- matrix(0, nrow=J, ncol=P)
   temp.x 			<- data.x1[1,,]
   if(m1>1){for(i in 2:m1){
	temp.x 		<- rbind(temp.x, data.x1[i,,])
   }}
   if(m2>0){for(i in 1:m2){
	temp.x 		<- rbind(temp.x, data.x2[i,,])
   }}
   temp <- solve(t(temp.x) %*% temp.x) %*% t(temp.x)
   for(j in 1:J){
	beta.hat[j,] 	<- temp %*% c(as.vector(t(data.y1[j,,])),as.vector(t(data.y2[j,,])))
   }

   # pick group centers randomly if "start>1"
   if(start>1) {
	temp 			<- sample(J, n.clst)
	zeta.hat 		<- beta.hat[temp,]
   }

   # group gene-specific beta.hat's by PAM if "start<=1"
   if(start<=1) {
	temp 			<- pam(beta.hat, n.clst)
	# assign group centers to zeta.hat (a KxP matrix)
	zeta.hat 		<- temp$medoids
   }

   # covariance matrix for random effects
   D.hat 			<- rep(1, n.clst) 

   # measurement error
   sigma2.hat 		<- rep(1, n.clst)

   # frequency of each cluster
   pi.hat 			<- rep(1/n.clst, n.clst)

   return(list(zeta.hat=zeta.hat, D.hat=D.hat, sigma2.hat=sigma2.hat, pi.hat=pi.hat))
}


### EM algorithm to fit the CLMM
###
## data.x[i,l,p]
## data.y[j,i,l]
## data.z[i,l,q]
###
## zeta.hat, D.hat, sigma2.hat, and pi.hat are the starting values for the parameters
###

fit.CLMM.simple.EM.2 	<- function(data.x1, data.y1, data.z1, data.x2, data.y2, data.z2, data.x.x.sum, theta.hat){
   # number of genes, samples, covariates, and clusters
   J 	<- dim(data.y1)[1]
   m1 <- dim(data.y1)[2]
   L1 <- dim(data.y1)[3]
   m2 <- dim(data.y2)[2]
   L2 <- dim(data.y2)[3]
   P 	<- dim(data.x1)[3]
   Q 	<- dim(data.z1)[3]
   K 	<- length(theta.hat$pi.hat)

   # "log likelihood"
   llh.old 			<- -9999999999
   llh 			<- -9999999990

   ### iterate btw E- and M- steps
   while(llh-llh.old>0.01){
	# E-step
	hats 			<- compute.Ehats.CLMM.simple.2(data.x1, data.y1, data.z1, data.x2, data.y2, data.z2, 
						theta.hat, J, m1, L1, m2, L2, K, Q)

	# M-step
	temp 			<- compute.theta.hat.CLMM.simple.2(data.x1, data.y1, data.z1, data.x2, data.y2, data.z2,
						data.x.x.sum, hats, J, m1, L1, m2, L2, K, P, Q)

	# update only when llh increases; when some cluster disappears, llh decreases
	llh.old 		<- llh
	if(!is.na(temp$llh) & temp$llh > llh){
	   theta.hat 	<- temp
	   llh 		<- temp$llh
	   print(llh)
	}
   } 
   return(list(u.hat=hats$u.hat, b.hat.1=hats$b.hat.1, b.hat.2=hats$b.hat.2, theta.hat=theta.hat))
}


### compute "u.hat" - the expected clustering indicator
compute.Ehats.CLMM.simple.2 <- function(data.x1, data.y1, data.z1, data.x2, data.y2, data.z2, 
					theta.hat, J, m1, L1, m2, L2, K, Q){
   # delist "theta.hat"
   zeta.hat 		<- theta.hat$zeta.hat
   D.hat 			<- theta.hat$D.hat
   sigma2.hat 		<- theta.hat$sigma2.hat
   pi.hat 			<- theta.hat$pi.hat

   # compute the inverse of V
   V.inv.1 			<- compute.V.inv.simple.2(data.z1, D.hat, sigma2.hat, J, m1, L1)
   V.inv.2 			<- compute.V.inv.simple.2(data.z2, D.hat, sigma2.hat, J, m2, L2)

   # compute the partial residuals "gene by gene" - pResid[j, i, k, l]
   pResid.1 		<- compute.pResid.CLMM.simple.2(data.x1, data.y1, zeta.hat, J, m1, L1, K)
   pResid.2 		<- compute.pResid.CLMM.simple.2(data.x2, data.y2, zeta.hat, J, m2, L2, K)

   # compute "u.hat"
   u.hat 			<- compute.u.hat.CLMM.simple.2(pResid.1, V.inv.1, pResid.2, V.inv.2, pi.hat, J, m1, m2, K)

   # compute "b.hat", bhat[j, i, k, q]
   b.hat.1 			<- compute.b.hat.CLMM.simple.2(data.z1, pResid.1, D.hat, V.inv.1, J, m1, Q, K)
   b.hat.2 			<- compute.b.hat.CLMM.simple.2(data.z2, pResid.2, D.hat, V.inv.2, J, m2, Q, K)

   # compute "b2.hat", b2.hat[j, k, q, q]
   b2.hat 			<- compute.b2.hat.CLMM.simple.2(data.z1, data.z2, u.hat, b.hat.1, b.hat.2, D.hat, 
						V.inv.1, V.inv.2, J, m1, m2, Q, K)

   # compute "e.hat", e.hat[j, i, k, l]
   e.hat.1 			<- compute.e.hat.CLMM.simple.2(data.z1, b.hat.1, pResid.1, J, m1, L1, K)
   e.hat.2 			<- compute.e.hat.CLMM.simple.2(data.z2, b.hat.2, pResid.2, J, m2, L2, K)

   # compute "b2.hat", b2.hat[j, k]
   e2.hat 			<- compute.e2.hat.CLMM.simple.2(data.z1, e.hat.1, data.z2, e.hat.2, 
						D.hat, V.inv.1, V.inv.2, sigma2.hat, J, m1, L1, m2, L2, K)

   return(list(u.hat=u.hat, b.hat.1=b.hat.1, b.hat.2=b.hat.2, b2.hat=b2.hat, e2.hat=e2.hat))
}


### M-step for fitting CLMM
###
### INPUT
##
#   data.x.x[j,,] = t(data.x[j,,])%*%(data.x[j,,])
#   data.x.y[j,] = t(data.x[j,,])%*%(data.y[j,])
##
#   u.hat is a J*K matrix of cluster membership probabilities
##
###
compute.theta.hat.CLMM.simple.2 <- function(data.x1, data.y1, data.z1, data.x2, data.y2, data.z2,
						data.x.x.sum, hats, J, m1, L1, m2, L2, K, P, Q){
   # delist "hats"
   u.hat 		<- hats$u.hat
   b.hat.1		<- hats$b.hat.1
   b.hat.2		<- hats$b.hat.2
   b2.hat 		<- hats$b2.hat
   e2.hat 		<- hats$e2.hat

   # frequency of each cluster
   pi.hat 		<- apply(u.hat, FUN=sum, MARGIN=2)/J

   # cluster-specific regression parameters
   zeta.hat 		<- matrix(0, nrow=K, ncol=P)
   D.hat		<- rep(0, K)
   sigma2.hat	<- rep(0, K)
   for(k in 1:K){
	zeta.num 	<- 0
	for(j in 1:J){
	   for(i in 1:m1){
		zeta.num<- zeta.num + u.hat[j,k]*t(data.x1[i,,])%*%(data.y1[j,i,]-as.matrix(data.z1[i,,])%*%b.hat.1[j,i,,k])
	   }  
	   for(i in 1:m2){
		zeta.num<- zeta.num + u.hat[j,k]*t(data.x2[i,,])%*%(data.y2[j,i,]-as.matrix(data.z2[i,,])%*%b.hat.2[j,i,,k])
	   }  
	}
	zeta.den 		<- sum(u.hat[,k])*data.x.x.sum
	zeta.hat[k,] 	<- solve(zeta.den) %*% zeta.num

	D.hat[k]		<- sum(u.hat[,k]*b2.hat[,k])/((m1+m2)*Q*sum(u.hat[,k]))

	sigma2.hat[k] 	<- sum(u.hat[,k]*e2.hat[,k])/((m1*L1+m2*L2)*sum(u.hat[,k]))
   }

   # compute the inverse of V
   V.inv.1 			<- compute.V.inv.simple.2(data.z1, D.hat, sigma2.hat, J, m1, L1)
   V.inv.2 			<- compute.V.inv.simple.2(data.z2, D.hat, sigma2.hat, J, m2, L2)

   # compute the "log likelihood" given this MLE
   pResid.1 		<- compute.pResid.CLMM.simple.2(data.x1, data.y1, zeta.hat, J, m1, L1, K)
   pResid.2 		<- compute.pResid.CLMM.simple.2(data.x2, data.y2, zeta.hat, J, m2, L2, K)
   llh 			<- compute.llh.CLMM.2(pResid.1, V.inv.1, pResid.2, V.inv.2, pi.hat, J, m1, m2, K)

   return(list(zeta.hat=zeta.hat, pi.hat=pi.hat, D.hat=D.hat, sigma2.hat=sigma2.hat, llh=llh))
}

### data.x[i,l,p], data.z[i,l,q]
compute.x.z.sum.CLMM.simple.2 <- function(data.x, data.z){
   m 	<- dim(data.x)[1]
   P 	<- dim(data.x)[3]
   Q 	<- dim(data.z)[3]

   data.x.z.sum 		<- matrix(0, nrow=P, ncol=Q)
   for(i in 1:m){
	data.x.z.sum 	<- data.x.z.sum + t(data.x[i,,]) %*% data.z[i,,]
   }
   return(data.x.z.sum)
}

### compute the inverse of V.inv[j, i, l, l, k]
compute.V.inv.simple.2 	<- function(data.z, D.hat, sigma2.hat, J, m, L){
   K				<- length(sigma2.hat)
   V.inv 			<- array(0, dim=c(J, m, L, L, K))
   for(k in 1:K){
	temp1			<- diag(sigma2.hat[k], L)
	for(i in 1:m){
	   temp2 		<- data.z[i,,]
	   for(j in 1:J){
	   	temp 		<- temp1 + D.hat[k]*(temp2%*%t(temp2))
	   	V.inv[j,i,,,k] 	<- solve(temp, tol=1e-50)
	   }
	}
   }
   return(V.inv)
}

### compute the partial residuals "gene by gene" - pResid[j,i,l,k]
compute.pResid.CLMM.simple.2 	<- function(data.x, data.y, zeta.hat, J, m, L, K){
   pResid 				<- array(0, dim=c(J, m, L, K))
   for(k in 1:K){
	for(i in 1:m){
	   y.hat 			<- data.x[i,,] %*% zeta.hat[k,] 	# L x 1
	   pResid[,i,,k]		<- t(t(data.y[,i,]) - as.vector(y.hat))	# J x L
	}
   }
   return(pResid)
}

### compute "b.hat", bhat[j, i, q, k]
compute.b.hat.CLMM.simple.2 	<- function(data.z, pResid, D.hat, V.inv, J, m, Q, K){
   b.hat 				<- array(0, dim=c(J, m, Q, K))
   for(j in 1:J)
	for(i in 1:m)
	   for(k in 1:K)
	   	b.hat[j,i,,k] 	<- D.hat[k]*t(data.z[i,,])%*%V.inv[j,i,,,k]%*%pResid[j,i,,k]

   return(b.hat)
}

### compute "b2.hat", b2.hat[j, k]
### if the variance matrix is a Diagonal matrix
compute.b2.hat.CLMM.simple.2	<- function(data.z1, data.z2, u.hat, b.hat.1, b.hat.2, D.hat, V.inv.1, V.inv.2, 
						J, m1, m2, Q, K){
	b2.hat 			<- array(0, dim=c(J, K))
	for(j in 1:J){
	   for(k in 1:K){
		temp1 		<- 0
	   	temp2 		<- 0
		tau2		<- D.hat[k]
	   	for(i in 1:m1){
		   temp1 		<- temp1 + sum(b.hat.1[j,i,,k]^2)
		   temp2 		<- temp2 + sum(diag(t(data.z1[i,,])%*%V.inv.1[j,i,,,k]%*%data.z1[i,,]))
	   	}
	   	for(i in 1:m2){
		   temp1 		<- temp1 + sum(b.hat.2[j,i,,k]^2)
		   temp2 		<- temp2 + sum(diag(t(data.z2[i,,])%*%V.inv.2[j,i,,,k]%*%data.z2[i,,]))
	   	}
	   	b2.hat[j,k]		<- temp1 + tau2*Q*(m1+m2) - tau2^2*temp2
	   }
	}
   return(b2.hat)
}

### compute "e.hat", e.hat[j, i, l, k]
compute.e.hat.CLMM.simple.2 	<- function(data.z, b.hat, pResid, J, m, L, K){
   e.hat 				<- array(0, dim=c(J,m,L,K))
   for(i in 1:m)
	for(j in 1:J){
	   e.hat[j,i,,] 		<- pResid[j,i,,] - as.matrix(data.z[i,,])%*%b.hat[j,i,,]
	}
   return(e.hat)
}

### compute "e2.hat", e2.hat[j, k]
compute.e2.hat.CLMM.simple.2 	<- function(data.z1, e.hat.1, data.z2, e.hat.2, 
					D.hat, V.inv.1, V.inv.2, sigma2.hat, J, m1, L1, m2, L2, K){
	e2.hat 			<- array(0, dim=c(J, K))
	for(j in 1:J){
	   for(k in 1:K){
		temp1 		<- 0
	   	temp2 		<- 0
	   	for(i in 1:m1){
		   temp1 		<- temp1 + sum(e.hat.1[j,i,,k]^2)
		   temp2 		<- temp2 + sum(diag(V.inv.1[j,i,,,k]))
	   	}
	   	for(i in 1:m2){
		   temp1 		<- temp1 + sum(e.hat.2[j,i,,k]^2)
		   temp2 		<- temp2 + sum(diag(V.inv.2[j,i,,,k]))
	   	}
	   	e2.hat[j,k]		<- temp1 + sigma2.hat[k]*(L1*m1+L2*m2) - (sigma2.hat[k])^2*temp2
	   }
	}
   return(e2.hat)
}

### compute the log-transformed numerator for "u.hat"
compute.u.hat.CLMM.simple.2	<- function(pResid.1, V.inv.1, pResid.2, V.inv.2, pi.hat, J, m1, m2, K){
   # compute log(u.hat.numerator)
   log.u.hat.num 			<- matrix(0, nrow=J, ncol=K)
   for(k in 1:K){
	for(j in 1:J){
	   temp 			<- 0
	   for(i in 1:m1){
		temp 			<- temp + mvn.dnorm.log.2(pResid.1[j,i,,k], V.inv.1[j,i,,,k])
	   }
	   for(i in 1:m2){
		temp 			<- temp + mvn.dnorm.log.2(pResid.2[j,i,,k], V.inv.2[j,i,,,k])
	   }
	   log.u.hat.num[j,k] 	<- log(pi.hat[k]) + temp
	}
   }
   u.hat.num 			<- exp(t(apply(log.u.hat.num, MARGIN=1, FUN=all.ceiling.2)))
   u.hat.den 			<- apply(u.hat.num, FUN=sum, MARGIN=1)
   u.hat 				<- u.hat.num/u.hat.den

   return(u.hat)
}

### substract a constant from a vector to make its max = cutoff
all.ceiling.2 	<- function(aVector, cutoff=600){
    xx 		<- max(aVector)
    aVector 	<- aVector - xx + cutoff
    return(aVector)
}

### compute the density of a vector according to a MVN distribution
mvn.dnorm.log.2 	<- function(aVector, var.inv){
   L 			<- length(aVector)
   y 			<- matrix(aVector, nrow=L)
   # dens 		<- (abs(det(var.inv))^(1/2)) * ((2*pi)^(-L/2)) * exp(-(1/2) * t(y) %*% var.inv %*% y)

   log.dens 	<- log(abs(det(var.inv)))/2 + (log(2*pi)*(-L/2)) + ((-1/2) * t(y) %*% var.inv %*% y)

   return(log.dens)
}


### compute the "log likelihood" given this MLE
compute.llh.CLMM.2 		<- function(pResid.1, V.inv.1, pResid.2, V.inv.2, pi.hat, J, m1, m2, K){
   llh 			<- 0
   for(j in 1:J){
	temp.log 		<- rep(0, K)
	for(k in 1:K){
	   temp.ind 	<- 0
	   for(i in 1:m1){
	       temp.ind 	<- temp.ind + mvn.dnorm.log.2(pResid.1[j,i,,k], V.inv.1[j,i,,,k])
	   }
	   for(i in 1:m2){
	       temp.ind 	<- temp.ind + mvn.dnorm.log.2(pResid.2[j,i,,k], V.inv.2[j,i,,,k])
	   }
	   temp.log[k] 	<- temp.ind  # exp(temp.ind) = "inf" if temp.ind < -700
	}
	temp.log.max 	<- max(temp.log)
	temp 			<-  exp(temp.log - temp.log.max) %*% pi.hat
	llh 			<- llh + log(temp) + temp.log.max
   }
   return(llh)
}

############			the end 			##############

#############################################################################
#####							      				#####
###									  			  ###
#	      fit a CLMM model with ONE common clustering	  	    	    #
#	       	  and SAMPLE-specific covariates x_i	  	  	    #
#								  	    			    #
#								  	    			    #
#	    UPDATE: 1. compute "u.hat" slightly differently	  	    	    #
#		    2. D can ONLY be "Diagonal"	  	    		    	    #
#	    	    3. ALLOW FOR CLUSTER-SPECIFIC VARIANCE		    	    #
###									  			  ###
#####							      	        		#####
#############################################################################

##### INPUT
###
##	data.y - matrix of observations,
##		  data.y[j, i, l] for sample i, gene j, and measurement l;
###
## 	data.x - design matrix for fixed effects (common for all genes),
##		  data.x[i, l, p] for sample i, measurement l, and covariate p,
###
##	data.z - design matrix for random effects (common for all genes),	
##		  data.z[i, l, q] for sample i, measurement l, and covariate p,
###
##	n.clst - number of clusters for the beta associated with data.x;
###
##	var.type - type of the variance matrix for random effects; 
##		    either "D" for "Diagonal" or "G" for "General"
###
##	n.run - fit the CLMM "n.run" times, each with a different staring value
###
##### OUTPUT (a list of)
###
##	theta.hat - regression parameters estimated via EM algorithm;
###
##	data.u - clustering associated with data.x
##
###

fit.CLMM.simple.2data.NA <- function(data.y1,data.x1,data.z1,data.y2,data.x2,data.z2,na.flag1,na.flag2,n.clst,n.run=1){
   ### this will be repeatedly used by M-step
   data.x.x.sum 		<- compute.x.z.sum.CLMM.simple.NA.2(data.x=data.x1, data.z=data.x1, na.flag=na.flag1)+compute.x.z.sum.CLMM.simple.NA.2(data.x=data.x2,data.z=data.x2, na.flag=na.flag2)	

   ### try different starting values
   llh 			<- -9999999999
   for(s in 1:n.run){
	### get "start values"
	theta.hat 		<-  fit.CLMM.simple.start.NA.2(data.x1, data.y1, data.z1, data.x2, data.y2, data.z2, n.clst, start=s)

	### iterate btw E- and M- steps
	est.hat.new 	<- fit.CLMM.simple.EM.NA.2(data.x1, data.y1, data.z1,data.x2, data.y2, data.z2, na.flag1, na.flag2, data.x.x.sum, theta.hat)

	if(est.hat.new$theta.hat$llh > llh){
	   est.hat 		<- est.hat.new
	   llh 		<- est.hat.new$theta.hat$llh
	}
   }
   return(est.hat)
}


### find a starting value for "zeta" in CLMM
library(cluster)
fit.CLMM.simple.start.NA.2 	<- function(data.x1, data.y1, data.z1, data.x2, data.y2, data.z2, n.clst, start=1){
   # number of genes and covariates
   J 	<- dim(data.y1)[1]
   m1 <- dim(data.x1)[1]
   m2 <- dim(data.x2)[1]
   P 	<- dim(data.x1)[3]
   Q 	<- dim(data.z1)[3]

   # get beta.hat for each gene-specific model
   beta.hat 		<- matrix(0, nrow=J, ncol=P)
   temp.x 			<- data.x1[1,,]
   if(m1>1){for(i in 2:m1){
	temp.x 		<- rbind(temp.x, data.x1[i,,])
   }}
   if(m2>0){for(i in 1:m2){
	temp.x 		<- rbind(temp.x, data.x2[i,,])
   }}
   temp <- solve(t(temp.x) %*% temp.x) %*% t(temp.x)
   for(j in 1:J){
	beta.hat[j,] 	<- temp %*% c(as.vector(t(data.y1[j,,])),as.vector(t(data.y2[j,,])))
   }

   # pick group centers randomly if "start>1"
   if(start>1) {
	temp 			<- sample(J, n.clst)
	zeta.hat 		<- beta.hat[temp,]
   }

   # group gene-specific beta.hat's by PAM if "start<=1"
   if(start<=1) {
	temp 			<- pam(beta.hat, n.clst)
	# assign group centers to zeta.hat (a KxP matrix)
	zeta.hat 		<- temp$medoids
   }

   # covariance matrix for random effects
   D.hat 		<- rep(1, n.clst) 

   # measurement error
   sigma2.hat 		<- rep(1, n.clst)

   # frequency of each cluster
   pi.hat 			<- rep(1/n.clst, n.clst)

   return(list(zeta.hat=zeta.hat, D.hat=D.hat, sigma2.hat=sigma2.hat, pi.hat=pi.hat))
}


### EM algorithm to fit the CLMM
###
## data.x[i,l,p]
## data.y[j,i,l]
## data.z[i,l,q]
###
## zeta.hat, D.hat, sigma2.hat, and pi.hat are the starting values for the parameters
###

fit.CLMM.simple.EM.NA.2 	<- function(data.x1, data.y1, data.z1, data.x2, data.y2, data.z2, 
						na.flag1, na.flag2, data.x.x.sum, theta.hat){
   # number of genes, samples, covariates, and clusters
   J 	<- dim(data.y1)[1]
   m1 <- dim(data.y1)[2]
   L1 <- dim(data.y1)[3]
   m2 <- dim(data.y2)[2]
   L2 <- dim(data.y2)[3]
   P 	<- dim(data.x1)[3]
   Q 	<- dim(data.z1)[3]
   K 	<- length(theta.hat$pi.hat)

   # "log likelihood"
   llh.old 			<- -9999999999
   llh 			<- -9999999990

   ### iterate btw E- and M- steps
   while(llh-llh.old>0.1){
	# E-step
	hats <- compute.Ehats.CLMM.simple.NA.2(data.x1, data.y1, data.z1, data.x2, data.y2, data.z2, na.flag1, na.flag2, 
						theta.hat, J, m1, L1, m2, L2, K, Q)

	# M-step
	temp <- compute.theta.hat.CLMM.simple.NA.2(data.x1, data.y1, data.z1, data.x2, data.y2, data.z2, na.flag1, na.flag2,
						data.x.x.sum, hats, J, m1, L1, m2, L2, K, P, Q)

	# update only when llh increases; when some cluster disappears, llh decreases
	llh.old 		<- llh
	if(!is.na(temp$llh) & temp$llh > llh){
	   theta.hat 	<- temp
	   llh 		<- temp$llh
	   print(llh)
	}
   } 
   return(list(u.hat=hats$u.hat, b.hat.1=hats$b.hat.1, b.hat.2=hats$b.hat.2, theta.hat=theta.hat))
}


### compute "u.hat" - the expected clustering indicator
compute.Ehats.CLMM.simple.NA.2 <- function(data.x1, data.y1, data.z1, data.x2, data.y2, data.z2, na.flag1, na.flag2, 
					theta.hat, J, m1, L1, m2, L2, K, Q){
   # delist "theta.hat"
   zeta.hat 	<- theta.hat$zeta.hat
   D.hat 		<- theta.hat$D.hat
   sigma2.hat 	<- theta.hat$sigma2.hat
   pi.hat 		<- theta.hat$pi.hat

   # compute the inverse of V
   V.1 		<- compute.V.simple.NA.2(data.z1, D.hat, sigma2.hat, J, m1, L1)
   V.2 		<- compute.V.simple.NA.2(data.z2, D.hat, sigma2.hat, J, m2, L2)

   # compute the partial residuals "gene by gene" - pResid[j, i, k, l]
   pResid.1 	<- compute.pResid.CLMM.simple.NA.2(data.x1, data.y1, zeta.hat, J, m1, L1, K)
   pResid.2 	<- compute.pResid.CLMM.simple.NA.2(data.x2, data.y2, zeta.hat, J, m2, L2, K)

   # compute "u.hat"
   u.hat 		<- compute.u.hat.CLMM.simple.NA.2(pResid.1, V.1, pResid.2, V.2, pi.hat, J, m1, m2, K, na.flag1, na.flag2)

   # compute "b.hat", bhat[j, i, k, q]
   b.hat.1 		<- compute.b.hat.CLMM.simple.NA.2(data.z1, pResid.1, D.hat, V.1, J, m1, Q, K, na.flag1)
   b.hat.2 		<- compute.b.hat.CLMM.simple.NA.2(data.z2, pResid.2, D.hat, V.2, J, m2, Q, K, na.flag2)

   # compute "b2.hat", b2.hat[j, k, q, q]
   b2.hat 		<- compute.b2.hat.CLMM.simple.NA.2(data.z1, data.z2, u.hat, b.hat.1, b.hat.2, D.hat, V.1, V.2, 
					J, m1, m2, Q, K, na.flag1, na.flag2)

   # compute "e.hat", e.hat[j, i, k, l]
   e.hat.1 		<- compute.e.hat.CLMM.simple.NA.2(data.z1, b.hat.1, pResid.1, J, m1, L1, K)
   e.hat.2 		<- compute.e.hat.CLMM.simple.NA.2(data.z2, b.hat.2, pResid.2, J, m2, L2, K)

   # compute "b2.hat", b2.hat[j, k]
   e2.hat 		<- compute.e2.hat.CLMM.simple.NA.2(data.z1, e.hat.1, data.z2, e.hat.2, 
					D.hat, V.1, V.2, sigma2.hat, J, m1, L1, m2, L2, K, na.flag1, na.flag2)

   return(list(u.hat=u.hat, b.hat.1=b.hat.1, b.hat.2=b.hat.2, b2.hat=b2.hat, e2.hat=e2.hat))
}


### M-step for fitting CLMM
###
### INPUT
##
#   data.x.x[j,,] = t(data.x[j,,])%*%(data.x[j,,])
#   data.x.y[j,] = t(data.x[j,,])%*%(data.y[j,])
##
#   u.hat is a J*K matrix of cluster membership probabilities
##
###
compute.theta.hat.CLMM.simple.NA.2 <- function(data.x1, data.y1, data.z1, data.x2, data.y2, data.z2, na.flag1, na.flag2,
						data.x.x.sum, hats, J, m1, L1, m2, L2, K, P, Q){
   # delist "hats"
   u.hat 			<- hats$u.hat
   b.hat.1			<- hats$b.hat.1
   b.hat.2			<- hats$b.hat.2
   b2.hat 			<- hats$b2.hat
   e2.hat 			<- hats$e2.hat

   # frequency of each cluster
   pi.hat 			<- apply(u.hat, FUN=sum, MARGIN=2)/J

   # cluster-specific regression parameters
   zeta.hat 		<- matrix(0, nrow=K, ncol=P)
   D.hat			<- rep(0, K)
   sigma2.hat		<- rep(0, K)
   for(k in 1:K){
	zeta.num 		<- 0
	zeta.den		<- 0
	for(j in 1:J){
	   for(i in 1:m1){
		temp1		<- !(na.flag1[j,i,])
		temp2		<- t(data.x1[i,temp1,])%*%(data.y1[j,i,temp1]-as.matrix(data.z1[i,temp1,])%*%b.hat.1[j,i,,k])
		zeta.num	<- zeta.num + u.hat[j,k]*temp2
	   }  
	   for(i in 1:m2){
		temp1		<- !(na.flag2[j,i,])
		temp2		<- t(data.x2[i,temp1,])%*%(data.y2[j,i,temp1]-as.matrix(data.z2[i,temp1,])%*%b.hat.2[j,i,,k])
		zeta.num	<- zeta.num + u.hat[j,k]*temp2
	   }  
	   zeta.den 	<- zeta.den + u.hat[j,k]*data.x.x.sum[j,,]
	}

	zeta.hat[k,] 	<- solve(zeta.den, tol=1e-50) %*% zeta.num

	D.hat[k]		<- sum(u.hat[,k]*b2.hat[,k])/((m1+m2)*Q*sum(u.hat[,k]))

	sigma2.hat[k] 	<- sum(u.hat[,k]*e2.hat[,k])/((m1*L1+m2*L2)*sum(u.hat[,k]))
   }

   # compute the inverse of V
   V.1 			<- compute.V.simple.NA.2(data.z1, D.hat, sigma2.hat, J, m1, L1)
   V.2 			<- compute.V.simple.NA.2(data.z2, D.hat, sigma2.hat, J, m2, L2)

   # compute the "log likelihood" given this MLE
   pResid.1 		<- compute.pResid.CLMM.simple.NA.2(data.x1, data.y1, zeta.hat, J, m1, L1, K)
   pResid.2 		<- compute.pResid.CLMM.simple.NA.2(data.x2, data.y2, zeta.hat, J, m2, L2, K)
   llh 			<- compute.llh.CLMM.NA.2(pResid.1, V.1, pResid.2, V.2, pi.hat, J, m1, m2, K, na.flag1, na.flag2)

   return(list(zeta.hat=zeta.hat, pi.hat=pi.hat, D.hat=D.hat, sigma2.hat=sigma2.hat, llh=llh))
}

### data.x[i,l,p], data.z[i,l,q]
compute.x.z.sum.CLMM.simple.NA.2 	<- function(data.x, data.z, na.flag){
   m <- dim(data.x)[1]
   P <- dim(data.x)[3]
   Q <- dim(data.z)[3]
   J <- dim(na.flag)[1]

   data.x.z.sum 			<- array(0, dim=c(J,P,Q))
   for(j in 1:J){
	for(i in 1:m){
	   temp			<- !(na.flag[j,i,])
	   data.x.z.sum[j,,] 	<- data.x.z.sum[j,,] + t(data.x[i,temp,]) %*% data.z[i,temp,]
	}
   }
   return(data.x.z.sum)
}

### compute V[j, i, l, l, k]
compute.V.simple.NA.2 		<- function(data.z, D.hat, sigma2.hat, J, m, L){
   K				<- length(sigma2.hat)
   V 				<- array(0, dim=c(J, m, L, L, K))
   for(k in 1:K){
	temp1			<- diag(sigma2.hat[k], L)
	for(i in 1:m){
	   temp2 		<- data.z[i,,]
	   for(j in 1:J){
	   	V[j,i,,,k] 	<- temp1 + D.hat[k]*(temp2%*%t(temp2))
	   }
	}
   }
   return(V)
}

### compute the partial residuals "gene by gene" - pResid[j,i,l,k]
compute.pResid.CLMM.simple.NA.2 	<- function(data.x, data.y, zeta.hat, J, m, L, K){
   pResid 				<- array(0, dim=c(J, m, L, K))
   for(k in 1:K){
	for(i in 1:m){
	   y.hat 			<- data.x[i,,] %*% zeta.hat[k,] 	# L x 1
	   pResid[,i,,k]		<- t(t(data.y[,i,]) - as.vector(y.hat))	# J x L
	}
   }
   return(pResid)
}

### compute "b.hat", bhat[j, i, q, k]
compute.b.hat.CLMM.simple.NA.2 	<- function(data.z, pResid, D.hat, V, J, m, Q, K, na.flag){
   b.hat 				<- array(0, dim=c(J, m, Q, K))
   for(j in 1:J){
	for(i in 1:m){
	   temp			<- !(na.flag[j,i,])
	   for(k in 1:K){
	   	b.hat[j,i,,k] 	<- D.hat[k]*t(data.z[i,temp,])%*%solve(V[j,i,temp,temp,k],tol=1e-50)%*%pResid[j,i,temp,k]
	   }
	}
   }
   return(b.hat)
}

### compute "b2.hat", b2.hat[j, k]
### if the variance matrix is a Diagonal matrix
compute.b2.hat.CLMM.simple.NA.2	<- function(data.z1, data.z2, u.hat, b.hat.1, b.hat.2, D.hat, V.1, V.2, 
						J, m1, m2, Q, K, na.flag1, na.flag2){
	b2.hat 			<- array(0, dim=c(J, K))
	for(j in 1:J){
	   for(k in 1:K){
		temp1 		<- 0
	   	temp2 		<- 0
		tau2		<- D.hat[k]
	   	for(i in 1:m1){
		   temp		<- !(na.flag1[j,i,])
		   temp1 		<- temp1 + sum(b.hat.1[j,i,,k]^2)
		   temp2 <- temp2 + sum(diag(t(data.z1[i,temp,])%*%solve(V.1[j,i,temp,temp,k],tol=1e-50)%*%data.z1[i,temp,]))
	   	}
	   	for(i in 1:m2){
		   temp		<- !(na.flag2[j,i,])
		   temp1 		<- temp1 + sum(b.hat.2[j,i,,k]^2)
		   temp2 <- temp2 + sum(diag(t(data.z2[i,temp,])%*%solve(V.2[j,i,temp,temp,k],tol=1e-50)%*%data.z2[i,temp,]))
	   	}
	   	b2.hat[j,k]		<- temp1 + tau2*Q*(m1+m2) - tau2^2*temp2
	   }
	}
   return(b2.hat)
}

### compute "e.hat", e.hat[j, i, l, k]
compute.e.hat.CLMM.simple.NA.2 	<- function(data.z, b.hat, pResid, J, m, L, K){
   e.hat 				<- array(0, dim=c(J,m,L,K))
   for(i in 1:m)
	for(j in 1:J){
	   e.hat[j,i,,] 		<- pResid[j,i,,] - as.matrix(data.z[i,,])%*%b.hat[j,i,,]
	}
   return(e.hat)
}

### compute "e2.hat", e2.hat[j, k]
compute.e2.hat.CLMM.simple.NA.2 	<- function(data.z1, e.hat.1, data.z2, e.hat.2, D.hat, V.1, V.2, sigma2.hat, 
							J, m1, L1, m2, L2, K, na.flag1, na.flag2){
	e2.hat 			<- array(0, dim=c(J, K))
	for(j in 1:J){
	   for(k in 1:K){
		temp1 		<- 0
	   	temp2 		<- 0
	   	for(i in 1:m1){
		   temp		<- !(na.flag1[j,i,])
		   temp1 		<- temp1 + sum(e.hat.1[j,i,temp,k]^2)
		   temp2 		<- temp2 + sum(diag(solve(V.1[j,i,temp,temp,k],tol=1e-50)))
	   	}
	   	for(i in 1:m2){
		   temp		<- !(na.flag2[j,i,])
		   temp1 		<- temp1 + sum(e.hat.2[j,i,temp,k]^2)
		   temp2 		<- temp2 + sum(diag(solve(V.2[j,i,temp,temp,k],tol=1e-50)))
	   	}
	   	e2.hat[j,k]		<- temp1 + sigma2.hat[k]*(L1*m1+L2*m2) - (sigma2.hat[k])^2*temp2
	   }
	}
   return(e2.hat)
}

### compute the log-transformed numerator for "u.hat"
compute.u.hat.CLMM.simple.NA.2	<- function(pResid.1, V.1, pResid.2, V.2, pi.hat, J, m1, m2, K, na.flag1, na.flag2){
   # compute log(u.hat.numerator)
   log.u.hat.num 			<- matrix(0, nrow=J, ncol=K)
   for(k in 1:K){
	for(j in 1:J){
	   temp1 			<- 0
	   for(i in 1:m1){
		temp			<- !(na.flag1[j,i,])
		temp1			<- temp1 + mvn.dnorm.log.NA.2(pResid.1[j,i,temp,k], solve(V.1[j,i,temp,temp,k],tol=1e-50))
	   }
	   for(i in 1:m2){
		temp			<- !(na.flag2[j,i,])
		temp1			<- temp1 + mvn.dnorm.log.NA.2(pResid.2[j,i,temp,k], solve(V.2[j,i,temp,temp,k],tol=1e-50))
	   }
	   log.u.hat.num[j,k] 	<- log(pi.hat[k]) + temp1
	}
   }
   u.hat.num 			<- exp(t(apply(log.u.hat.num, MARGIN=1, FUN=all.ceiling.NA.2)))
   u.hat.den 			<- apply(u.hat.num, FUN=sum, MARGIN=1)
   u.hat 				<- u.hat.num/u.hat.den

   return(u.hat)
}

### substract a constant from a vector to make its max = cutoff
all.ceiling.NA.2 	<- function(aVector, cutoff=600){
    xx 		<- max(aVector)
    aVector 	<- aVector - xx + cutoff
    return(aVector)
}

### compute the density of a vector according to a MVN distribution
mvn.dnorm.log.NA.2 	<- function(aVector, var.inv){
   L 			<- length(aVector)
   y 			<- matrix(aVector, nrow=L)

   # dens 		<- (abs(det(var.inv))^(1/2)) * ((2*pi)^(-L/2)) * exp(-(1/2) * t(y) %*% var.inv %*% y)
   log.dens 	<- log(abs(det(var.inv)))/2 + (log(2*pi)*(-L/2)) + ((-1/2) * t(y) %*% var.inv %*% y)

   return(log.dens)
}


### compute the "log likelihood" given this MLE
compute.llh.CLMM.NA.2 		<- function(pResid.1, V.1, pResid.2, V.2, pi.hat, J, m1, m2, K, na.flag1, na.flag2){
   llh 			<- 0
   for(j in 1:J){
	temp.log 		<- rep(0, K)
	for(k in 1:K){
	   temp.ind 	<- 0
	   for(i in 1:m1){
		temp		<- !(na.flag1[j,i,])
		temp.ind	<- temp.ind + mvn.dnorm.log.NA.2(pResid.1[j,i,temp,k], solve(V.1[j,i,temp,temp,k],tol=1e-50))
	   }
	   for(i in 1:m2){
		temp		<- !(na.flag2[j,i,])
		temp.ind	<- temp.ind + mvn.dnorm.log.NA.2(pResid.2[j,i,temp,k], solve(V.2[j,i,temp,temp,k],tol=1e-50))
	   }
	   temp.log[k] 	<- temp.ind  # exp(temp.ind) = "inf" if temp.ind < -700
	}
	temp.log.max 	<- max(temp.log)
	temp 			<- exp(temp.log - temp.log.max) %*% pi.hat
	llh 			<- llh + log(temp) + temp.log.max
   }
   return(llh)
}

############			the end 			##############
