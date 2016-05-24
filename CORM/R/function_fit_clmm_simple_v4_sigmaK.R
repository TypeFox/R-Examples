#wrapper function
fit.CLMM <- function(data.y, data.x, data.z, n.clst, n.start=1){
	na.flag<-data.y;
	na.flag[!is.na(na.flag)]=FALSE;
	na.flag[is.na(na.flag)]=TRUE;
	if(TRUE%in%na.flag){
		time.NA=apply(data.y,c(1,2),function(x){if(FALSE%in%is.na(x)){y=0}else{y=1};return(y)});
		sample.NA=apply(data.y,c(1,3),function(x){if(FALSE%in%is.na(x)){y=0}else{y=1};return(y)});
		gene.NA=apply(data.y,c(2,3),function(x){if(FALSE%in%is.na(x)){y=0}else{y=1};return(y)});
		if((1%in%time.NA)|(1%in%sample.NA)|(1%in%gene.NA)){
			stop("All NAs are found in a pair of observations. Please recheck your input.");
		}else{
			index.NA=which(is.na(data.y),arr.ind=TRUE);
			for(i in 1:nrow(index.NA)){
				temp=data.y[index.NA[i,1],,index.NA[i,3]];
				data.y[index.NA[i,1],index.NA[i,2],index.NA[i,3]]=mean(temp[!is.na(temp)]);
			};
			fit.CLMM.simple.diag.NA(data.y, data.x, data.z, na.flag, n.clst, n.start);
		}
	}else{
		fit.CLMM.simple.diag(data.y, data.x, data.z, n.clst, n.start);
	}
}
#############################################################################
#####							      				#####
###									  			  ###
#	      fit a CLMM model with ONE common clustering	  	    	    #
#	       	  and SAMPLE-specific covariates x_i	  	  	    #
#								  	    			    #
#								  	    			    #
#	    UPDATE: 1. compute "u.hat" slightly differently	  	    	    #
#		    	2. D can ONLY be "Diagonal"	  	    		    	    #
#	    	    	3. CLUSTER-SPECIFIC MEASUREMENT ERROR	     		    #
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

fit.CLMM.simple.diag 	<- function(data.y, data.x, data.z, n.clst, n.start=1){
   temp			<- dim(data.y)
   if(length(temp)==2) {
	temp1			<- array(0, dim=c(temp[1], 1, temp[2]))
	temp1[,1,]		<- data.y
	data.y		<- temp1
   }
   temp			<- dim(data.x)
   if(length(temp)==2) {
	temp1			<- array(0, dim=c(1, temp))
	temp1[1,,]		<- data.x
	data.x		<- temp1
   }
   temp			<- dim(data.z)
   if(length(temp)==2) {
	temp1			<- array(0, dim=c(1, temp))
	temp1[1,,]		<- data.z
	data.z		<- temp1
   }

   ### this will be repeatedly used by M-step
   data.x.x.sum 		<- compute.x.z.sum.CLMM.simple(data.x, data.x)	# dim is PxP

   ### try different starting values
   llh 			<- -9999999999
   for(s in 1:n.start){
	### get "start values"
	theta.hat 		<-  fit.CLMM.simple.start(data.x, data.y, data.z, n.clst, start=s)$theta.hat

	### iterate btw E- and M- steps
	est.hat.new 	<- fit.CLMM.simple.EM(data.x, data.y, data.z, data.x.x.sum, theta.hat)
	if(est.hat.new$theta.hat$llh > llh){
	   est.hat 		<- est.hat.new
	   llh		<- est.hat.new$theta.hat$llh
	}
   }
   return(est.hat)
}


### find a starting value for "zeta" in CLMM
library(cluster)
fit.CLMM.simple.start 	<- function(data.x, data.y, data.z, n.clst, start=1){
   # number of genes and covariates
   J 				<- dim(data.y)[1]
   m 				<- dim(data.x)[1]
   P 				<- dim(data.x)[3]
   Q 				<- dim(data.z)[3]
   u.hat			<- matrix(0, nrow=J, ncol=n.clst)

   # get beta.hat for each gene-specific model
   beta.hat 		<- matrix(0, nrow=J, ncol=P)
   temp.x 			<- data.x[1,,]
   if(m>1){for(i in 2:m){
	temp.x 		<- rbind(temp.x, data.x[i,,])
   }}
   temp 			<- solve(t(temp.x) %*% temp.x) %*% t(temp.x)
   for(j in 1:J)
	beta.hat[j,] 	<- temp %*% as.vector(t(data.y[j,,]))

   # group gene-specific beta.hat's by PAM if "start==1"
   if(start==1) {
	temp 			<- pam(beta.hat, n.clst)
	# assign group centers to zeta.hat (a KxP matrix)
	zeta.hat 		<- temp$medoids
	temp			<- temp$clustering
   }

   # group gene-specific beta.hat's by PAM if "start==2"
   if(start==2) {
	temp 			<- kmeans(beta.hat, n.clst)
	# assign group centers to zeta.hat (a KxP matrix)
	zeta.hat 		<- temp$centers
 	temp			<- temp$cluster
  }

   # pick group centers randomly if "start>1"
   if(start>2) {
	temp 			<- sample(J, n.clst)
	zeta.hat 		<- beta.hat[temp,]
	temp			<- sample(n.clst, J)
   }

   for(k in 1:n.clst)
      u.hat[temp==k, k]	<- 1

  # covariance matrix for random effects
   D.hat 			<- rep(1, n.clst) 

   # measurement error
   sigma2.hat 		<- rep(1, n.clst)

   # frequency of each cluster
   pi.hat 			<- rep(1/n.clst, n.clst)

   return(list(theta.hat=list(zeta.hat=zeta.hat, D.hat=D.hat, sigma2.hat=sigma2.hat, pi.hat=pi.hat), u.hat=u.hat))
}


### EM algorithm to fit the CLMM
###
## data.x[i,l,p]
## data.y[j,i,l]
## data.z[i,l,q]
###
## zeta.hat, D.hat, sigma2.hat, and pi.hat are the starting values for the parameters
###

fit.CLMM.simple.EM 	<- function(data.x, data.y, data.z, data.x.x.sum, theta.hat){
   # number of genes, samples, covariates, and clusters
   J <- dim(data.y)[1]
   m <- dim(data.y)[2]
   L <- dim(data.y)[3]
   P <- dim(data.x)[3]
   Q <- dim(data.z)[3]
   K <- length(theta.hat$pi.hat)

   # "log likelihood"
   llh.old 			<- -9999999999
   llh 			<- -9999999990

   ### iterate btw E- and M- steps
   while(llh-llh.old>0.1){
	# E-step
	hats 			<- compute.Ehats.CLMM.simple(data.x, data.y, data.z, theta.hat, J, m, L, K, Q)

	# M-step
	temp 			<- compute.theta.hat.CLMM.simple(data.x, data.y, data.z, data.x.x.sum, hats, J, m, L, K, P, Q)

	# update only when llh increases; when some cluster disappears, llh decreases
	llh.old 		<- llh
	if(!is.na(temp$llh) & temp$llh > llh){
	   theta.hat 	<- temp
	   llh 		<- temp$llh
#	   print(llh)
	}
   } 
   return(list(u.hat=hats$u.hat, b.hat=hats$b.hat, theta.hat=theta.hat))
}


### compute "u.hat" - the expected clustering indicator
compute.Ehats.CLMM.simple <- function(data.x, data.y, data.z, theta.hat, J, m, L, K, Q){
   # delist "theta.hat"
   zeta.hat 	<- theta.hat$zeta.hat
   D.hat 		<- theta.hat$D.hat
   sigma2.hat 	<- theta.hat$sigma2.hat
   pi.hat 		<- theta.hat$pi.hat

   # compute the inverse of V
   V.inv 		<- compute.V.inv.simple(data.z, D.hat, sigma2.hat, J, m, L)

   # compute the partial residuals "gene by gene" - pResid[j, i, k, l]
   pResid 		<- compute.pResid.CLMM.simple(data.x, data.y, zeta.hat, J, m, L, K)

   # compute "u.hat"
   u.hat 		<- compute.u.hat.CLMM.simple(pResid, V.inv, pi.hat, J, m, K)

   # compute "b.hat", bhat[j, i, k, q]
   b.hat 		<- compute.b.hat.CLMM.simple(data.z, pResid, D.hat, V.inv, J, m, Q, K)

   # compute "b2.hat", b2.hat[j, k, q, q]
   b2.hat 		<- compute.b2.hat.CLMM.simple(data.z, u.hat, b.hat, D.hat, V.inv, J, m, Q, K)

   # compute "e.hat", e.hat[j, i, k, l]
   e.hat 		<- compute.e.hat.CLMM.simple(data.z, b.hat, pResid, J, m, L, K)

   # compute "b2.hat", b2.hat[j, k]
   e2.hat 		<- compute.e2.hat.CLMM.simple(data.z, u.hat, e.hat, V.inv, sigma2.hat, J, m, L, K)
						
   return(list(u.hat=u.hat, b.hat=b.hat, b2.hat=b2.hat, e2.hat=e2.hat))
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
compute.theta.hat.CLMM.simple <- function(data.x, data.y, data.z, data.x.x.sum, hats, J, m, L, K, P, Q){
   # delist "hats"
   u.hat 			<- hats$u.hat
   b.hat 			<- hats$b.hat
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
	for(j in 1:J){
	   for(i in 1:m){
		zeta.num	<- zeta.num + u.hat[j,k]*t(data.x[i,,])%*%(data.y[j,i,]-as.matrix(data.z[i,,])%*%b.hat[j,i,,k])
	   }  
	}
	zeta.den 		<- sum(u.hat[,k])*data.x.x.sum
	zeta.hat[k,] 	<- solve(zeta.den) %*% zeta.num

	D.hat[k]		<- sum(u.hat[,k]*b2.hat[,k])/(m*Q*sum(u.hat[,k]))

	sigma2.hat[k] 	<- sum(u.hat[,k]*e2.hat[,k])/(m*L*sum(u.hat[,k]))
   }

   # compute the inverse of V
   V.inv 			<- compute.V.inv.simple(data.z, D.hat, sigma2.hat, J, m, L)

   # compute the "log likelihood" given this MLE
   pResid 			<- compute.pResid.CLMM.simple(data.x, data.y, zeta.hat, J, m, L, K)
   llh 			<- compute.llh.CLMM(pResid, V.inv, pi.hat, J, m, K)

   return(list(zeta.hat=zeta.hat, pi.hat=pi.hat, D.hat=D.hat, sigma2.hat=sigma2.hat, llh=llh))
}

### data.x[i,l,p], data.z[i,l,q]
compute.x.z.sum.CLMM.simple 	<- function(data.x, data.z){
   m <- dim(data.x)[1]
   P <- dim(data.x)[3]
   Q <- dim(data.z)[3]

   data.x.z.sum 			<- matrix(0, nrow=P, ncol=Q)
   for(i in 1:m){
	data.x.z.sum 		<- data.x.z.sum + t(data.x[i,,]) %*% data.z[i,,]
   }
   return(data.x.z.sum)
}

### compute the inverse of V.inv[j, i, l, l, k]
compute.V.inv.simple 		<- function(data.z, D.hat, sigma2.hat, J, m, L){
   K					<- length(sigma2.hat)
   V.inv 				<- array(0, dim=c(J, m, L, L, K))
   for(k in 1:K){
	temp1				<- diag(sigma2.hat[k], L)
	for(i in 1:m){
	   temp2 			<- data.z[i,,]
	   for(j in 1:J){
	   	temp 			<- temp1 + D.hat[k]*(temp2%*%t(temp2))
	   	V.inv[j,i,,,k] 	<- solve(temp, tol=1e-50)
	   }
	}
   }
   return(V.inv)
}

### compute the partial residuals "gene by gene" - pResid[j,i,l,k]
compute.pResid.CLMM.simple 	<- function(data.x, data.y, zeta.hat, J, m, L, K){
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
compute.b.hat.CLMM.simple 	<- function(data.z, pResid, D.hat, V.inv, J, m, Q, K){
   b.hat 				<- array(0, dim=c(J, m, Q, K))
   for(j in 1:J)
	for(i in 1:m)
	   for(k in 1:K)
	   	b.hat[j,i,,k] 	<- D.hat[k]*t(data.z[i,,])%*%V.inv[j,i,,,k]%*%pResid[j,i,,k]

   return(b.hat)
}

### compute "b2.hat", b2.hat[j, k]
### if the variance matrix is a Diagonal matrix
compute.b2.hat.CLMM.simple 	<- function(data.z, u.hat, b.hat, D.hat, V.inv, J, m, Q, K){
	b2.hat 			<- array(0, dim=c(J, K))
	for(j in 1:J){
	   for(k in 1:K){
		temp1 		<- 0
	   	temp2 		<- 0
		tau2			<- D.hat[k]
	   	for(i in 1:m){
		   temp1 		<- temp1 + sum(b.hat[j,i,,k]^2)
		   temp2 		<- temp2 + sum(diag(t(data.z[i,,])%*%V.inv[j,i,,,k]%*%data.z[i,,]))
	   	}
	   	b2.hat[j,k]		<- temp1 + tau2*Q*m - tau2^2*temp2
	   }
	}
   return(b2.hat)
}

### compute "e.hat", e.hat[j, i, l, k]
compute.e.hat.CLMM.simple 	<- function(data.z, b.hat, pResid, J, m, L, K){
   e.hat 				<- array(0, dim=c(J,m,L,K))
   for(i in 1:m)
	for(j in 1:J){
	   e.hat[j,i,,] 		<- pResid[j,i,,] - as.matrix(data.z[i,,])%*%b.hat[j,i,,]
	}
   return(e.hat)
}

### compute "e2.hat", e2.hat[j, k]
compute.e2.hat.CLMM.simple 	<- function(data.z, u.hat, e.hat, V.inv, sigma2.hat, J, m, L, K){
	e2.hat 			<- array(0, dim=c(J, K))
	for(j in 1:J){
	   for(k in 1:K){
		temp1 		<- 0
	   	temp2 		<- 0
	   	for(i in 1:m){
		   temp1 		<- temp1 + sum(e.hat[j,i,,k]^2)
		   temp2 		<- temp2 + sum(diag(V.inv[j,i,,,k]))
	   	}
	   	e2.hat[j,k]		<- temp1 + sigma2.hat[k]*L*m - (sigma2.hat[k])^2*temp2
	   }
	}
   return(e2.hat)
}

### compute the log-transformed numerator for "u.hat"
compute.u.hat.CLMM.simple 	<- function(pResid, V.inv, pi.hat, J, m, K){
   # compute log(u.hat.numerator)
   log.u.hat.num 			<- matrix(0, nrow=J, ncol=K)
   for(k in 1:K){
	for(j in 1:J){
	   temp 			<- 0
	   for(i in 1:m){
		temp 			<- temp + mvn.dnorm.log(pResid[j,i,,k], V.inv[j,i,,,k])
	   }
	   log.u.hat.num[j,k] 	<- log(pi.hat[k]) + temp
	}
   }
   u.hat.num 			<- exp(t(apply(log.u.hat.num, MARGIN=1, FUN=all.ceiling)))
   u.hat.den 			<- apply(u.hat.num, FUN=sum, MARGIN=1)
   u.hat 				<- u.hat.num/u.hat.den

   return(u.hat)
}

### substract a constant from a vector to make its max = cutoff
all.ceiling 	<- function(aVector, cutoff=600){
    xx 		<- max(aVector)
    aVector 	<- aVector - xx + cutoff
    return(aVector)
}

### compute the density of a vector according to a MVN distribution
mvn.dnorm.log 	<- function(aVector, var.inv){
   L 			<- length(aVector)
   y 			<- matrix(aVector, nrow=L)
   # dens 		<- (abs(det(var.inv))^(1/2)) * ((2*pi)^(-L/2)) * exp(-(1/2) * t(y) %*% var.inv %*% y)

   log.dens 	<- log(abs(det(var.inv)))/2 + (log(2*pi)*(-L/2)) + ((-1/2) * t(y) %*% var.inv %*% y)

   return(log.dens)
}


### compute the "log likelihood" given this MLE
compute.llh.CLMM 		<- function(pResid, V.inv, pi.hat, J, m, K){
   llh 			<- 0
   for(j in 1:J){
	temp.log 		<- rep(0, K)
	for(k in 1:K){
	   temp.ind 	<- 0
	   for(i in 1:m){
	       temp.ind 	<- temp.ind + mvn.dnorm.log(pResid[j,i,,k], V.inv[j,i,,,k])
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
#		    2. D can ONLY be "Diagonal"	  	    		          #
#	    	    3. CLUSTER-SPECIFIC MEASUREMENT ERROR	    	    		    #
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

fit.CLMM.simple.diag.NA	<- function(data.y, data.x, data.z, na.flag, n.clst, n.start=1){
   temp			<- dim(data.y)
   if(length(temp)==2) {
	temp1			<- array(0, dim=c(temp[1], 1, temp[2]))
	temp1[,1,]		<- data.y
	data.y		<- temp1
   }
   temp			<- dim(data.x)
   if(length(temp)==2) {
	temp1			<- array(0, dim=c(1, temp))
	temp1[1,,]		<- data.x
	data.x		<- temp1
   }
   temp			<- dim(data.z)
   if(length(temp)==2) {
	temp1			<- array(0, dim=c(1, temp))
	temp1[1,,]		<- data.z
	data.z		<- temp1
   }

   ### this will be repeatedly used by M-step
   data.x.x.sum 		<- compute.x.z.sum.CLMM.simple.NA(data.x, data.x, na.flag)	# dim is JxPxP

   ### try different starting values
   llh 			<- -9999999999
   for(s in 1:n.start){
	### get "start values"
	theta.hat 		<-  fit.CLMM.simple.start.NA(data.x, data.y, data.z, n.clst, start=s)$theta.hat

	### iterate btw E- and M- steps
	est.hat.new 	<- fit.CLMM.simple.EM.NA(data.x, data.y, data.z, na.flag, data.x.x.sum, theta.hat)
	if(est.hat.new$theta.hat$llh > llh){
	   est.hat 		<- est.hat.new
	   llh		<- est.hat.new$theta.hat$llh
	}
   }
   return(est.hat)
}


### find a starting value for "zeta" in CLMM
library(cluster)
fit.CLMM.simple.start.NA 	<- function(data.x, data.y, data.z, n.clst, start=1){
   # number of genes and covariates
   J 		<- dim(data.y)[1]
   m 		<- dim(data.x)[1]
   P 		<- dim(data.x)[3]
   Q 		<- dim(data.z)[3]
   u.hat	<- matrix(0, nrow=J, ncol=n.clst)

   # get beta.hat for each gene-specific model
   beta.hat 		<- matrix(0, nrow=J, ncol=P)
   temp.x 			<- data.x[1,,]
   if(m>1){for(i in 2:m){
	temp.x 		<- rbind(temp.x, data.x[i,,])
   }}
   temp 			<- solve(t(temp.x) %*% temp.x) %*% t(temp.x)
   for(j in 1:J)
	beta.hat[j,] 	<- temp %*% as.vector(t(data.y[j,,]))

   # group gene-specific beta.hat's by PAM if "start==1"
   if(start==1) {
	temp 			<- pam(beta.hat, n.clst)
	# assign group centers to zeta.hat (a KxP matrix)
	zeta.hat 		<- temp$medoids
	temp			<- temp$clustering
   }

   # group gene-specific beta.hat's by PAM if "start==2"
   if(start==2) {
	temp 			<- kmeans(beta.hat, n.clst)
	# assign group centers to zeta.hat (a KxP matrix)
	zeta.hat 		<- temp$centers
 	temp			<- temp$cluster
  }

   # pick group centers randomly if "start>1"
   if(start>2) {
	temp 			<- sample(J, n.clst)
	zeta.hat 		<- beta.hat[temp,]
	temp			<- sample(n.clst, J)
   }

   for(k in 1:n.clst)
      u.hat[temp==k, k]	<- 1

  # covariance matrix for random effects
   D.hat 			<- rep(1, n.clst) 

   # measurement error
   sigma2.hat 		<- rep(1, n.clst)

   # frequency of each cluster
   pi.hat 			<- rep(1/n.clst, n.clst)

   return(list(theta.hat=list(zeta.hat=zeta.hat, D.hat=D.hat, sigma2.hat=sigma2.hat, pi.hat=pi.hat), u.hat=u.hat))
}


### EM algorithm to fit the CLMM
###
## data.x[i,l,p]
## data.y[j,i,l]
## data.z[i,l,q]
###
## zeta.hat, D.hat, sigma2.hat, and pi.hat are the starting values for the parameters
###

fit.CLMM.simple.EM.NA 	<- function(data.x, data.y, data.z, na.flag, data.x.x.sum, theta.hat){
   # number of genes, samples, covariates, and clusters
   J <- dim(data.y)[1]
   m <- dim(data.y)[2]
   L <- dim(data.y)[3]
   P <- dim(data.x)[3]
   Q <- dim(data.z)[3]
   K <- length(theta.hat$pi.hat)

   # "log likelihood"
   llh.old 			<- -9999999999
   llh 			<- -9999999990

   ### iterate btw E- and M- steps
   while(llh-llh.old>0.1){
	# E-step
	hats 			<- compute.Ehats.CLMM.simple.NA(data.x, data.y, data.z, na.flag, theta.hat, J, m, L, K, Q)

	# M-step
	temp 			<- compute.theta.hat.CLMM.simple.NA(data.x, data.y, data.z, na.flag, data.x.x.sum,hats,J,m,L,K,P,Q)

	# update only when llh increases; when some cluster disappears, llh decreases
	llh.old 		<- llh
	if(!is.na(temp$llh) & temp$llh > llh){
	   theta.hat 	<- temp
	   llh 		<- temp$llh
	   #print(llh)
	}
   } 
   return(list(u.hat=hats$u.hat, b.hat=hats$b.hat, theta.hat=theta.hat))
}


### compute "u.hat" - the expected clustering indicator
compute.Ehats.CLMM.simple.NA <- function(data.x, data.y, data.z, na.flag, theta.hat, J, m, L, K, Q){
   # delist "theta.hat"
   zeta.hat 	<- theta.hat$zeta.hat
   D.hat 		<- theta.hat$D.hat
   sigma2.hat 	<- theta.hat$sigma2.hat
   pi.hat 		<- theta.hat$pi.hat

   # compute V", V[j,i,l,l,k]
   V 			<- compute.V.simple.NA(data.z, D.hat, sigma2.hat, J, m, L)

   # compute the partial residuals "gene by gene" - pResid[j, i, k, l]
   pResid 		<- compute.pResid.CLMM.simple.NA(data.x, data.y, zeta.hat, J, m, L, K)

   # compute "u.hat"
   u.hat 		<- compute.u.hat.CLMM.simple.NA(pResid, V, pi.hat, J, m, K, na.flag)

   # compute "b.hat", bhat[j, i, k, q]
   b.hat 		<- compute.b.hat.CLMM.simple.NA(data.z, pResid, D.hat, V, J, m, Q, K, na.flag)

   # compute "b2.hat", b2.hat[j, k, q, q]
   b2.hat 		<- compute.b2.hat.CLMM.simple.NA(data.z, u.hat, b.hat, D.hat, V, J, m, Q, K, na.flag)

   # compute "e.hat", e.hat[j, i, k, l]
   e.hat 		<- compute.e.hat.CLMM.simple.NA(data.z, b.hat, pResid, J, m, L, K)

   # compute "e2.hat", e2.hat[j, k]
   e2.hat 		<- compute.e2.hat.CLMM.simple.NA(data.z, u.hat, e.hat, V, sigma2.hat, J, m, L, K, na.flag)
						
   return(list(u.hat=u.hat, b.hat=b.hat, b2.hat=b2.hat, e2.hat=e2.hat))
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
compute.theta.hat.CLMM.simple.NA <- function(data.x, data.y, data.z, na.flag, data.x.x.sum, hats, J, m, L, K, P, Q){
   # delist "hats"
   u.hat 				<- hats$u.hat
   b.hat 				<- hats$b.hat
   b2.hat 				<- hats$b2.hat
   e2.hat 				<- hats$e2.hat

   # frequency of each cluster
   pi.hat 				<- apply(u.hat, FUN=sum, MARGIN=2)/J

   # cluster-specific regression parameters
   zeta.hat 			<- matrix(0, nrow=K, ncol=P)
   D.hat				<- rep(0, K)
   sigma2.hat			<- rep(0, K)
   for(k in 1:K){
	zeta.num 			<- 0
	zeta.den			<- 0
	for(j in 1:J){
	   for(i in 1:m){
		temp1			<- !(na.flag[j,i,])
		temp2			<- t(data.x[i,temp1,])%*%(data.y[j,i,temp1]-as.matrix(data.z[i,temp1,])%*%b.hat[j,i,,k])
		zeta.num		<- zeta.num + u.hat[j,k]*temp2
	   }  
	   zeta.den 		<- zeta.den + u.hat[j,k]*data.x.x.sum[j,,]
	}
	zeta.hat[k,] 		<- solve(zeta.den) %*% zeta.num

	D.hat[k]			<- sum(u.hat[,k]*b2.hat[,k])/(m*Q*sum(u.hat[,k]))

	sigma2.hat[k] 		<- sum(u.hat[,k]*e2.hat[,k])/(m*L*sum(u.hat[,k]))
   }

   # compute V
   V 					<- compute.V.simple.NA(data.z, D.hat, sigma2.hat, J, m, L)

   # compute the "log likelihood" given this MLE
   pResid 				<- compute.pResid.CLMM.simple.NA(data.x, data.y, zeta.hat, J, m, L, K)
   llh 				<- compute.llh.CLMM.NA(pResid, V, pi.hat, J, m, K, na.flag)

   return(list(zeta.hat=zeta.hat, pi.hat=pi.hat, D.hat=D.hat, sigma2.hat=sigma2.hat, llh=llh))
}

### data.x[i,l,p], data.z[i,l,q]
compute.x.z.sum.CLMM.simple.NA 	<- function(data.x, data.z, na.flag){
   m <- dim(data.x)[1]
   P <- dim(data.x)[3]
   Q <- dim(data.z)[3]
   J <- dim(na.flag)[1]

   data.x.z.sum 	<- array(0, dim=c(J,P,Q))
   for(j in 1:J){
	for(i in 1:m){
	   temp			<- !(na.flag[j,i,])
	   data.x.z.sum[j,,] 	<- data.x.z.sum[j,,] + t(data.x[i,temp,]) %*% data.z[i,temp,]
	}
   }
   return(data.x.z.sum)
}

### compute V[j, i, l, l, k]
compute.V.simple.NA 	<- function(data.z, D.hat, sigma2.hat, J, m, L){
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
compute.pResid.CLMM.simple.NA 	<- function(data.x, data.y, zeta.hat, J, m, L, K){
   pResid 				<- array(0, dim=c(J, m, L, K))
   for(k in 1:K){
	for(i in 1:m){
	   y.hat 			<- data.x[i,,] %*% zeta.hat[k,] 		# L x 1
	   pResid[,i,,k]		<- t(t(data.y[,i,]) - as.vector(y.hat))	# J x L
	}
   }
   return(pResid)
}

### compute "b.hat", bhat[j, i, q, k]
compute.b.hat.CLMM.simple.NA 	<- function(data.z, pResid, D.hat, V, J, m, Q, K, na.flag){
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
compute.b2.hat.CLMM.simple.NA 	<- function(data.z, u.hat, b.hat, D.hat, V, J, m, Q, K, na.flag){
	b2.hat 			<- array(0, dim=c(J, K))
	for(j in 1:J){
	   for(k in 1:K){
		temp1 		<- 0
	   	temp2 		<- 0
		tau2			<- D.hat[k]
	   	for(i in 1:m){
		   temp		<- !(na.flag[j,i,])
		   temp1 		<- temp1 + sum(b.hat[j,i,,k]^2)
		   temp2 		<-temp2+sum(diag(t(data.z[i,temp,])%*%solve(V[j,i,temp,temp,k],tol=1e-50)%*%data.z[i,temp,]))
	   	}
	   	b2.hat[j,k]		<- temp1 + tau2*Q*m - (tau2^2)*temp2
	   }
	}
   	return(b2.hat)
}

### compute "e.hat", e.hat[j, i, l, k]
compute.e.hat.CLMM.simple.NA 	<- function(data.z, b.hat, pResid, J, m, L, K){
   e.hat 				<- array(0, dim=c(J,m,L,K))
   for(j in 1:J)
	for(i in 1:m){
	   e.hat[j,i,,] 		<- pResid[j,i,,] - as.matrix(data.z[i,,])%*%b.hat[j,i,,]
	}
   return(e.hat)
}

### compute "e2.hat", e2.hat[j, k]
compute.e2.hat.CLMM.simple.NA 	<- function(data.z, u.hat, e.hat, V, sigma2.hat, J, m, L, K, na.flag){
	e2.hat 			<- array(0, dim=c(J, K))
	for(j in 1:J){
	   for(k in 1:K){
		temp1 		<- 0
	   	temp2 		<- 0
	   	for(i in 1:m){
		   temp		<- !(na.flag[j,i,])
		   temp1 		<- temp1 + sum(e.hat[j,i,temp,k]^2)
		   temp2 		<- temp2 + sum(diag(solve(V[j,i,temp,temp,k],tol=1e-50)))
	   	}
	   	e2.hat[j,k]		<- temp1 + sigma2.hat[k]*(L*m-sum(na.flag[j,,])) - (sigma2.hat[k])^2*temp2
	   }
	}
   	return(e2.hat)
}

### compute the log-transformed numerator for "u.hat"
compute.u.hat.CLMM.simple.NA 	<- function(pResid, V, pi.hat, J, m, K, na.flag){
   # compute log(u.hat.numerator)
   log.u.hat.num 			<- matrix(0, nrow=J, ncol=K)
   for(k in 1:K){
	for(j in 1:J){
	   temp1 			<- 0
	   for(i in 1:m){
		temp			<- !(na.flag[j,i,])
		temp1			<- temp1 + mvn.dnorm.log.NA(pResid[j,i,temp,k], solve(V[j,i,temp,temp,k],tol=1e-50))
	   }
	   log.u.hat.num[j,k] 	<- log(pi.hat[k]) + temp1
	}
   }
   u.hat.num 			<- exp(t(apply(log.u.hat.num, MARGIN=1, FUN=all.ceiling.NA)))
   u.hat.den 			<- apply(u.hat.num, FUN=sum, MARGIN=1)
   u.hat 				<- u.hat.num/u.hat.den

   return(u.hat)
}

### substract a constant from a vector to make its max = cutoff
all.ceiling.NA 	<- function(aVector, cutoff=600){
    xx 		<- max(aVector)
    aVector 	<- aVector - xx + cutoff
    return(aVector)
}

### compute the density of a vector according to a MVN distribution
mvn.dnorm.log.NA 	<- function(aVector, var.inv){
   L 			<- length(aVector)
   y 			<- matrix(aVector, nrow=L)
   # dens 		<- (abs(det(var.inv))^(1/2)) * ((2*pi)^(-L/2)) * exp(-(1/2) * t(y) %*% var.inv %*% y)

   log.dens 	<- log(abs(det(var.inv)))/2 + (log(2*pi)*(-L/2)) + ((-1/2) * t(y) %*% var.inv %*% y)

   return(log.dens)
}


### compute the "log likelihood" given this MLE
compute.llh.CLMM.NA 		<- function(pResid, V, pi.hat, J, m, K, na.flag){
   llh 			<- 0
   for(j in 1:J){
	temp.log 		<- rep(0, K)
	for(k in 1:K){
	   temp.ind 	<- 0
	   for(i in 1:m){
		temp		<- !(na.flag[j,i,])
		temp.ind	<- temp.ind + mvn.dnorm.log.NA(pResid[j,i,temp,k], solve(V[j,i,temp,temp,k],tol=1e-50))
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
