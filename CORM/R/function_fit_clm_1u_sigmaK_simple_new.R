###################################################################
#####							      		  #####
###										    ###
#	      fit a CLM model with ONE common clustering	      #
#	       	  and SAMPLE-specific covariates x_i	      #
#								  			#
#	    ALLOW FOR CLUSTER-SPECIFIC MEASUREMENT ERROR	  	#
#								  			#
#			simplified M-step for zeta.hat  			#
###										    ###
#####							     			  #####
###################################################################

##### INPUT
###
##	data.y - matrix of observations,
##		  data.y[j, i] for sample i and gene j;
###
## 	data.x - matrix of covariates,
##		  data.x[i, p] for sample i, gene j, and covariate p,
###
##	data.z - NULL (place holder)
###
##	n.clst - number of clusters for the beta associated with data.x;
###
##	type.x - type for data.x1, where it takes value
##			"sample" for sample-specific covariates,
##		  	"sample-gene" for sample-gene-specific covariates;
###
##### OUTPUT (a list of)
###
##	theta.hat - regression parameters estimated via EM algorithm;
###
##	data.u - clustering associated with data.x
###

#The name of main function is changed from "fit.CLM.1u.sigmaK.simple" to "fit.CLM".

fit.CLM 	<- function(data.y, data.x, n.clst, n.start=1){
   ### this will be repeatedly used by M-step
   data.x.x 		<- t(data.x) %*% data.x 	# dim is PxP

   ### try different starting values
   llh 			<- -9999999999
   for(s in 1:n.start){
	### get "start values"
	theta.hat 		<-  fit.CLM.1u.simple.start(data.x, data.y, K=n.clst, start=s)$theta.hat

	### iterate btw E- and M- steps
	est.hat.new 	<- fit.CLM.1u.simple.EM(data.x, data.y, data.x.x, theta.hat)
	if(est.hat.new$theta.hat$llh > llh){
	    est.hat 	<- est.hat.new
	    llh		<- est.hat.new$theta.hat$llh
	}
	print(llh)
   }
   return(est.hat)
}


### find a starting value for "zeta" in CLM
library(cluster)
fit.CLM.1u.simple.start <- function(data.x, data.y, K, start){
   # number of genes and covariates
   J 				<- nrow(data.y)
   P 				<- ncol(data.x)
   u.hat			<- matrix(0, nrow=J, ncol=K)

   # K-means if "start<=4"
   if(start<=4) {
	temp 		<- kmeans(data.y, centers=K)
	zeta.hat	<- temp$centers %*% data.x %*% t(solve(t(data.x) %*% data.x)) 	# KxP matrix
	temp		<- temp$cluster
   }

   # PAM if "start<=8"
   if(start>4 & start<=8) {
	temp 		<- pam(data.y, K)
	zeta.hat	<- temp$medoids %*% data.x %*% t(solve(t(data.x) %*% data.x))
	temp		<- temp$clustering
   }

  # group gene-specific beta.hat's by K-means if "start==9"
   if(start==9) {
   	beta.hat 	<- data.y %*% data.x %*% t(solve(t(data.x) %*% data.x))		# JxP matrix
	temp 		<- kmeans(beta.hat, K)
	zeta.hat	<- temp$centers
	temp		<- temp$cluster
   }

   # pick group centers randomly if "start>11"
   if(start>9) {
	temp 		<- sample(J, K)
   	beta.hat 	<- data.y %*% data.x %*% t(solve(t(data.x) %*% data.x))		# JxP matrix
	zeta.hat	<- beta.hat[temp,]
	temp		<- sample(K, J, replace=TRUE)
   }

   for(k in 1:K)   
	u.hat[temp==k, k]	<- 1

   # measurement error
   sigma2.hat 		<- rep(10, K)

   # frequency of each cluster
   pi.hat 			<- rep(1/K, K)

   return(list(theta.hat=list(zeta.hat=zeta.hat, sigma2.hat=sigma2.hat, pi.hat=pi.hat), u.hat=u.hat))
}


### EM algorithm to fit the CLM
###
## data.x[i,p]
## data.y[j,i]
###
## zeta.hat, sigma2.hat, and pi.hat are the starting values for the parameters
###

fit.CLM.1u.simple.EM 	<- function(data.x, data.y, data.x.x, theta.hat){
   # number of genes, samples, covariates, and clusters
   J <- nrow(data.y) 
   N <- ncol(data.y)
   P <- ncol(data.x)
   K <- length(theta.hat$pi.hat)

   # "log likelihood"
   llh.old 		<- -9999999999
   llh 		<- -9999999990

   ### iterate btw E- and M- steps
   while(llh-llh.old>0.01){
	# E-step
	u.hat 	<- compute.u.hat.simple(data.x, data.y, theta.hat, J, N, K)

	# M-step
	temp 		<- compute.theta.hat.simple(data.x, data.y, data.x.x, u.hat, J, N, K, P)

	# update only when llh increases; when some cluster disappears, llh decreases
	llh.old 	<- llh
	if(temp$llh > llh){
	   theta.hat<- temp
	   llh 	<- temp$llh
	}
   } 
   return(list(u.hat=u.hat, theta.hat=theta.hat))
}


### compute "u.hat" - the expected clustering indicator
compute.u.hat.simple 	<- function(data.x, data.y, theta.hat, J, N, K){
   # delist "theta.hat"
   zeta.hat 		<- theta.hat$zeta.hat
   sigma2.hat 		<- theta.hat$sigma2.hat
   pi.hat 			<- theta.hat$pi.hat

   # compute the residuals "gene by gene"
   residuals 		<- compute.residuals.simple(data.x, data.y, zeta.hat, J, N, K)

   # compute the numerator for "u.hat"
   u.hat.num 		<- matrix(0, nrow=J, ncol=K)
   for(k in 1:K){
	temp 			<- dnorm(residuals[,,k], sd=sqrt(sigma2.hat[k]))
	u.hat.num[,k] 	<- pi.hat[k] * apply(temp, FUN=prod, MARGIN=1)
   }

   # compute the denominator for "u.hat"
   u.hat.den 		<- apply(u.hat.num, FUN=sum, MARGIN=1)

   return(u.hat.num/u.hat.den)
}


### M-step for fitting CLM
###
### INPUT
##
#   data.x.x[j,,] = t(data.x[j,,])%*%(data.x[j,,])
#   data.x.y[j,] = t(data.x[j,,])%*%(data.y[j,])
##
#   u.hat is a J*K matrix of cluster membership probabilities
##
###
compute.theta.hat.simple <- function(data.x, data.y, data.x.x, u.hat, J, N, K, P){
   # estimtate "zeta" for each cluster
   zeta.hat 		<- matrix(0, nrow=K, ncol=P)
   for(k in 1:K){
	zeta.hat.num 	<- apply(u.hat[,k]*data.y, FUN=sum, MARGIN=2)
	zeta.hat.den 	<- sum(u.hat[,k]) * data.x.x
	zeta.hat[k,] 	<- solve(zeta.hat.den) %*% (t(data.x) %*% zeta.hat.num)
   }

   # estimate the measurement error
   sigma2.hat 		<- rep(0, K)
   residuals 		<- compute.residuals.simple(data.x, data.y, zeta.hat, J, N, K)
   for(k in 1:K){
	res.sum 		<- apply(residuals[,,k]^2, FUN=sum, MARGIN=1)%*%u.hat[,k]	# numerator
	res.den		<- N*sum(u.hat[,k])							# denominator
	sigma2.hat[k]	<- res.sum/res.den
   }
   
   # frequency of each cluster
   pi.hat 			<- apply(u.hat, FUN=sum, MARGIN=2)/J

   # compute the "log likelihood" given this MLE
   llh 			<- compute.llh(residuals, sigma2.hat, pi.hat, J, K)

   return(list(zeta.hat=zeta.hat, pi.hat=pi.hat, sigma2.hat=sigma2.hat, llh=llh))
}


### compute "residuals"
compute.residuals.simple <- function(data.x, data.y, zeta.hat, J, N, K){
   # compute the fitted values
   y.hat 			<- data.x %*% t(zeta.hat)

   # compute the residuals
   residuals 		<- array(0, dim=c(J,N,K))
   for(k in 1:K) {
	residuals[,,k] 	<- t(t(data.y) - y.hat[,k])
   }
   return(residuals)
}


### compute the "log likelihood" given this MLE
compute.llh 	<- function(residuals, sigma2.hat, pi.hat, J, K){
   llh 		<- 0
   for(j in 1:J){
	temp 		<- 0
	for(k in 1:K){
	   temp 	<- temp + pi.hat[k]*prod(dnorm(residuals[j,,k], sd=sqrt(sigma2.hat[k])))
	}
	llh 		<- llh + log(temp)
   }
   return(llh)
}


############			the end 			##############
