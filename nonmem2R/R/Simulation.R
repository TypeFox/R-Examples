#' Combine fix grid and simulated grid based on multivariate normal distribution
#' @description
#' grid.sim produce comparable output to mvnorm however with one column having a fixed range of values rather than a
#' simulated range of values. As the number of simulations (n) goes to infinity grid.sim is identical to mvnorm in
#' that the covariance of the output's both will converge to the input covariance matrix sigma.
#' The advantage with grid.sim is that the function provides more stable results, hence number of simulations can
#' be reduced and still have equally stable results when used to represent parameter uncertainty and or population variability
#' in model predictions.
#' @param n
#' Number of simulations
#' @param means
#' vector of mean values
#' @param sigma
#' covariance matrix
#' @param grid.param
#' the index of the parameter for which a fix grid (from qnorm) is used instead of a simulated grid (rnorm).
#' If grid.param=NULL ( default) the fix grid will be used for the parameter with largest variance.
#' @param pure.sim
#' do pure simulation for all parameters, default=FALSE
#' @return
#' row-matrix of parameters
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats cov
#' @importFrom stats qnorm
#' @importFrom stats rnorm
#' @importFrom stats sd
#' @examples
#' sigma<-matrix(c(1,0.5,0.5,2),ncol=2)
#' sim1<-grid.sim(1000,sigma=sigma)
#' pairs(sim1)
#' cov(sim1)
grid.sim<-function(n,means=NULL,sigma,grid.param=NULL,pure.sim=FALSE) {


  if(!is.matrix(sigma)){
		sigma<-matrix(sigma,ncol=floor(sqrt(length(sigma))))
	}

	if(length(means)==0) {
		means=rep(0,ncol(sigma))
	}

	res<-matrix(means,byrow=T,ncol=length(means),nrow=n)

	## Check for any variance==0
	zero.var<-(diag(sigma)==0)

	## Covariance and mean of those with variance>0
	sigma2<-sigma[!zero.var,!zero.var]
	means2<-means[!zero.var]


	## for grid & simulation combination, select parameter for grid by maximum variance
	if(is.null(grid.param)){
		grid.param<-which.max(diag(sigma2))
	}

	X<-NULL
	if(length(means2)>1){
		X<-rmvnorm(n,mean=means2,sigma=sigma2)
	}
	if(length(means2)==1){
		X<-matrix(rnorm(n,mean=means2,sd=sqrt(sigma2)))
  }

	if(!pure.sim & length(means2)>1){
		## Formulas for conditional MV from
		## https://en.wikipedia.org/wiki/Multivariate_normal_distribution
		X2<-qnorm(seq(0.5/n,1-0.5/n,length=n))
		X2<-X2*sqrt(diag(sigma2)[grid.param])/sd(X2)+means2[grid.param]

		mu1<-means2[-grid.param]
	    mu2<-means2[ grid.param]
		SIGMA12<-sigma2[-grid.param, grid.param]
		SIGMA22<-sigma2[ grid.param, grid.param]
		SIGMA11<-sigma2[-grid.param,-grid.param]
		mubar<-NULL
		for(i in 1:length(mu1)){
			mubar<-cbind(mubar,mu1[i]+SIGMA12[i]*(X2-mu2)/SIGMA22)
		}
		SIGMAbar<-SIGMA11-SIGMA12%*%t(SIGMA12)/SIGMA22
		X1<-rmvnorm(n,sigma=SIGMAbar)+mubar
		## Put X1 and X2 togheter
		#X<-cbind(X1,X2)
		X[, grid.param]<-X2
		X[,-grid.param]<-X1
	}
	if(!pure.sim & length(means2)==1){
		X<-qnorm(seq(0.5/n,1-0.5/n,length=n))
		X<-matrix(X*sqrt(sigma2)/sd(X) + means2)
	}

	res[,!zero.var]<-X
	res

}

#' Function for testing grid.sim and compare with mvnorm
#' @description
#' Test grid.sim
#' @param n
#' number of simulations
#' @param k
#' subset of parameters from a 4X4 sigma to use
#' @return
#' grapics
#' @export
#' @importFrom lattice xyplot
#' @examples
#' \dontrun{
#' require(lattice)
#' test.grid.sim(n=1000)
#' }
test.grid.sim<-function(n=1000,k=1:4){
	sigma<-matrix(ncol=4,byrow=T,c(
		0,0  ,0 ,0,
		0,1  ,0.5,0,
		0,0.5,1.0,0.25,
		0,0  ,0.25,2
	))
	means<-0:3
	sigma<-sigma[k,k]
	means<-means[k]
	a<-grid.sim(n,means,sigma)
	if(length(means)>1){
		b<-rmvnorm(n,mean=means,sigma=sigma)
	}
	else{
		b<-matrix(rep(means,length=n))
		if(sigma>0){
			b<-matrix(rnorm(n,mean=means,sd=sqrt(sigma)))
		}
	}

	#b<-rmvnorm(n,means,sigma)
	r<-NULL
	q<-qnorm(seq(0.5/n,1-0.5/n,length=n))
	for(i in 1:length(means)){
		r<-rbind(r,data.frame(i=k[i],q=q,grid.sim=sort(a[,i]),mvnorm=sort(b[,i])))
	}
	p1<-xyplot(grid.sim+mvnorm~q|factor(i),cex=c(0.8,0.6),data=r,as.table=T,auto.key=list(columns=2))
	cat("Covariance from grid.sim output:\n")
	print(cov(a))
	cat("Covariance from mvnorm output:\n")
	print(cov(b))
	p1
}



#test.grid.sim(n=1000,k=1:4)
#test.grid.sim(n=1000,k=1:3)
#test.grid.sim(n=1000,k=2:4)
#test.grid.sim(n=1000,k=2:3)
#test.grid.sim(n=1000,k=1:2)
#test.grid.sim(n=1000,k=2)
#test.grid.sim(n=1000,k=1)



