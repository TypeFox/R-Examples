EmpiricalSemivariogram <-
function(ssn.object, varName,
	nlag = 20, directions = c(0,45,90,135),
	tolerance = 22.5, inc = 0, maxlag = 1e32, nlagcutoff = 1,
	EmpVarMeth = "MethMoment")
# EMPIRICAL SEMIVARIOGRAM FUNCTION
# var1 is a matrix or data frame with x-coord in the first column
#                                     y-coord in the second column
#                                     z (response) in the third column
{

  if(class(ssn.object) == "influenceSSN") ssn.object <- ssn.object$ssn.object
	
  #if(class(ssn.object) == "SpatialStreamNetwork") {
    data <- ssn.object@obspoints@SSNPoints[[1]]@point.data
    var <- varName
    x <- ssn.object@obspoints@SSNPoints[[1]]@point.coords[,1]
    y <- ssn.object@obspoints@SSNPoints[[1]]@point.coords[,2]
  #}
	n1 <- length(data[,1])
   # distance matrix among locations
	distance <- sqrt( ( matrix(x,nrow=n1,ncol=1) %*%
		matrix(rep(1,times=n1),nrow=1,ncol=n1) -
		matrix(rep(1,times=n1),nrow=n1,ncol=1) %*%
		matrix(x,nrow=1,ncol=n1) )^2 +
		( matrix(y,nrow=n1,ncol=1) %*%
		matrix(rep(1,times=n1),nrow=1,ncol=n1) -
		matrix(rep(1,times=n1),nrow=n1,ncol=1) %*%
		matrix(y,nrow=1,ncol=n1) )^2 )
	difx <- -(matrix(y,nrow=n1,ncol=1) %*%
		matrix(rep(1,times=n1),nrow=1,ncol=n1) -
		matrix(rep(1,times=n1),nrow=n1,ncol=1) %*%
		matrix(y,nrow=1,ncol=n1))
	signind <- -(matrix(x,nrow=n1,ncol=1) %*%
		matrix(rep(1,times=n1),nrow=1,ncol=n1) -
		matrix(rep(1,times=n1),nrow=n1,ncol=1) %*%
		matrix(x,nrow=1,ncol=n1)) < 0
	distance <- distance*1.0000000001
	theta.deg <- acos(difx/distance)*180/pi
   # matrix of degrees clockwise from north between locations
	theta.deg[signind] <- 360-theta.deg[signind]
	diff2 <- ( matrix(data[,var],nrow=n1,ncol=1) %*%
		matrix(rep(1,times=n1),nrow=1,ncol=n1) -
		matrix(rep(1,times=n1),nrow=n1,ncol=1) %*%
		matrix(data[,var],nrow=1,ncol=n1) )^2
	sqrtdiff <- sqrt(abs( matrix(data[,var],nrow=n1,ncol=1) %*%
		matrix(rep(1,times=n1),nrow=1,ncol=n1) -
		matrix(rep(1,times=n1),nrow=n1,ncol=1) %*%
		matrix(data[,var],nrow=1,ncol=n1) ) )
	if(EmpVarMeth == "CovMean") temp4cov <- data[,var] - mean(data[,var])
	else temp4cov <- data[,var]
	covprod <- (matrix(temp4cov,nrow=n1,ncol=1) %*%
		matrix(rep(1,times=n1),nrow=1,ncol=n1)) *
		(matrix(rep(1,times=n1),nrow=n1,ncol=1) %*%
		matrix(temp4cov,ncol=n1,nrow=1))
# convert to vectors
	distance <- matrix(distance, ncol = 1)
	theta.deg <- matrix(theta.deg, ncol = 1)
	diff2 <- matrix(diff2, ncol = 1)
	sqrtdiff <- matrix(sqrtdiff, ncol = 1)
	covprod <- matrix(covprod, ncol = 1)
# trim off values greater than maxlag
	indmax <- distance <= maxlag
	distance <- distance[indmax,]
	theta.deg <- theta.deg[indmax,]
	diff2 <- diff2[indmax,]
	sqrtdiff <- sqrtdiff[indmax,]
	covprod <- covprod[indmax,]

	maxd<-max(distance)
	if( inc <= 0) inc <- maxd/nlag
	ind <- distance==0
	ndir <- length(directions)
	store.results <- matrix(data = NA, ncol = 6,
		dimnames = list(NULL, c("distance", "gamma", "np", "azimuth", "hx", "hy")))
	for (j in 1:ndir) {
		for ( i in 1:nlag){
			if( (directions[j]-tolerance)<0 && (directions[j]+tolerance)>0 )
				ind1 <- theta.deg >= 360+directions[j]-tolerance |
					theta.deg < directions[j]+tolerance
			else if( (directions[j]+tolerance)>360 && (directions[j]-tolerance)<360 )
				ind1 <- theta.deg < directions[j]+tolerance-360 |
					theta.deg >= directions[j]-tolerance
			else
				ind1 <- theta.deg >= directions[j]-tolerance &
					theta.deg < directions[j]+tolerance
			ind<-distance>(i-1)*inc & distance<=i*inc &
				!is.na(theta.deg) & ind1
			nclass <- sum(ind)
			if(EmpVarMeth == "MethMoment") cv <- mean(diff2[ind])
			if(EmpVarMeth == "RobustMean") cv <- ((mean(sqrtdiff[ind]))^4)/(.457 + .494/sum(ind))
			if(EmpVarMeth == "RobustMedian") cv <- (median(sqrtdiff[ind]))^4/.457
			if(EmpVarMeth == "Covariance" | EmpVarMeth == "CovMean") cv <- mean(covprod[ind])
			mean.dis <- mean(distance[ind])
			if(nclass > 0) store.results<-rbind(store.results,
				c(mean.dis,cv,nclass,directions[j],0,0))
		}
	}
	store.results[,"hx"]<-store.results[,"distance"]*sin(store.results[,"azimuth"]*pi/180)
	store.results[,"hy"]<-store.results[,"distance"]*cos(store.results[,"azimuth"]*pi/180)
	store.results[,"gamma"]<-store.results[,"gamma"]/2
	ind <- store.results[,"np"] >= nlagcutoff
	store.results <- store.results[ind,]
	ind <- !is.na(store.results[,"hx"])
	store.results <- store.results[ind,]
	store.results <- as.data.frame(store.results)
	class(store.results) <- "EmpiricalSemivariogramSSN"
	store.results
}

