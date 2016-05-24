`Stem.Kriging` <-
function(StemModel,coord.newlocations,covariates.newlocations,
					K.newlocations,time.point,cov.spat=Sigmastar.exp){

z 		= StemModel$data$z
p 		= StemModel$skeleton$p
n 		= StemModel$data$n
d 		= StemModel$data$d
r 		= StemModel$data$r
K		= StemModel$skeleton$K

covariates 	= StemModel$data$covariates
covariates  = changedimension_covariates(covariates,d=d,r=r,n=n)

if (length(StemModel$estimates$phi.hat) > 0 ){ 
  phi = StemModel$estimates$phi.hat
  }
else { 
  (phi = StemModel$skeleton$phi)
  }
   
phi$logb 	= log(phi$sigma2eps/phi$sigma2omega)
phi$logtheta= log(phi$theta)

y.smoothed	= StemModel$estimates$y.smoothed

dist    	= as.matrix(dist(StemModel$data$coordinates,diag=TRUE)) #distance matrix between original locations
colnames(coord.newlocations) = colnames(StemModel$data$coordinates) 

m = nrow(coord.newlocations)

if(time.point<1 | time.point>StemModel$data$n) stop(paste("time.point must be between 1 and",StemModel$data$n))
if(!ncol(coord.newlocations==2)) stop("coord.newlocations must have 2 columns")
if(!ncol(covariates.newlocations==r)) stop(paste("covariates.newlocations must have",StemModel$data$r,"columns"))
if(!(m==nrow(covariates.newlocations))) stop("coord.newlocations and covariates.newlocations must have the same number or rows")


new.distancematrix = matrix(NA, m, d)
for(s in 1:m){
		new.distancematrix[s,] = as.matrix(dist(rbind(StemModel$data$coordinates,coord.newlocations[s,]),diag=TRUE))[(d+1),-(d+1)]
	}
	
multi.pred = spatial.pred(mu1 = covariates[,,time.point]%*%phi$beta + K%*%y.smoothed[time.point,],
				mu2	= covariates.newlocations %*% phi$beta + K.newlocations %*% y.smoothed[time.point,],
				Sigma11 = phi$sigma2omega * cov.spat(d=d , logb=phi$logb , logtheta=phi$logtheta , dist=dist),
				Sigma12 = phi$sigma2omega * exp(-phi$theta * new.distancematrix),
				Sigma22 = diag(phi$sigma2omega,m),
				X1	= z[time.point,]
		)	
	

data.newlocations =  list(coordinates = coord.newlocations, covariates=covariates.newlocations,
					K = K.newlocations,z=multi.pred$pred,se.pred=multi.pred$se.pred)

return(list(data.newlocations=data.newlocations,time.point=time.point))

}

