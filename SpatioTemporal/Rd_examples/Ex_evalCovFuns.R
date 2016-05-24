##vector of distances
d <- seq(0,10,length.out=1e4);
##just the simplest case (exponential, range=2, sill=0.7)
plot(d, evalCovFuns("exp", c(2,0.7), d), type="l")
     
##create list of ranges
range <- c(1, 2, 3.5, 5);
##list names
name <- list("exp", "exp2", "cubic", "spherical", "cauchy", "cauchy",
             "matern", "matern")
##and list of shapes
shape <- c(vector("list",4), list(1, 5, .25, 5))

##matrix holding results
covf <- array(NA,c(length(d),length(name),length(range)))

##compute a few covariance functions
for(i in 1:length(name)){
  for(j in 1:length(range)){
    pars <- c(range[j],1,shape[[i]])
    covf[,i,j] <- evalCovFuns(name[[i]], pars, d)
  }
}

##plot the covariance function for comparison
par(mfrow=c(2,2))
for(j in 1:length(range)){
  plot(0, 0, type="n", main=range[j],
       xlim=range(d), ylim=range(covf[,,j],na.rm=TRUE))
  for(i in 1:length(name)){
    lines(d, covf[,i,j], col=i)
  }
  abline(v=range[j])
  if(j==1){
    legend("topright", lty=1, col=1:length(name),
           legend=paste("covf:", sapply(name,as.character),
			sapply(shape, function(x){
				      if(is.null(x)){""}else{as.character(x)} })))
  }
}
