# calculates the skew of a distribution for a vector of observations
# (not counts or densities)

skew<-function(x) {
 # get rid of NAs, worry about the effects later
 x<-x[!is.na(x)]
 nx<-length(x)
 meanx<-mean(x)
 devx<-x-meanx
 sampskew<-(sum(devx^3)/nx)/(sum(devx^2)/nx)^1.5
 popskew<-sampskew*sqrt(nx*(nx-1))/(nx-2)
 return(list(sample=sampskew,population=popskew))
}
