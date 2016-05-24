##load data
data(mesa.model)

##Mark 30% of observations
I <- runif(dim(mesa.model$obs)[1])<.3
##drop these observations
mesa.model.new <- dropObservations(mesa.model, I)

##This reduces the remaining number of observations
print(mesa.model)
print(mesa.model.new)

\dontshow{
  if( (length(mesa.model$obs$obs)<=length(mesa.model.new$obs$obs)) ||
     (length(mesa.model.new$obs$obs)!=sum(!I)) ){
    stop("drop.observations 1: Observations not dropped")
  }
}
##create cross validation structure
Icv <- createCV(mesa.model, groups=10)

##drop observations from the second CV group
mesa.model.new <- dropObservations(mesa.model, Icv==2)

##This reduces the remaining number of observations (and locations)
print(mesa.model)
print(mesa.model.new)

\dontshow{
  if( (length(mesa.model$obs$obs)<=length(mesa.model.new$obs$obs)) ||
     (length(mesa.model.new$obs$obs)!=sum(Icv!=2)) ||
     (dim(mesa.model$locations)[1] <= dim(mesa.model.new$locations)[1]) ){
    stop("drop.observations 2: Observations not dropped")
  }
}

