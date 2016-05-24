observation.obj.func <- function(residual,sample.step.res,initial.subsample.size){
N <- dim(residual)[1]
obj <- apply(residual^2,c(1,3),sum)
min.exclude <- integer(0)
max.include <- integer(0)
m.obj <- integer(0)
mplus1.obj <- integer(0)
for(i in initial.subsample.size:(N-1)){
  sub.set <- sample.step.res[which(sample.step.res[,i]>0),i]
  max.include[i] <- sort(obj[sub.set,i])[i]
  min.exclude[i] <- sort(obj[-sub.set,i])[1]
  m.obj[i] <- sort(obj[,i])[i]
  mplus1.obj[i] <- sort(obj[,i])[i+1]
  }
obj.list <- list(observation.obj=obj, min.excl=min.exclude, max.incl=max.include, m.obj=m.obj, mplus1.obj=mplus1.obj)
class(obj.list) <- "fs.class"
return(obj.list)
}
