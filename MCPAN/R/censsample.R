`censsample` <-
function(n, scale.m, shape.m, scale.t, shape.t=3, tmax )
{

if( any(c(length(n), length(scale.m), length(shape.m), length(scale.t), length(shape.t))!=1) )
 {stop("input single numbers")}

# T.m time of mortality
# T.t time of tumor onset

T.m <- rweibull(n=n, shape=shape.m, scale=scale.m)      ###### Zeit bis Tod
T.t <- rweibull(n=n, shape=shape.t, scale=scale.t)      ###### Zeit bis Tumor = Tod

# at time of death of an animal dies: is a tumor present or not:
# time at which animals die or are sacrificed min(T.t, T.m)
# status=TRUE if before alternative death or final sacrifice a tumor developed
# status is FALSE if time of tumor onset T.t is greater than tmax or T.m

mat<-cbind(T.t, T.m, tmax)
time <- apply(X=mat, MARGIN=1, FUN=function(x){min(x[2:3])})
status <- apply( X=mat, MARGIN=1, FUN=function(x){ x[1] < min(x[2:3]) } )
return(data.frame(time=time, status=status, T.t=T.t, T.m=T.m, tmax = tmax))

}

