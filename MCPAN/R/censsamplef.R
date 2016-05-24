`censsamplef` <-
function(n, scale.m, shape.m, scale.t, shape.t=3, tmax )

{
k <- length(n)

if( any(c(length(scale.m), length(shape.m), length(scale.t), length(shape.t))!=k))
 {stop("n, scale.m, shape.m, scale.t, shape.t must be vectors of equal length")}

names <- as.factor(rep( 1:k , times=c(n)))


dat<-data.frame()

for(i in 1:k)
{
 dat <- rbind(dat, 
 censsample(n=n[i], scale.m=scale.m[i], shape.m=shape.m[i], scale.t=scale.t[i], shape.t=shape.t[i], tmax=tmax ) )
}
dat$f <- names

return(dat)

}

