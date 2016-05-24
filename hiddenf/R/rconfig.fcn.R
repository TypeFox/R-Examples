rconfig.fcn <-
function(rows)
{
avalues <- as.vector(names(table(rows)))
a <- length(table(rows))
if(a > 20){stop("package not yet ready for a>20")}
cc <- 2^(a-1)-1
rconfig.mtx <- matrix(NA,nrow=length(rows),ncol=cc)
counter <- 1
#if(is.even(b))
if((a %% 2)==0)
{
for(g in 1:((a/2)-1))
{
ulim <- choose(a,g)
for(j in 1:ulim)
{
combo <- combn(a,g)
rconfig.mtx[,counter] <- 1*is.element(rows,combo[,j])
counter <- counter+1
}
}
g <- g+1
ulim <- choose(a,g)/2
for(j in 1:ulim)
{
combo <- combn(a,g)[,1:ulim]
rconfig.mtx[,counter] <- 1*is.element(rows,combo[,j])
counter <- counter+1
}
}
#else if (is.odd(b))
else if((a %% 2)==1)
{
for(g in 1:((a-1)/2))
{
ulim <- choose(a,g)
for(j in 1:ulim)
{
combo <- combn(a,g)
rconfig.mtx[,counter] <- 1*is.element(rows,combo[,j])
counter <- counter+1
}
}
}
rconfig.mtx
}
