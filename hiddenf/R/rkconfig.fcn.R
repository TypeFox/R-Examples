rkconfig.fcn <-
function(rows)
{
avalues <- as.vector(names(table(rows)))
a <- length(table(rows))
cc <- 2^(a-1)-1-a
rkconfig.mtx <- matrix(NA,nrow=length(rows),ncol=cc)
counter <- 1
#if(is.even(b))
if(1-(a%%2))
{
if(a==4)
{
for(j in 1:3){
combo <- combn(a,2)
rkconfig.mtx[,j] <- 1*is.element(rows,combo[,j])
}
}
else
{
for(g in 2:((a/2)-1))
{
ulim <- choose(a,g)
for(j in 1:ulim)
{
# print(c(b,g,j,ulim))
combo <- combn(a,g)
rkconfig.mtx[,counter] <- 1*is.element(rows,combo[,j])
counter <- counter+1
}
}
g <- g+1
ulim <- choose(a,g)/2
for(j in 1:ulim)
{
combo <- combn(a,g)[,1:ulim]
rkconfig.mtx[,counter] <- 1*is.element(rows,combo[,j])
counter <- counter+1
}
}
}
#else if (is.odd(b))
else if (a%%2)
{
for(g in 2:((a-1)/2))
{
ulim <- choose(a,g)
for(j in 1:ulim)
{
combo <- combn(a,g)
rkconfig.mtx[,counter] <- 1*is.element(rows,combo[,j])
counter <- counter+1
}
}
}
rkconfig.mtx
}
