nntsloglikSymmetric <-
function(cpars=c(0,0),M=0,data){

size<-length(cpars)-1
if (size<M)
{
temp<-size+1
cparscorr<-c(cpars[1:size],array(0,(M-temp+1)),cpars[temp])
cpars<-cparscorr
cat("Warning: Missing parameters set to 0
")
}

y <- 0
for (k in 1:length(data))
{
y <- y - log(nntsSymmetricDensity(cpars,M,data[k]))
}
res <- Re(y)
return(res)
}

