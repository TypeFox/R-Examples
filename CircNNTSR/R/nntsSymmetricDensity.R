nntsSymmetricDensity <-
function(cpars=c(0,0),M=0,theta){

if (M==0) 
return(sqrt(1/(2*pi)))
else 
{
size<-length(cpars)-1
if (size<M)
{
temp<-size+1
cparscorr<-c(cpars[1:size],array(0,(M-temp+1)),cpars[temp])
cpars<-cparscorr
cat("Warning: Missing parameters set to 0
")
}

if (1/(2*pi) - sum(cpars[1:M])<0) 
return("sum of componentes greater than condition") 
else
{
cpars0 <- sqrt(1/(2*pi) - sum(cpars[1:M]))
cparsnew <- c(cpars0,sqrt(cpars[1:M]))
aux <- complex(M+1)
for (k in 0:M)
{
aux[k+1] <- exp(1i*k*(theta - cpars[(M+1)]))
}
aux <- cparsnew*aux
res <- Re(sum(aux)*Conj(sum(aux)))
return(res)
}
}
}

