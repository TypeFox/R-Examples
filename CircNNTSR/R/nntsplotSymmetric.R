nntsplotSymmetric <-
function(cpars=c(0,0),M=0,...){
if (M==0)
{
x<-rep(1/(2*pi),2)
return(plot(c(0,2*pi),x,type="l",xlab="theta"))
}

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

cparsSymmetric<-c(cpars[1:M],array(0,M))

nntsplotint<-function(theta){
res <- nntsABDensitySymmetric(cparsSymmetric,M,theta-cpars[(M+1)]) 
}
return(curve(nntsplotint,0,2*pi,xlab="theta",...))
}

