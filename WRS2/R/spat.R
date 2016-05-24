spat<-function(x){
#
# compute spatial median
# x is an n by p matrix
#
if(!is.matrix(x))stop("x must be a matrix")
x<-elimna(x)
START<-apply(x,2,median)
#val<-nelder(x,ncol(x),spat.sub,START=START)
val=optim(START,spat.sub,x=x,method='BFGS')$par             
val
}
