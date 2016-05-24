rbiNormGibbs=function(initx=2,inity=-2,rho,burnin=100,R=500)
{
#
# revision history:
#     P. Rossi 1/05
#
# purpose:
#    illustrate the function of bivariate normal gibbs sampler
#
# arguments:
#   initx,inity  initial values for draw sequence
#   rho  correlation
#   burnin draws to be discarded in final paint
#   R -- number of draws
#
# output:
#   opens graph window and paints all moves and normal contours
#   list containing draw matrix
#
# model:
#  theta is bivariate normal with zero means, unit variances and correlation rho
#
# define needed functions
#
kernel=
function(x,mu,rooti){
# function to evaluate -.5*log of MV NOrmal density kernel with  mean mu, var Sigma
# and with sigma^-1=rooti%*%t(rooti)   
# rooti is in the inverse of upper triangular chol root of sigma
#          note: this is the UL decomp of sigmai not LU!
#                Sigma=root'root   root=inv(rooti)
z=as.vector(t(rooti)%*%(x-mu))
(z%*%z)
}

#
# check input arguments
#
if(missing(rho)) {pandterm("Requires rho argument ")}
#
# print out settings
#
cat("Bivariate Normal Gibbs Sampler",fill=TRUE)
cat("rho= ",rho,fill=TRUE)
cat("initial x,y coordinates= (",initx,",",inity,")",fill=TRUE)
cat("burn-in= ",burnin," R= ",R,fill=TRUE)
cat(" ",fill=TRUE)
cat(" ",fill=TRUE)

sd=(1-rho**2)**(.5)
sigma=matrix(c(1,rho,rho,1),ncol=2)
rooti=backsolve(chol(sigma),diag(2))
mu=c(0,0)

x=seq(-3.5,3.5,length=100)
y=x
z=matrix(double(100*100),ncol=100)
for (i in 1:length(x)) 
{
   for(j in 1:length(y))
   {
   z[i,j]=kernel(c(x[i],y[j]),mu,rooti)
   }
}
prob=c(.1,.3,.5,.7,.9,.99)
lev=qchisq(prob,2)


par(mfrow=c(1,1))
contour(x,y,z,levels=lev,labels=prob,
   xlab="theta1",ylab="theta2",drawlabels=TRUE,col="green",labcex=1.3,lwd=2.0)
title(paste("Gibbs Sampler with Intermediate Moves: Rho =",rho))

points(initx,inity,pch="B",cex=1.5)

oldx=initx
oldy=inity
continue="y"
r=0
draws=matrix(double(R*2),ncol=2)
draws[1,]=c(initx,inity)
cat(" ")
cat("Starting Gibbs Sampler ....",fill=TRUE)
cat("(hit enter or y to display moves one-at-a-time)",fill=TRUE)
cat("('go' to paint all moves without stopping to prompt)",fill=TRUE)
cat(" ",fill=TRUE)
while(continue != "n"&& r < R)
{
  if(continue != "go") continue=readline("cont?")
  newy=sd*rnorm(1) + rho*oldx
  lines(c(oldx,oldx),c(oldy,newy),col="magenta",lwd=1.5)
  newx=sd*rnorm(1)+rho*newy
  lines(c(oldx,newx),c(newy,newy),col="magenta",lwd=1.5)	
  oldy=newy
  oldx=newx
  r=r+1
  draws[r,]=c(newx,newy)
}
continue=readline("Show Comparison to iid Sampler?")
if(continue != "n" & continue != "No" & continue != "no"){
   par(mfrow=c(1,2))
   contour(x,y,z,levels=lev,
      xlab="theta1",ylab="theta2",drawlabels=TRUE,labels=prob,labcex=1.1,col="green",lwd=2.0)
   title(paste("Gibbs Draws: Rho =",rho))
   points(draws[(burnin+1):R,],pch=20,col="magenta",cex=.7)

   idraws=t(chol(sigma))%*%matrix(rnorm(2*(R-burnin)),nrow=2)
   idraws=t(idraws)
   contour(x,y,z,levels=lev,
      xlab="theta1",ylab="theta2",drawlabels=TRUE,labels=prob,labcex=1.1,col="green",lwd=2.0)
   title(paste("IID draws: Rho =",rho))
   points(idraws,pch=20,col="magenta",cex=.7)
}
attributes(draws)$class=c("bayesm.mat","mcmc")
attributes(draws)$mcpar=c(1,R,1)
return(draws)
}
