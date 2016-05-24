############################################################################
#                               Demo 1                                     # 
#                              1-dim SDE                                   #
############################################################################        
set.seed(1234)
f <- expression((2-x)/(1-t))
g <- expression(x)
res1 <- snssde1d(type="str",drift=f,diffusion=g,M=10,x0=1,N=1000)
res1
summary(res1)
moment(res1,order=c(2,3,4))[which(time(res1)==1),]
plot(res1,plot.type="single")
lines(time(res1),mean(res1),col=2,lwd=2)
lines(time(res1),bconfint(res1,level=0.95)[,1],col=4,lwd=2)
lines(time(res1),bconfint(res1,level=0.95)[,2],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of", 95,"% confidence")),inset = .01,
       col=c(2,4),lwd=2,cex=0.8)

############################################################################
#                               Demo 2                                     # 
#                              2-dim SDE                                   #
############################################################################        
set.seed(1234)
fx <- expression( y )
gx <- expression( 0 )
fy <- expression( (4*( 1-x^2 )* y - x) )
gy <- expression( 0.2)

res1 <- snssde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,type="str",T=100,
                 ,N=10000)
res1
plot(res1,pos=2)
dev.new()
plot(res1,union = FALSE)
dev.new()
plot2d(res1,type="n") ## in plane (O,X,Y)
dev.new()
points2d(res1,col=rgb(0,100,0,50,maxColorValue=255), pch=16)

############################################################################
#                               Demo 3                                     # 
#                              3-dim SDE                                   #
############################################################################        
set.seed(1234)
fx <- expression(4*(-1-x)*y)
gx <- expression(0.2)
fy <- expression(4*(1-y)*x)
gy <- expression(0.2)
fz <- expression(4*(1-z)*y)
gz <- expression(0.2)

res <- snssde3d(x0=2,y0=-2,z0=-2,driftx=fx,diffx=gx,drifty=fy,diffy=gy,
                driftz=fz,diffz=gz,N=1000,M=100)
plot(res,pos=2)
dev.new()
plot3D(res,display="persp")

############################################################################
#                               Demo 4                                     # 
#                          2-dim Bridge SDE                                #
############################################################################
set.seed(1234)
fx <- expression(4*(-1-x)*y)
gx <- expression(0.2)
fy <- expression(4*(1-y)*x)
gy <- expression(0.2)

res <- bridgesde2d(x0=c(0,-1),y=c(1,0),driftx=fx,diffx=gx,drifty=fy,diffy=gy,M=50)
res
plot(res)
dev.new()
plot2d(res,type="n")
points2d(res,col=rgb(0,100,0,50,maxColorValue=255), pch=16)

############################################################################
#                               Demo 5                                     # 
#                              1-dim FPT                                   #
############################################################################        
set.seed(1234)
f <- expression( 0.5*x*t )
g <- expression( sqrt(1+x^2) )
St <- expression(-0.5*sqrt(t)+exp(t^2))
res <- fptsde1d(drift=f,diffusion=g,boundary=St,x0=2)
res
plot(res)
summary(res)
dev.new()
plot(density(res$fpt[!is.na(res$fpt)]),main="Kernel Density of a First-Passage-Time")

############################################################################
#                               Demo 6                                     # 
#                              1-dim RN's SDE                              #
############################################################################        
set.seed(1234)
f <- expression( -3*(1+x) )
g <- expression( 0.5*x )
res <- rsde1d(drift=f,diffusion=g,M=100,N=1000,tau=0.5412)
summary(res)
bconfint(res,level=0.95)
moment(res,order=c(2,3,4,5))
plot(res)
dev.new()
plot(density(res$x))


############################################################################
#                               Demo 5                                     # 
#                            Fiting 1-dim SDE                              #
############################################################################  
set.seed(1234)
true <- c(1,-11,2,1,0.5)
pmle <- eval(formals(fitsde.default)$pmle)

fx <- expression(theta[1] + theta[2]*x + theta[3]*x^2)
gx <- expression(theta[4]*x^theta[5])

fres <- lapply(1:4, function(i) fitsde(mydata1,drift=fx,diffusion=gx,
	             pmle=pmle[i],start = list(theta1=1,theta2=1,theta3=1,theta4=1,
				 theta5=1),optim.method = "L-BFGS-B"))
Coef <- data.frame(true,do.call("cbind",lapply(1:4,function(i) coef(fres[[i]]))))
names(Coef) <- c("True",pmle)
Summary <- data.frame(do.call("rbind",lapply(1:4,function(i) logLik(fres[[i]]))),
                      do.call("rbind",lapply(1:4,function(i) AIC(fres[[i]]))),
                      do.call("rbind",lapply(1:4,function(i) BIC(fres[[i]]))),
                      row.names=pmle)
names(Summary) <- c("logLik","AIC","BIC")
Coef	
Summary

############################################################################
#                               Demo 6                                     # 
#                            Bridge 2-dim SDE                              #
############################################################################  
set.seed(1234)
fx <- expression(4*(-1-x)*y)
gx <- expression(0.2)
fy <- expression(4*(1-y)*x)
gy <- expression(0.2)

res <- bridgesde2d(x0=c(0,-1),y=c(1,0),driftx=fx,diffx=gx,drifty=fy,diffy=gy,M=50)
res
plot(res)
dev.new()
plot2d(res,type="n")
points2d(res,col=rgb(0,100,0,50,maxColorValue=255), pch=16)

############################################################################
#                               Demo 7                                     # 
#                            Bridge 3-dim SDE                              #
############################################################################  
set.seed(1234)
fx <- expression(4*(-1-x)*y)
gx <- expression(0.2)
fy <- expression(4*(1-y)*x)
gy <- expression(0.2)
fz <- expression(4*(1-z)*y)
gz <- expression(0.2)

res <- bridgesde3d(x0=c(0,-1,0.5),y=c(0,-2,0.5),driftx=fx,diffx=gx,
                drifty=fy,diffy=gy,driftz=fz,diffz=gz,M=20)
res
plot(res,union=TRUE)
dev.new()
plot3D(res,display = "persp",main="3-dim bridge sde")


