library(Sim.DiffProc)


### 2-dim Ito bridge sde
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

