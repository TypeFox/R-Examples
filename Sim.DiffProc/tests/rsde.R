library(Sim.DiffProc)

## rsde 1-dim
## Itô sde
set.seed(1234)

f <- expression( 2*(3-x) )
g <- expression( 1 )
res1 <- rsde1d(drift=f,diffusion=g,M=100,N=1000,tau=0.5412)
res1
summary(res1)
plot(res1,pos=3)
dev.new()
plot(density(res1$x))

## Stratonovich sde 
set.seed(1234)

fx <- expression(4*(-2-x)*t)
gx <- expression(0.4*y)
fy <- expression(0*(y>0)-2*(y<=0))
gy <- expression(x)
res2 <- rsde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,tau=0.6,M=100,
               N=1000,type="str")
res2
summary(res2)
plot(res2,union=FALSE)


## rsde 3-dim
## Itô sde
set.seed(1234)

fx <- expression(2*(3-x))
gx <- expression(y+z)
fy <- expression(2*(3-y))
gy <- expression(x+z)
fz <- expression(2*(3-z))
gz <- expression(x+y)

res3 <- rsde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,
               diffz=gz,N=1000,M=100,Dt=0.01,tau=8)
res3
summary(res3)
plot(res3,pos=3)

