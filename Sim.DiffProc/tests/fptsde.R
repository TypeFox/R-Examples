library(Sim.DiffProc)

## Example 1: Ito SDE
## dX(t) = -4*X(t) *dt + 0.5*dW(t)
## S(t) = 0 (constant boundary)
set.seed(1234)

f <- expression( -4*x )
g <- expression( 0.5 )
St <- expression(0) 
res1 <- fptsde1d(drift=f,diffusion=g,boundary=St,x0=2)
res1
summary(res1)
plot(res1)
dev.new()
plot(density(res1$fpt[!is.na(res1$fpt)]),main="Kernel Density of a First-Passage-Time")

## Example 2: Ito SDE
## X(t) Brownian motion
## S(t) = 0.3+0.2*t (time-dependent boundary)

f <- expression( 0 )
g <- expression( 1 )
St <- expression(0.5-0.5*t) 
res2 <- fptsde1d(drift=f,diffusion=g,boundary=St)
res2
summary(res2)
plot(res2)
dev.new()
plot(density(res2$fpt[!is.na(res2$fpt)]),main="Kernel Density of a First-Passage-Time")

