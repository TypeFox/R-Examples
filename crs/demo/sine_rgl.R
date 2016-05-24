## This illustration considers a product sine function and plots the
## results by constructing a 3D real-time rendering plot using OpenGL.

require(crs)
require(rgl)

set.seed(42)

## Interactively request number of observations, whether to do NOMAD
## or exhaustive search, and if NOMAD the number of multistarts

n <- as.numeric(readline(prompt="Input the number of observations desired: "))
cv <- as.numeric(readline(prompt="Input the cv method (0=nomad, 1=exhaustive): "))
cv <- ifelse(cv==0,"nomad","exhaustive")
if(cv=="nomad") nmulti <- as.numeric(readline(prompt="Input the number of multistarts desired (e.g. 10): "))
num.eval <- as.numeric(readline(prompt="Input the number of evaluation observations desired (e.g. 50): "))

x1 <- runif(n,0,1)
x2 <- runif(n,0,1)

dgp <- sin(pi*(x1+x2))^4*sin(pi*x1)^2

y <- dgp + rnorm(n,sd=.1)

model <- crs(y~x1+x2,
             cv=cv,
             complexity="degree-knots",
             knots="uniform",
             deriv=1,
             cv.func="cv.aic",
             nmulti=nmulti)

summary(model)

## Create a 3D rgl perspective plot (need to also assign colors)

x1.seq <- seq(min(x1),max(x1),length=num.eval)
x2.seq <- seq(min(x2),max(x2),length=num.eval)
x.grid <- expand.grid(x1.seq,x2.seq)
newdata <- data.frame(x1=x.grid[,1],x2=x.grid[,2])
z <- matrix(predict(model,newdata=newdata),num.eval,num.eval)

## Number of colors from color palette
num.colors <- 1000
colorlut <- topo.colors(num.colors) 
col <- colorlut[ (num.colors-1)*(z-min(z))/(max(z)-min(z)) + 1 ]

## Open an rgl 3d window and use `persp3d()', a high-level function
## for 3D surfaces (and define the size of the window to be
## 640x640). The function par3d() passes in a window size (the default
## is 256x256 which is quite small), the function rgl.viewpoint()
## allows you to modify the `field of view' to get more of a
## `perspective' feel to the plot, while the function grid3d() adds a
## grid to the plot.

open3d()

par3d(windowRect=c(900,100,900+640,100+640))
rgl.viewpoint(theta = 0, phi = -70, fov = 80)

persp3d(x=x1.seq,y=x2.seq,z=z,
        xlab="X1",ylab="X2",zlab="Y",
        ticktype="detailed",      
        border="red",
        color=col,
        alpha=.7,
        back="lines",
        main="Conditional Mean")

grid3d(c("x", "y+", "z"))

## You can also add other surfaces to the plot (e.g. error bounds) via
## rgl.surface(x, y, z.ub, color="grey", alpha=.7, back="lines")
## rgl.surface(x, y, z.lb, color="grey", alpha=.7, back="lines")
## where z.up and z.lb are the lower and upper bounds

## You could animate the results for 15 seconds using the line
## play3d(spin3d(axis=c(0,0,1), rpm=5), duration=15)
## By default you can manually rotate the figure by dragging the plot
## via your mouse/keypad

## Note - to plot an rgl figure first get it oriented how you want
## (i.e. resize, rotate etc.) and then call rgl.postscript to create,
## i.e., a PDF of your graphic as in
## rgl.postscript("foo.pdf","pdf"). Or better still,
## rgl.snapshot("foo.png") for a png that can be called directly in
## LaTeX via \includegraphics[scale=.6]{foo.png}

## Note also that Sweave support exists as of v0.92.858 and can be
## incorporated per the following illustration:
## <<fig=true, grdevice=rgl.Sweave, pdf=false, stayopen=TRUE>>= 
## x <- rnorm(100); y <- rnorm(100); z <- rnorm(100) 
## plot3d(x, y, z) 
## @ 
