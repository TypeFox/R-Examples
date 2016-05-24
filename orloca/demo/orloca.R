# Loading required graphics packages
require(grDevices)
require(graphics)

# Creation of new loca.p object
o <- loca.p(x = c(-1, 1, 0), y = c(0, 0, 1))

# Summaring and printing the object
print(o)

# Plot the demand points
plot(o)

# Evaluation at point (0, .5)
zsum(o, x=0, y=.5)

# Contour plot the objective function
contour(o)

# 3D plot of the objective function
persp(o)

# 3D nice plot
persp(o, col=cm.colors(10000), border=FALSE, shade=TRUE, theta=50, phi=5, ltheta=135)

# Another 3D plot
persp(o, col=cm.colors(10000), border=FALSE, shade=TRUE, theta=50, phi=5, ltheta=135, lphi=90)

# Plots with a background image
if (require('png')) {  
  file = system.file('img', 'spain_provinces.png', package='orloca')
  img = readPNG(file)
  plot(loca.p(x=.55, y=.62), img=img,  xlim=c(0,1), ylim=c(0,1), xleft=0, ybottom=0, xright=1, ytop=1)
  contour(loca.p(x=.55, y=.62), img=img,  xmin=0, ymin=0, xmax=1, ymax=1, xleft=0, ybottom=0, xright=1, ytop=1)
}

# Find the minimum
zsummin(o)

# New random loca.p object with 10 demand point
p <- rloca.p(10)

# Find the minimun 
sol <- zsummin(p)

# Show it
sol

# Eval the function at minimun
zsum(p, sol[1], sol[2])

# Timing the algorithm
system.time(zsummin(p))
