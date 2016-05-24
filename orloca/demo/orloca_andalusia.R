# Adjusting graphics parameters
require(grDevices)
require(graphics)
require(png)

# Load data of cities
data(andalusia)

# Creation of new loca.p object
o <- loca.p(x=andalusia$x[1:8], y=andalusia$y[1:8])

# Compute the limits of the map
xmin <- min(andalusia$x)
ymin <- min(andalusia$y)
xmax <- max(andalusia$x)
ymax <- max(andalusia$y)

# Plot of capitals of Andalusia
file = system.file('img', 'andalusian_provinces.png', package='orloca')
img = readPNG(file)
plot(o, img=img, main=gettext('Andalusia'), xleft=xmin, ybottom=ymin, xright=xmax, ytop=ymax)

# Plot of capistals of Andalusia and contour plot
contour(o, img=img, main=gettext('Andalusia'), xleft=xmin, ybottom=ymin, xright=xmax, ytop=ymax)

# Look for optimal location
andalusia.loca.p <- loca.p(andalusia$x[1:8], andalusia$y[1:8])
sol <- zsummin(andalusia.loca.p)
# The optimal solution is located 35 Km northest from Antequera
# Antequera is usually considered as the geographical center of Andalusia
sol
points(sol[1], sol[2], type='p', col='red')

