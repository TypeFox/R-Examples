# http://www.r-bloggers.com/using-colorized-png-pictograms-in-r-base-plots/

library(png)
img <- readPNG(system.file("img", "Rlogo.png", package="png"))

rimg <- as.raster(img) # raster multilayer object
r <- nrow(rimg) / ncol(rimg) # image ratio
plot(c(0,1), c(0,r), type = "n", xlab = "", ylab = "", asp=1)
rasterImage(rimg, 0, 0, 1, r)

plot(c(0,1), c(0,.6), type = "n", xlab = "", ylab = "", asp=1)
rasterImage(rimg[40:50, 1:6], 0, 0, 1, .6)

# brightness = sqrt(0.299 * R^2 + 0.587 * G^2 + 0.114 * B^2)
#function to calculate brightness values
brightness <- function(hex) {
  v <- col2rgb(hex)
  sqrt(0.299 * v[1]^2 + 0.587 * v[2]^2 + 0.114 * v[3]^2) /255
}

# given a color ramp, map brightness to ramp also taking into account
# the alpha level. The defaul color ramp is grey
#
img_to_colorramp <- function(img, ramp=grey) {
  cv <- as.vector(img)
  b <- sapply(cv, brightness)
  g <- ramp(b)
  a <- substr(cv, 8,9)     # get alpha values
  ga <- paste0(g, a)       # add alpha values to new colors
  img.grey <- matrix(ga, nrow(img), ncol(img), byrow=TRUE)
}

# read png and modify
img <- readPNG(system.file("img", "Rlogo.png", package="png"))
img <- as.raster(img)           # raster multilayer object
r <- nrow(img) / ncol(img)      # image ratio
s <- 3.5                        # size

plot(c(0,10), c(0,3.5), type = "n", xlab = "", ylab = "", asp=1)

rasterImage(img, 0, 0, 0+s/r, 0+s)  # original
img2 <- img_to_colorramp(img)       # modify using grey scale
rasterImage(img2, 5, 0, 5+s/r, 0+s)

plot(c(0,10),c(0,8.5), type = "n", xlab = "", ylab = "", asp=1)

img1 <- img_to_colorramp(img)
rasterImage(img1, 0, 5, 0+s/r, 5+s)

reds <- function(x)
  rgb(colorRamp(c("darkred", "white"))(x), maxColorValue = 255)
img2 <- img_to_colorramp(img, reds)
rasterImage(img2, 5, 5, 5+s/r, 5+s)

greens <- function(x)
  rgb(colorRamp(c("darkgreen", "white"))(x), maxColorValue = 255)
img3 <- img_to_colorramp(img, greens)
rasterImage(img3, 0, 0, 0+s/r, 0+s)

single_color <- function(...) "#0000BB"
img4 <- img_to_colorramp(img, single_color)
rasterImage(img4, 5, 0, 5+s/r, 0+s)
