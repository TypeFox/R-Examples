library("pixmap")

x <- read.pnm(system.file("pictures/logo.ppm", package="pixmap")[1])
dx <- x@size
chx <- getChannels(x)
stopifnot(is(x, "pixmap"),
          is(x, "pixmapRGB"),
          dx == c(77, 101),
          dim(chx) == c(dx, 3))
par(mfrow = c(2,2))
plot(x, main = "R logo pixmap")
for(j in 1:3)
    plot(pixmapGrey(chx[,, j]), main=paste("channel",j))


x1 <- as(x, "pixmapGrey")
z <- getChannels(x1)
stopifnot(is(z, "matrix"),
          dim(z) == dx)

x2 <- as(x1, "pixmapRGB")
z <- getChannels(x2)
stopifnot(is(z, "array"),
          dim(z) == c(dx, 3))
