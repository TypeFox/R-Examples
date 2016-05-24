library(maps)
library(RColorBrewer)
library(pbdDEMO, quietly = TRUE)

lon <- TREFHT$def$dim[[1]]$vals               # longitude
lat <- TREFHT$def$dim[[2]]$vals               # latitude
da <- TREFHT$data                             # surface temperature

# Define Axes.
x <- c(lon[lon > 180] -360, lon[lon <= 195])  # adjustment for maps
y <- lat
z <- rbind(da[lon > 180,], da[lon <= 195,])   # adjustment for maps
xlim <- range(x)
ylim <- range(y)
zlim <- range(z)
col.z <- c(colorRampPalette(c("#0000FF", "#2BFCD3"))(100),
           colorRampPalette(c("#2BFCD3", "#5300AB"))(100),
           colorRampPalette(c("#5300AB", "#7CFA82"))(100),
           colorRampPalette(c("#7CFA82", "#A90055"))(100),
           colorRampPalette(c("#A90055", "#D6FC28"))(100),
           colorRampPalette(c("#D6FC28", "#FE0001"))(100))

# Plot
layout(matrix(c(1, 2), ncol = 1), heights = c(3, 1))
par(mar = c(4, 4, 3, 1), mgp = c(2, 1, 0))
plot(NULL, NULL, xlim = xlim, ylim = ylim, type = "n", axes = FALSE,
     xlab = "Longitude", ylab = "Latitude", main = "TREFHT (Jan. 2004)")
image(x, y, z, zlim = zlim, xlim = xlim, ylim = ylim,
      col = col.z, add = TRUE)

# Add Map.
map(add = TRUE)
abline(h = c(-23.5, 0, 23.5), v = 0, lty = 2)
xtickets <- seq(-180, 180, by = 30)
ytickets <- seq(-90, 90, by = 30)
box()
axis(1, at = xtickets, labels = xtickets)
axis(2, at = ytickets, labels = ytickets)

# Add Legend.
z.temp <- matrix(seq(zlim[1], zlim[2], length = 500), ncol = 1)
ztickets <- seq(230, 300, by = 10)
par(mar = c(5, 4, 1, 1), mgp = c(2, 1, 0))
plot(NULL, NULL, xlim = zlim, ylim = c(0, 1), type = "n", axes = FALSE,
     xlab = "TREFHT (Kelvin)", ylab = "")
image(z.temp, 0, z.temp, zlim = zlim, xlim = zlim, ylim = c(0, 1),
      col = col.z, add = TRUE)
axis(1, at = ztickets, labels = ztickets)
