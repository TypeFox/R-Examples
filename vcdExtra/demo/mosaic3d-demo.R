### mosaic3d-demo: proof-of-concept for exploring 3D mosaic plots


##
## split a 3D object along a given dimension, dim, into copies whose
## extent along that dimension are given by the proportions in vector p
## (rescaled to proportions if they are not already so).
##
## The objects are slightly separated along that dimension, allowing
## a total inter-object space = space
##

split3d <- function(obj, p, dim, space=.10) {
  range <-range3d(obj)
  min <- range[1,]
  p <- p/sum(p)                 # assure proportions
  uspace <- space/(length(p)-1) # unit space between objects
  scales <- p * (1-space)
  shifts <- c(0, cumsum(p)[-length(p)])*diff(range[,dim])
  result <- list()
  for (i in seq_along(p)) {
     xscale <- yscale <- zscale <- 1
     xshift <- yshift <- zshift <- 0

     if (dim == 1) {
       xscale <- scales[i]
       xshift <- shifts[i] + min[1]*(1-xscale) + (uspace * (i-1))
     } else if (dim == 2) {
       yscale <- scales[i]
       yshift <- shifts[i] + min[2]*(1-yscale) + (uspace * (i-1))
     } else if (dim == 3) {
       zscale <- scales[i]
       zshift <- shifts[i] + min[3]*(1-zscale) + (uspace * (i-1))
     }

     result[[i]] <- translate3d(scale3d(obj, xscale, yscale, zscale),
                                xshift, yshift, zshift)

   }
   result
}


range3d <- function(obj) {
	if (!"vb" %in% names(obj)) stop("Not a mesh3d or shape3d object")
  x <- with(obj, range(vb[1,]/vb[4,]))
  y <- with(obj, range(vb[2,]/vb[4,]))
  z <- with(obj, range(vb[3,]/vb[4,]))
  result <- cbind(x,y,z)
  rownames(result)<- c('min', 'max')
  result
}

label3d <- function(objlist, dim, text, offset=.1, adj=c(0.5, 1), ...) {
	ranges <- lapply(objlist, range3d)
	loc <- t(sapply(ranges, colMeans))   # positions of labels on dimension dim
	min <- t(sapply(ranges, function(x) x[1,]))  # other dimensions at min values
	xyz <- min - offset
	xyz[,dim] <- loc[,dim]
	texts3d(xyz, texts=text, adj=adj, ...)
	
}


library(rgl)
# level 1
open3d()

# use transparent colors for the side walls of mosaic cubes
crgb <- col2rgb(c("red", "gray90", "blue"))/255
clr <-rbind(crgb, alpha=0.5)
col <- rgb(clr[1,], clr[2,], clr[3,], clr[4,])
#col <- c("#FF000080", "#E5E5E580", "#0000FF80")

sl0 <- cube <- cube3d(alpha=0.3)
sl1 <- split3d(sl0, c(.2, .3, .5), 1)
shapelist3d(sl1, col=col)
label3d(sl1, 1, c("A1", "A2", "A3"))

# level 2
open3d()
sl2 <- list()
for (i in seq_along(sl1)) {
	p <- runif(1, .2, .8)
	sl2 <- c(sl2, split3d(sl1[[i]], c(p, 1-p), 2, space=.1))
	}
shapelist3d(sl2, col=col)
label3d(sl1, 1, c("A1", "A2", "A3"))
label3d(sl2[1:2], 2, c("B1", "B2"))

# level 3
open3d()
sl3 <- list()
for (i in seq_along(sl2)) {
	p <- runif(1, .2, .8)
	sl3 <- c(sl3, split3d(sl2[[i]], c(p, 1-p), 3, space=.05))
	}
shapelist3d(sl3, col=col)
label3d(sl1, 1, c("A1", "A2", "A3"))
label3d(sl2[1:2], 2, c("B1", "B2"))
#label3d(sl3[1:2], 3, c("C1", "C2"))
label3d(sl3[1:2], 3, c("C1", "C2"), adj=rev(c(0.5, 1)))

