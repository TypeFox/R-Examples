zip2image <- 
function( zip, line ) {
     im <- zip[line, ]
     print(paste("digit ", im[1], " taken"))
     im <- im[-1]
     im <- t(matrix(im, 16, 16, byrow=TRUE))
     im <- im[, 16:1]
     return(im) }

# To be plotted by:
# image( im, col=gray(256:0/256), zlim=c(0,1) )
#
