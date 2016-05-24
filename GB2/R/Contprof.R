# Contour plot of the profile log-likelihood

contprof.gb2 <- function(z, w=rep(1,length(z)), resol, low=0.1, high=20){
  
# Initial values of a and b under Fisk
x0 <- fisk(z, w)[1:2]    
aa <- seq(low, high, length.out = resol)*x0[1]
bb <- seq(low, high, length.out = resol)*x0[2]
d <- aa %o% bb
for (i in 1:resol){
		for (j in 1:resol){
		d[i,j] <- proflogl.gb2(z, aa[i], bb[j], w)
		                  }
	                }
dlim = range(d, finite = TRUE)
image(aa, bb, d, xlab = "a", ylab = "b", col = heat.colors(16), 
      breaks = seq(dlim[1], dlim[2], length.out=17), main = "Profile log-likelihood")
contour(aa, bb, d, levels = seq(dlim[1], dlim[2], length.out=17), add = TRUE)
# The initial Fisk estimate is added as point "F".
text(x0[1], x0[2], "F")
box()
}