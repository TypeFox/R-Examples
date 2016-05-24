# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function:
# plotR2pre() - a function to draws the prefractal in R^2.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# l - a list with prefractal ($pre) and
#     protofractal ($proto) points & indexes ($index);
# s - a string for the main title.
# Local functions & variables:
# r - a color vector is used to draw the protofractal points; 
# di2str() - a function to describe the protofractal distribution.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
plotR2pre <- function(l=preRIFS(),
                      s="Prefractal points for 3-gon: k=3; p=1/3; mu=1") {
  di2str <- function(di=l$distr) {
    str <- rep("", times=nrow(di))
    for (i in seq(nrow(di)))
      str[i] <- sprintf("k=%.0f; p=%.3g; mu=%.3g", i, 
                        di[i,"p"], di[i,"mu"])    
    return(str)
  }
  if (is.null(l$index)) {
    plot(l$proto, asp=1, main=s)
    points(l$pre, pch=46, col=gray(0.9))
    legend("topright", pch=1, legend=di2str(), bty="n")
  } else {
    r <- rainbow(n=nrow(l$proto), v=0.9)
    plot(l$proto, col=r, asp=1, main=s)
    points(l$pre, pch=46, col=r[l$index])
    legend("topright", pch=1, col=r, legend=di2str(), bty="n")
  }
}