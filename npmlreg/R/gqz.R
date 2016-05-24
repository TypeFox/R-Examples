#  gqglim.R				Nick Sofroniou (July 13, 2005)

gqz <- function(numnodes=20, minweight=0.000001){
#  Calculate Gaussian Quadrature points for the Normal distribution
#  using the abscissas and weights for Hermite integration. The
#  conversion of the locations and weights is given in Lindsey (1992,
#  page 169:3) and Skrondal & Rabe-Hesketh (2004, page 165:1).
#  The argument numnodes is the theoretical number of quadrature points,
#  locations with weights that are less than the argument minweight will
#  be omitted.
#  The default vale of minweight=0.000001 returns 14 masspoints for the
#  default numnodes=20 as in Aitkin, Francis & Hinde (2005).
    out <- gauss.quad(numnodes, "hermite")  # from statmod
    h <- rbind(out$nodes*sqrt(2), out$weights/sum(out$weights))
#  Sort the locations and weights into columns in decending order of the
#  location vector.
    ord<-order(h[1,], decreasing = TRUE)
    h <- h[,ord]
    h <- cbind(h[1,], h[2,])
    h <- subset(as.data.frame(h), (h[,2] >= minweight))
    names(h) <- c("location","weight")
    h
    }



