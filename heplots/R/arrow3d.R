## Original by Barry Rowlingson, R-help, 1/10/2010
## Modified by MF: inlined extprod3d, added barblen (to specify absolute barb length)

arrow3d <- function(p0=c(0,0,0), p1=c(1,1,1), barblen, s=0.05, theta=pi/6, n=3, ...){
 ##      p0: start point
 ##      p1: end point
 ## barblen: length of barb
 ##       s: length of barb as fraction of line length (unless barblen is specified)
 ##   theta: opening angle of barbs
 ##       n: number of barbs
 ##     ...: args passed to lines3d for line styling
 ##
 ## Returns (invisibly): integer ID of the shape added to the scene

 # require(geometry)  - inlined extprod3d
 #require(rgl)
 if (!requireNamespace("rgl")) stop("rgl package is required.")    
 

	# cross product of 3D vectors
  extprod3d <- function (x, y) 
  {
      x = matrix(x, ncol = 3)
      y = matrix(y, ncol = 3)
      drop(cbind(x[, 2] * y[, 3] - x[, 3] * y[, 2], x[, 3] * y[, 
          1] - x[, 1] * y[, 3], x[, 1] * y[, 2] - x[, 2] * y[, 
          1]))
  }


 ## rotational angles of barbs
 phi=seq(0,2*pi,len=n+1)[-1]

 ## length of line
 lp = sqrt(sum((p1-p0)^2))

 if (missing(barblen)) {
 	barblen <- s*lp
 }
 
 ## point down the line where the barb ends line up
 cpt=(1-(barblen*cos(theta)))*(p1-p0)

 ## draw the main line
 line = rgl::lines3d(c(p0[1],p1[1]),c(p0[2],p1[2]),c(p0[3],p1[3]),...)

 ## need to find a right-angle to the line. So create a random point:
 rpt = jitter(c(
   runif(1,min(p0[1],p1[1]),max(p0[1],p1[1])),
   runif(1,min(p0[2],p1[2]),max(p0[2],p1[2])),
   runif(1,min(p0[3],p1[3]),max(p0[3],p1[3]))
   ))

 ## and if it's NOT on the line the cross-product gives us a vector at right angles:
 r = extprod3d(p1-p0,rpt)
 ## normalise it:
 r = r / sqrt(sum(r^2))

 ## now compute the barb end points and draw:
 pts = list()
 for(i in 1:length(phi)){
   ptb=rgl::rotate3d(r,phi[i],(p1-p0)[1],(p1-p0)[2],(p1-p0)[3])
   rgl::lines3d(
           c(p1[1],cpt[1]+p0[1]+barblen*sin(theta)*ptb[1]),
           c(p1[2],cpt[2]+p0[2]+barblen*sin(theta)*ptb[2]),
           c(p1[3],cpt[3]+p0[3]+barblen*sin(theta)*ptb[3]),
           ...
           )
 }
 invisible(line)
}
