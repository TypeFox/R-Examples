
# two dimensional, isotropic
x <- matrix( rnorm(8000), ncol=2 )
histDirectional( x, k=2 )  # default plot.type="radial"
histDirectional( x, k=2, col='red',lwd=4 ) # tinker with color and line width
histDirectional( x, k=2, plot.type="index" )
histDirectional( x, k=2, plot.type="none" )

histRectangular( x, breaks=5, col='blue' ) # default plot.type="pillars"
histRectangular( x, breaks=5, plot.type="counts" )
histRectangular( x, breaks=5, plot.type="index" )
histRectangular( x, breaks=5, plot.type="none" )

histSimplex( x, 4*SolidSimplex( n=2, k=3 )$S, col="green", lwd=3 ) # default plot.type="pillars"
histSimplex( x, 4*SolidSimplex( n=2, k=3 )$S, plot.type="counts" )
histSimplex( x, 4*UnitBall( n=2, k=2 )$S, col="red" )
histSimplex( x, 4*UnitBall( n=2, k=2 )$S, plot.type="counts", col="purple" )

# Specify simplices explicitly to get specific region, e.g. restrict to x[1] >= 0
S1 <- matrix( c(0,0,  10,0,  0,10, 10,10), ncol=2, byrow=TRUE )  # first quadrant (bounded)
S2 <- matrix( c(0,0,  10,0,  0,-10,  10,-10), ncol= 2,, byrow=TRUE ) # fourth quadrant (bounded)
S <- array( c(S1,S2), dim=c(4,2,2) )
simp <- histSimplex( x, S, plot.type="counts" )
text(2,9, paste("nrejects=",simp$nrejects), col='red' )

# check behavior with rejects and ties
r <- histSimplex( x, S, plot.type="counts" )
sum(r$counts,r$nrejects)

x[1:100,1] <- 0.0
r2 <- histSimplex( x, S, plot.type="counts" )
sum(r2$counts,r2$nrejects)

r3 <- histSimplex( x, S, report="all", plot.type="counts" )
sum(r3$counts,r3$nrejects)

par(mfrow=c(2,1),mar=c(2,2,2,1))
plot(r$counts,type='h',ylim=c(0,600)) # should be roughly uniform
plot(r3$counts,type='h',ylim=c(0,600)) # bin 1 count should be elevated



# three dimensional positive data
x <- matrix( abs(rnorm(9000)), ncol=3 )
histDirectional( x, k=3, positive.only=TRUE, col='blue', lwd=3 )
histDirectional( x, k=2, plot.type="orthogonal", positive.only=TRUE, p=1 )


# three dimensional omnidirectional
x <- matrix( rnorm(3000), ncol=3, nrow=1000 )

histDirectional( x, k=2 ) # default plot.type="radial"
histDirectional( x, k=2, plot.type="grayscale" )
histDirectional( x, k=2, plot.type="index" )
histDirectional( x, k=2, plot.type="none" )


# compare frequence, relative freq. and normalized histograms
n <- 20000; d <- 3; k <- 2
x <- matrix( rnorm( n*d ), nrow=n, ncol=d )
dev.new(); par(mfrow=c(3,1),mar=c(4,4,2,1))
histDirectional( x, k, plot.type="index" )
title("omnidirectional data: frequency")
histDirectional( x, k, freq=FALSE, plot.type="index" )
title("omnidirectional data: relative frequency")
histDirectional( x, k, plot.type="index", normalize.by.area=TRUE )
title("omnidirectional data: frequency/surface.area")

# first octant only
dev.new(); par(mfrow=c(3,1),mar=c(4,4,2,1))
histDirectional( abs(x), k, positive.only=TRUE, plot.type="index" )
title("positive data: frequency")
histDirectional( abs(x), k, positive.only=TRUE, freq=FALSE, plot.type="index" )
title("positive data: relative frequency")
histDirectional( abs(x), k, positive.only=TRUE, plot.type="index", normalize.by.area=TRUE )
title("positive data: frequency/surface.area")

histRectangular( x, breaks=3 ) # default plot.type="pillars", stacked 
histRectangular( x, breaks=3, plot.type="counts" ) 
histRectangular( x, breaks=3, plot.type="index"  ) 
histRectangular( x, breaks=3, plot.type="none"  ) 

histSimplex( x, 4*SolidSimplex( n=3, k=3 )$S ) # default plot.type="counts"
histSimplex( x, 4*SolidSimplex( n=3, k=3 )$S, plot.type="counts" )
histSimplex( x, 4*SolidSimplex( n=3, k=3 )$S, plot.type="index"  )
histSimplex( x, 4*SolidSimplex( n=3, k=3 )$S, plot.type="none"  )

histSimplex( x, 4*UnitBall( n=3, k=2 )$S  ) # default plot.type="counts"
histSimplex( x, 4*UnitBall( n=3, k=2 )$S, plot.type="index"  )

# higher dimensional data
d <- 5; n <- 10000
x <- matrix( rnorm(d*n), nrow=n, ncol=d )
histDirectional( x, k=1 )
histDirectional( x, k=1, normalize.by.area=TRUE )
histRectangular( x, breaks=2 )

