mesh1 <- SolidSimplex( n=2, k=5 )
print(str(mesh1))
mesh2 <- mvmeshFromSimplices(  mesh1$S )
print(str(mesh2))
mesh3 <-mvmeshFromSVI(  mesh1$V, mesh1$SVI )
print(str(mesh3))
mesh4 <-mvmeshFromVertices( mesh1$V )
print(str(mesh4))

# sample uniformly from the simplex
x <- rtessellation( n=1000, mesh1$S )
plot(mesh1)
points( x )

# drawing and simulating along a path
# trefoil knot
t <- seq(0,2*pi,length=101)
x <- sin(t) + 2*sin(2*t)
y <- cos(t) - 2*cos(2*t)
z <- sin(3*t)
# construct an mvmesh object from the points
V <- cbind( x,y,z)
SVI <- rbind( 1:length(t), 2:(length(t)+1) )
SVI[2,length(t)] <- 1  # join end to start
trefoil <- mvmeshFromSVI( V, SVI )
plot(trefoil, col="blue")
# sample randomly from the path
x <- rtessellation( 100, trefoil$S )
points3d( x, col='red', size=5 )
title3d( "trefoil knot" )

###############################################################
# drawing letters in 3d
###############################################################
make.path3d <- function( V, connect ){
stopifnot( nrow(V)==length(connect)+1 )

k <- 0
SVI <- matrix( 0L, nrow=2, ncol=0)
for (i in 1:length(connect)) {
  if(connect[i]) { SVI <- cbind( SVI, c(i,i+1) ) }
}
return( mvmeshFromSVI( V, SVI ) ) }
###############################################################
weights.by.length <- function( path ){

stopifnot( path$vps==2 )
nS <- dim(path$S)[3]
wt <- rep(0.0,nS)
for (i in 1:nS) {
  wt[i] <- sqrt(sum( (path$S[1,,i]-path$S[2,,i])^2 ) )
}
return(wt)}
###############################################################

V <- matrix( c(1,.2,0, 1,0,0, 0,0,0, 0,1,0, 1,1,0, 1,.8,0),byrow=TRUE, ncol=3 )
connect <- c(1,1,1,1,1)
letter.C <- make.path3d( V, connect )

V <- matrix( c(0,0,0, .5,1,0, 1,0,0, .25,.5,0, .75,.5,0 ),byrow=TRUE, ncol=3 )
connect <- c(1,1,0,1)
letter.A <- make.path3d( V, connect )

V <- matrix( c( 0,0,0, 1,0,0, 1,.5,0, 0,.5,0, 0,1,0, 1,1,0 ),byrow=TRUE, ncol=3 )
connect <- c(1,1,1,1,1)
letter.S <- make.path3d( V, connect )

V <- matrix( c( 0,0,0, 0,1,0, .5,0,0, 1,1,0,  1,0,0 ),byrow=TRUE, ncol=3 )
connect <- c(1,1,1,1)
letter.M <- make.path3d( V, connect )

V <- matrix( c( 0,1,0,  .5,0,0, 1,1,0 ),byrow=TRUE, ncol=3 )
connect <- c(1,1)
letter.V <- make.path3d( V, connect )

V <- matrix( c( 1,0,0, 0,0,0, 0,1,0, 1,1,0, 0,.5,0, .8,.5,0 ),byrow=TRUE, ncol=3 )
connect <- c(1,1,1,0,1)
letter.E <- make.path3d( V, connect )

V <- matrix( c( 0,0,0, 0,1,0, 1,1,0, 1,0,0,  0,.5,0, 1,.5,0 ),byrow=TRUE, ncol=3 )
connect <- c(1,0,1,0,1)
letter.H <- make.path3d( V, connect )

V <- matrix( c( 0,0,0, 0.9,0,0,  1,.1,0, 1,.9,0, .9,1,0, 0,1,0, 0,0,0  ),byrow=TRUE, ncol=3 )
connect <- c(1,1,1,1,1,1)
letter.D <- make.path3d( V, connect )

V <- matrix( c(0,.8,0,  0,1,0,  1,1,0,  1,.5,0, 0,.5,0, 0,0,0, 1,0,0),byrow=TRUE, ncol=3 )
number.2 <- make.path3d( V, c(1,1,1,1,1,1) )

V <- matrix( c(0,0,0,  1,0,0,  1,1,0, 0,1,0, 0,0,0),byrow=TRUE, ncol=3 )
number.0 <- make.path3d( V, c(1,1,1,1) )

V <- matrix( c(.25,0,0, .25,1,0, .1,.8,0),byrow=TRUE, ncol=3 )
number.1 <- make.path3d( V, c(1,1) )

V <- matrix( c(1,1,0,  0,1,0, 0,.5,0, 1,.5,0, 1,0,0,  0,0,0),byrow=TRUE, ncol=3 )
number.5 <- make.path3d( V, c(1,1,1,1,1) )



# combine to get one object
I <- diag(c(1,1,1)) 
word <- mvmeshCombine( letter.M,AffineTransform(letter.V,I,c(1.2,.25,-.5))  )
word <- mvmeshCombine( word, AffineTransform(letter.M,I,c(2.4,0,.25)) )
word <- mvmeshCombine( word, AffineTransform(letter.E,I,c(3.6,.25,-.25)) )
word <- mvmeshCombine( word, AffineTransform(letter.S,I,c(4.8,0,.25)) )
word <- mvmeshCombine( word, AffineTransform(letter.H,I,c(6,.25,-.25)) )

word <- mvmeshCombine( word, AffineTransform(number.2,I,c(1.2,-2,-.5)) )
word <- mvmeshCombine( word, AffineTransform(number.0,I,c(2.4,-2,.25)) )
word <- mvmeshCombine( word, AffineTransform(number.1,I,c(3.6,-2,0)) )
word <- mvmeshCombine( word, AffineTransform(number.5,I,c(4.4,-2,-.5)) )
plot(word )

# sample from the tessellation
x <- rtessellation( n=400, word$S,  weights.by.length( word ) )
points3d( x, col='purple', size=5 )
title3d("sampling from text")

