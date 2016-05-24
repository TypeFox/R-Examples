# 
# This demo illustrates the use of 'iprior' function that searches 
# for vertices of a convex polyhedron using a set of inequalities 
# or arbitrarily given random points
# 
# updated on 2015.11.23
#

## case 1. a set of linear inequalities in 2 dimensions
ui <- rbind(c(1,0),c(0,1),c(0,-1),c(1,1),c(-2,-1))
ci <- c(1,2,-8,5,-14)
op <- iprior(ui=ui, ci=ci) # (3,8), (1,8), (1,4), (3,2), (6,2)
op$vtx
plot(op)


## case 2. using a set of linear inequalities in 3 dimensions 
ui <- rbind(c(1,0,0),c(0,1,0),c(0,0,1),c(-1,0,0),c(0,-1,0),c(0,0,-1))
ci <- c(1,1,1,-2,-2,-2)
op <- iprior(ui=ui, ci=ci)
op$vtx
plot(op)


## case 3. using arbitrarily points in 2 dimensions
pmat <- matrix(rnorm(20), ncol=2)
op <- iprior(pmat=pmat)
plot(op)


## case 4. using arbitrarily points in 3 dimensions
pmat <- matrix(rnorm(60), ncol=3)
op <- iprior(pmat=pmat)
plot(op)


## case 4. Sphere using arbitrarily points in 3 dimensions
pmat <- matrix(rnorm(300), ncol=3)
pmat <- sqrt(3)*pmat/drop(sqrt((pmat^2) %*% rep(1, 3)))
op <- iprior(pmat=pmat)
plot(op)


## case 5. Convex polyhedron in n-dimensions
pmat <- matrix(rnorm(300), ncol=5)
op <- iprior(pmat=pmat)
op$vtx
