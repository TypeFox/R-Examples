require(BB)

# from projectLinear.Rd

fn <- function(x) (x[1] - 3/2)^2 + (x[2] - 1/8)^4

gr <- function(x) c(2 * (x[1] - 3/2) , 4 * (x[2] - 1/8)^3)

Amat <- matrix(c(1, -1, 1, 1, -1, 1, -1, -1), 4, 2, byrow=TRUE)
b <- c(-1, -1, -1, -1)
meq <- 0  # all 4 conditions are inequalities

p0 <-  c( -0.3599199, -1.2219309)

r1 <- spg(par=p0, fn=fn, gr=gr, project="projectLinear", 
        projectArgs=list(A=Amat, b=b, meq=meq))

if(any(1e-14 < r1$par - c(1,0))) stop("projectLinear test 1 failed.")


meq <- 1  # first condition is now an equality
r2 <- spg(par=p0, fn=fn, gr=gr, project="projectLinear", 
        projectArgs=list(A=Amat, b=b, meq=meq))

if(any(1e-14 < r2$par - c(0,1))) stop("projectLinear test 2 failed.")

# box-constraints 

Amat <- matrix(c(1, 0, 0, 1, -1, 0, 0, -1), 4, 2, byrow=TRUE)
b <- c(0, 0, -0.5, -0.5)

r3 <- spg(par=p0, fn=fn, gr=gr, project="projectLinear", 
            projectArgs=list(A=Amat, b=b, meq=0))

if(any(1e-14 < r3$par - c(0.5, 0.1146409327454269))
  ) stop("projectLinear test 3 failed.")

