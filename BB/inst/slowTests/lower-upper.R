################################################################
require(BB)

# Note that r0 may not converge with a different seed
set.seed(1234)
p0 <- rnorm(2)

fuzz <- 1e-6

# Extended from example in project.Rd

fn <- function(x) (x[1] - 3/2)^2 + (x[2] - 1/8)^4

gr <- function(x) c(2 * (x[1] - 3/2) , 4 * (x[2] - 1/8)^3)

Amat <- matrix(c(1, 0, 0, 1, -1, 0, 0, -1), 4, 2, byrow=TRUE)
b <- c(0, 0, -0.5, -0.5)

meq <- 0

r0 <- spg(par=p0, fn=fn, gr=gr, project="projectLinear", 
   projectArgs=list(A=Amat, b=b, meq=meq))

if (0 != r0$convergence) stop("lower-upper test 1 did not converge!")

if(fuzz < max(abs(r0$par - c( 0.500000000000000000, 0.136700820055768862)))){
   print(r0$par, digits=18)
   stop("lower-upper test 0 converged to different parameter values!")
   }

if(fuzz < max(abs(r0$value - 1.00000001874412625 ))){
   print(r0$value, digits=18)
   stop("lower-upper test 0 converged to different function value!")
   }

r0

# $par
# [1] 0.5000000 0.1367008
# $value
# [1] 1
# $gradient
# [1] 6.407799e-06
# $fn.reduction
# [1] 6.328745
# $iter
# [1] 10
# $feval
# [1] 11
# $convergence
# [1] 0
# $message
# [1] "Successful convergence"

#############

# Note that the above should be the same as all the following:

#############

r1 <- spg(par=p0, fn=fn, gr=gr, lower=0, upper=0.5, 
    project="projectLinear", 
    projectArgs=list(A=matrix(NA,0,2), b=vector("numeric", 0), meq=0))

if(fuzz < max(abs(r0$par - r1$par))){
   print(r1$par, digits=18)
   stop("lower-upper test 1 converged to different parameter values!")
   }

if(fuzz < max(abs(r0$value - r1$value ))){
   print(r1$value, digits=18)
   stop("lower-upper test 1 converged to different function value!")
   }


#############

r2 <- spg(par=p0, fn=fn, gr=gr, lower=0, upper=0.5)

if(fuzz < max(abs(r0$par - r2$par))){
   print(r2$par, digits=18)
   stop("lower-upper test 2 converged to different parameter values!")
   }

if(fuzz < max(abs(r0$value - r2$value ))){
   print(r2$value, digits=18)
   stop("lower-upper test 2 converged to different function value!")
   }

#############

r3 <- spg(par=p0, fn=fn, gr=gr, lower=c(0,0), upper=c(0.5, 0.5))

if(fuzz < max(abs(r0$par - r3$par))){
   print(r3$par, digits=18)
   stop("lower-upper test 3 converged to different parameter values!")
   }

if(fuzz < max(abs(r0$value - r3$value ))){
   print(r3$value, digits=18)
   stop("lower-upper test 3 converged to different function value!")
   }


#############

r4 <- BBoptim(par=p0, fn=fn, gr=gr, lower= 0, upper= 0.5)

if(fuzz < max(abs(r0$par - r4$par))){
   print(r4$par, digits=18)
   stop("lower-upper test 4 converged to different parameter values!")
   }

if(fuzz < max(abs(r0$value - r4$value ))){
   print(r4$value, digits=18)
   stop("lower-upper test 4 converged to different function value!")
   }

#############
set.seed(12345) # r0 above fails to converge with this seed

pmat <- matrix(rnorm(40), 20, 2)  # 20 starting values each of length 2 

r5  <- multiStart(par=pmat, fn=fn, gr=gr, 
           lower=c(0,0), upper=c(0.5, 0.5), action="optimize")

r5$par[r5$converged, ] #converged solutions

unique(r5$par[r5$converged, ] )

#      [,1]      [,2]
# [1,]  0.5 0.1141247
# [2,]  0.5 0.1146409
# [3,]  0.5 0.1122670
# [4,]  0.5 0.1145620
# [5,]  0.5 0.1135830
# [6,]  0.5 0.1355988
# [7,]  0.5 0.1116386

print(unique(r5$fvalue[r5$converged] ), digits=18)

# [1] 1.00000001398813110 1.00000001151578677 1.00000002628593299
# [4] 1.00000001187043175 1.00000001699065155 1.00000001261909199
# [7] 1.00000003187147524

# unique and choosing [1,] are too sensitive to seed setting, and 
# will not be robust (especially on different platforms)
# if(0.01 < max(abs(r0$par - unique(r5$par[r5$converged,])[1,]))){
#    print(r5$par[r5$converged,], digits=18)
#    stop("lower-upper test 5 converged to different parameter values!")
#    }

if(fuzz < max(abs(r0$value - r5$fvalue[r5$converged]  ))){
   print(r5$fvalue, digits=18)
   stop("lower-upper test 5 converged to different function value!")
   }

################################################################

## Rosenbrock Banana function  from project.Rd with 
##   additional lower and upper constraint

fr <- function(x) { 
x1 <- x[1] 
x2 <- x[2] 
100 * (x2 - x1 * x1)^2 + (1 - x1)^2 
} 
# Impose a constraint that sum(x) = 1 


p0 <- c(0.4, 0.94)    # need feasible starting point

r6 <- spg(par=p0, fn=fr,  lower=c(-0.6, -Inf), upper=c(0.6, Inf),
  project="projectLinear", projectArgs=list(A=matrix(1, 1, 2), b=1, meq=1)) 


if(fuzz < max(abs(r6$par - c( 0.599999994039535522, 0.400000005960464478)))){
   print(r6$par, digits=18)
   stop("lower-upper test 6 converged to different parameter values!")
   }

if(fuzz < max(abs(r6$value - 0.320000109672563315 ))){
   print(r6$value, digits=18)
   stop("lower-upper test 6 converged to different function value!")
   }

################################################################
