require(BB)
set.seed(123) 

rosbkext.f <- function(x){
n <- length(x)
sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
}

n <- 20
p0 <- rnorm(n)

# 2 constraints: parameters sum to 1 and the last parameter is non-negative
# x[1] + ... + x[n] = 1
# x[n] >= 0
Amat <- rbind(rep(1,n), c(rep(0,n-1),1))
b <- c(1, 0)

#  with projectLinear as in release 2014.1-1 next gave
#  Failure in initial projection!Error in solve.QP(dvec = rep(0, n), Amat = t(A), 
#  bvec = b - c(A %*% par),  : 
#  Error in projection

ans <- spg(par=p0, fn=rosbkext.f, project="projectLinear", projectArgs=list(A=Amat, b=b, meq=1)) 

fuzz <- 5e-7

if(fuzz < max(abs(ans$par -
   c(5.46001058136910467e-01,  3.00133145466337903e-01,  9.30077533280728036e-02,
     1.18680513875347674e-02,  3.38851273609307169e-03,  3.25955271864946869e-03,
     3.25870517822328398e-03,  3.25869904467344929e-03,  3.25870153610197805e-03,
     3.25869504035516781e-03,  3.25869964139831152e-03,  3.25870210072135361e-03,
     3.25869386108257157e-03,  3.25870148633357684e-03,  3.25869879980698884e-03,
     3.25869926113893892e-03,  3.25869852720016805e-03,  3.25856193090499902e-03,
     3.23766981846091914e-03, -6.09863722023096244e-19)))){	
   print(ans$par, digits=18)
   cat("difference:\n")
   print(ans$par -   
   c(5.46001058136910467e-01,  3.00133145466337903e-01,  9.30077533280728036e-02,
     1.18680513875347674e-02,  3.38851273609307169e-03,  3.25955271864946869e-03,
     3.25870517822328398e-03,  3.25869904467344929e-03,  3.25870153610197805e-03,
     3.25869504035516781e-03,  3.25869964139831152e-03,  3.25870210072135361e-03,
     3.25869386108257157e-03,  3.25870148633357684e-03,  3.25869879980698884e-03,
     3.25869926113893892e-03,  3.25869852720016805e-03,  3.25856193090499902e-03,
     3.23766981846091914e-03, -6.09863722023096244e-19),
    digits=18)
   stop("converged to different parameter values!")
   }

if(fuzz < max(abs(ans$value - 17.4152583859326917))){
   print(ans$value, digits=18)
   stop("converged to different function value!")
   }

ans
# $par
#  [1]  5.460011e-01  3.001331e-01  9.300775e-02  1.186805e-02  3.388513e-03
#  [6]  3.259553e-03  3.258705e-03  3.258699e-03  3.258702e-03  3.258695e-03
# [11]  3.258700e-03  3.258702e-03  3.258694e-03  3.258701e-03  3.258699e-03
# [16]  3.258699e-03  3.258699e-03  3.258562e-03  3.237670e-03 -6.098637e-19
# 
# $value
# [1] 17.41526
# 
# $gradient
# [1] 0.0001785687
# 
# $fn.reduction
# [1] 5362.369
# 
# $iter
# [1] 52
# 
# $feval
# [1] 53
# 
# $convergence
# [1] 0
# 
# $message
# [1] "Successful convergence"
