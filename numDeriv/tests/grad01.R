if(!require("numDeriv"))stop("this test requires numDeriv.")

###################################################################
#  sin test. scalar valued function with scalar arg
###################################################################

print(g.anal  <-  cos(pi))

print(g.calcR <-  grad(sin, pi, method="Richardson"))
cat("error: ", err <- max(abs(g.calcR - g.anal)),"\n")
if(err > 1e-11) stop("grad01 test 1 FAILED") # 1e-13 with d=0.01

print(g.calcS <-   grad(sin, pi, method="simple"))
cat("error: ", err <- max(abs(g.calcS - g.anal)),"\n")
if(err > 1e-8) stop("grad01 test 2 FAILED")


###################################################################
#  sin test. vector argument, scalar result
###################################################################

func2a <- function(x) sum(sin(x))

x <- (0:10)*2*pi/10
print(g.anal  <- cos(x))

print(g.calcR <- grad(func2a, x, method="Richardson"))
cat("error: ", err <- max(abs(g.calcR - g.anal)),"\n")
if(err > 1e-10) stop("grad01 test 3 FAILED")

print(g.calcS <- grad(func2a, x, method="simple"))
cat("error: ", err <- max(abs(g.calcS - g.anal)),"\n")
if(err > 1e-4) stop("grad01 test 4 FAILED")


###################################################################
#  sin test. vector argument, vector result
###################################################################

x <- (0:10)*2*pi/10
print(g.anal <-  cos(x))

print(g.calcR <-  grad(sin, x, method="Richardson"))
cat("error: ", err <- max(abs(g.calcR - g.anal)),"\n")
if(err > 1e-10) stop("grad01 test 5 FAILED")# 1e-12 with d=0.01

print(g.calcS <-   grad(sin, x, method="simple"))
cat("error: ", err <- max(abs(g.calcS - g.anal)),"\n")
if(err > 1e-4) stop("grad01 test 6 FAILED")


