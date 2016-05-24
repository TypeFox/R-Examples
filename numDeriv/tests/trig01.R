if(!require("numDeriv"))stop("this test requires numDeriv.")

###################################################################
# 3 test functions to test the accuracy of numerical derivatives
#  in "numDeriv" package written by Paul Gilbert
#  Author:  Ravi Varadhan
#  March 27, 2006
###################################################################
options(digits=12)


###################################################################
#         asin test
###################################################################
func1 <- function(x){asin(x)}

x <- c(0.9,0.99,0.999)
exact <- 1/sqrt(1 - x^2)

# With d = 0.0001
print(g.calcS <- grad(func1, x,method.args=list(d=0.0001)))
rel.err <- g.calcS/exact - 1
cbind(x, g.calcS, exact, rel.err)
if(any(rel.err > 1e-10)) stop("trig01 test 1 FAILED")


###################################################################
#         sin test
###################################################################
func2 <- function(x){sin(1/x)}

x <- c(0.1,0.01,0.001,0.0001)
exact <- cos(1/x) * (-1/x^2)

# With d = 0.0001
print(g.calcS <- grad(func2, x,method.args=list(d=0.0001)))

rel.err <- g.calcS/exact - 1
cbind(x, g.calcS, exact, rel.err)
if(any(rel.err > 1e-10)) stop("trig02 test 1 FAILED")



###################################################################
#         power test
###################################################################
func3 <- function(x){(x-100)^2 + 1.e-06 * (x - 300)^3}

x <- c(100.001,300.001)
exact <- 2*(x-100) + 3.e-06*(x-300)^2

# With d = 0.0001
print(g.calcS <- grad(func3, x,method.args=list(d=0.0001)))

rel.err <- g.calcS/exact - 1
cbind(x, g.calcS, exact, rel.err)
if(any(rel.err > 1e-10)) stop("trig03 test 1 FAILED")
