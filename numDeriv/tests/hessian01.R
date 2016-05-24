#  check hessian
if(!require("numDeriv"))stop("this test requires numDeriv.")


####################################################################

#   sin tests

####################################################################


x <- 0.25 * pi
print(calc.h <- hessian(sin, x) )
print(anal.h <-  sin(x+pi))
cat("error: ", err <- max(abs(calc.h - anal.h)),"\n")
if( err > 1e-4) stop("hessian test 1 FAILED") # 1e-8 with d=0.01


func1 <- function(x) sum(sin(x))

x <- (0:2)*2*pi/2
#x <- (0:10)*2*pi/10
print(anal.h <-  matrix(0, length(x), length(x)))
print(calc.h <- hessian(func1, x) )
cat("error: ", err <- max(abs(anal.h - calc.h)),"\n")
if( err > 1e-10) stop("hessian test 2 FAILED")

funcD1 <- function(x) grad(sin,x)
print(calc.j <- jacobian(funcD1, x) )
cat("error: ", err <- max(abs(calc.h - calc.j)),"\n")
if( err > 1e-5) stop("hessian test 3 FAILED") # 1e-8 with d=0.01

####################################################################

#   exp tests

####################################################################

fun1e <- function(x) exp(2*x)
funD1e <- function(x) 2*exp(2*x)

x <- 1
print(anal.h <- 4*exp(2*x) )
print(calc.h  <- hessian(fun1e, x) )
cat("\nerror: ", err <- max(abs(calc.h - anal.h)),"\n")
if( err > 1e-3) stop("hessian test 5 FAILED") # 1e-7 with d=0.01

print(calc.j <- jacobian(funD1e, x) )
cat("\nerror: ", err <- max(abs(calc.j - anal.h)),"\n")
if( err > 1e-9) stop("hessian test 6 FAILED") # 1e-10 with d=0.01

fun1e <- function(x) sum(exp(2*x))
funD1e <- function(x) 2*exp(2*x)
x <- c(1,3,5)
print(anal.h <- diag(4*exp(2*x)) )

cat("\n************  d=0.01 works better here.*********\n")
print(calc.h  <- hessian(fun1e, x, method.args=list(d=0.01)) )

cat("\n relative error: \n")
print( err <-  (calc.h - anal.h) /(anal.h+1e-4))
cat("\n max relative error: ", err <- max(abs(err)),"\n")
# above is 901.4512 with d=0.0001

cat("\n error: \n")
print( err <-  calc.h - anal.h)
cat("\n max error: ", err <- max(abs(err)),"\n")
# above is 0.1670381 with d=0.0001

if( err > 1e-5) stop("hessian test 7 FAILED")

print(calc.j <- jacobian(funD1e, x) )
cat("error: ", err <- max(abs(calc.j - anal.h)),"\n")
if( err > 1e-5) stop("hessian test 8 FAILED") # 1e-6 with d=0.01
