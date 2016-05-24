#  check jacobian

if(!require("numDeriv"))stop("this test requires numDeriv.")

x <- pi
print(j.calc <- jacobian(sin, x))
cat("error: ", err <- max(abs(j.calc - cos(x))),"\n")
if( err > 1e-11) stop("jacobian matrix test 1 FAILED") # 1e-13 with d=0.01

x <- (1:2)*2*pi/2
print(j.calc <- jacobian(sin, x))
cat("error: ", err <- max(abs(j.calc -  diag(cos(x)))),"\n")
if( err > 1e-11) stop("jacobian matrix test 2 FAILED") # 1e-13 with d=0.01

func2 <- function(x) c(sin(x), cos(x))

x <- (1:2)*2*pi/2
print(j.calc <- jacobian(func2, x))
cat("error: ", err <- max(abs(j.calc - rbind(diag(cos(x)), diag(-sin(x))))),"\n")
if( err > 1e-11) stop("jacobian matrix test 3 FAILED") # 1e-13 with d=0.01

x <- (0:1)*2*pi
print(j.calc <- jacobian(func2, x))
cat("error: ", err <- max(abs(j.calc - rbind(diag(cos(x)), diag(-sin(x))))),"\n")
if( err > 1e-11) stop("jacobian matrix test 4 FAILED") # 1e-13 with d=0.01


x <- (0:10)*2*pi/10
print(j.calc <- jacobian(func2, x))
cat("error: ", err <- max(abs(j.calc - rbind(diag(cos(x)), diag(-sin(x))))),"\n")
if( err > 1e-10) stop("jacobian matrix test 5 FAILED")# 1e-12 with d=0.01

func3 <- function(x) sum(sin(x)) # R^n -> R

x <- (1:2)*2*pi/2
print(j.calc <- jacobian(func3, x))
cat("error: ", err <- max(abs(j.calc - cos(x))),"\n")
if( err > 1e-11) stop("jacobian matrix test 6 FAILED")# 1e-13 with d=0.01
