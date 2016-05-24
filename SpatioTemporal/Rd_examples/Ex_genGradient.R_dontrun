#create a two variable function
f.test <- function(x){sin(x[1])*cos(x[2])}

#compute the gradient using forward difference
genGradient(c(.5,.5), f.test, diff.type=1)
#and central difference
genGradient(c(.5,.5), f.test, diff.type=0)
#compared to the true value
c(cos(.5)*cos(.5),-sin(.5)*sin(.5))

#Compute the Hessian
genHessian(c(.5,.5), f.test, h=1e-4)
#and compare to the true value
matrix(c(-sin(.5)*cos(.5),-cos(.5)*sin(.5),
         -cos(.5)*sin(.5),-sin(.5)*cos(.5)),2,2)

\dontshow{
  if( max(abs(genGradient(c(.5,.5), f.test, h=1e-5, diff.type=0) -
              c(cos(.5)*cos(.5),-sin(.5)*sin(.5)))) > 1e-10 ){
    stop("genGradient, diff.type=0: Results not equal")
  }
  if( max(abs(genGradient(c(.5,.5), f.test, h=1e-5, diff.type=1) -
              c(cos(.5)*cos(.5),-sin(.5)*sin(.5)))) > 1e-5 ){
    stop("genGradient, diff.type=1: Results not equal")
  }
  if( max(abs(genGradient(c(.5,.5), f.test, h=1e-5, diff.type=-1) -
              c(cos(.5)*cos(.5),-sin(.5)*sin(.5)))) > 1e-5 ){
    stop("genGradient, diff.type=-1: Results not equal")
  }
  if( max(abs(genHessian(c(.5,.5), f.test, h=1e-4) -
              matrix(c(-sin(.5)*cos(.5),-cos(.5)*sin(.5),
                       -cos(.5)*sin(.5),-sin(.5)*cos(.5)),2,2))) > 1e-5 ){
    stop("genHessian: Results not equal")
  }
}

