## solve G11 problem nrun times  and plot the results of all nrun runs
nrun=4
feval=30

## Defining the constraint problem: G11
fn <- function(x) {
  y<-x[1]*x[1]+((x[2]-1)^2)
  y<-as.numeric(y)

  g1 <- as.numeric(+(x[2] - x[1]^2))
  
  return(c(objective=y, g1=g1))
}
funcName="G11"
d=2
nConstraints<-1
lower<-c(-1,-1) 
upper<-c(+1,+1) 

## Initializing and running cobra
cobra <- cobraInit(xStart=c(0,0), fn=fn, fName=funcName, lower=lower, upper=upper,
                   nConstraints=nConstraints,  feval=feval,
                   initDesPoints=3*d, DOSAC=1,cobraSeed=1)

mres <- multiCOBRA(fn,lower,upper,nrun=nrun,feval=feval,optim=0.75
                  ,cobra=cobra,funcName=funcName
                  ,ylim=c(1e-12,1e-0),plotPDF=FALSE,startSeed=42)
  

#Two true solutions at x1* = c(-sqrt(0.5),0.5) and x2* = c(+sqrt(0.5),0.5)
#where the true optimum is f(x1*) = f(x2*) = -0.75
print(getXbest(mres$cobra))
print(getFbest(mres$cobra))

print(mres$z2)
