## solve G01 problem

## defining the constraint problem: G01
fn<-function(x){
  obj<- sum(5*x[1:4])-(5*sum(x[1:4]*x[1:4]))-(sum(x[5:13]))
  g1<- (2*x[1]+2*x[2]+x[10]+x[11] - 10)
  g2<- (2*x[1]+2*x[3]+x[10]+x[12] - 10)
  g3<- (2*x[2]+2*x[3]+x[11]+x[12] - 10)
  
  g4<- -8*x[1]+x[10]
  g5<- -8*x[2]+x[11]
  g6<- -8*x[3]+x[12]
  
  g7<- -2*x[4]-x[5]+x[10]
  g8<- -2*x[6]-x[7]+x[11]
  g9<- -2*x[8]-x[9]+x[12]
  
  res<-c(obj, g1 ,g2 , g3
         , g4 , g5 , g6
         , g7 , g8 , g9)
  return(res)
}
fName="G01"
d=13
lower=rep(0,d)
upper=c(rep(1,9),rep(100,3),1)
nConstraints=9
set.seed(1)
xStart<-runif(d,min=lower,max=upper)

## Initializing cobra
cobra <- cobraInit(xStart=xStart, fn=fn, fName=fName, lower=lower, upper=upper,
                   nConstraints=nConstraints,  feval=70,
                   initDesPoints=3*d, DOSAC=1,cobraSeed=1)

## Running cobra optimizer
cobra <- startCobra(cobra)
#The true solution is at x* = c(rep(1,9),rep(3,3),1)
#where the true optimum is f(x*) = -15
print(getXbest(cobra))
print(getFbest(cobra))

## Plot the resulting error (best-so-far feasible optimizer result - true optimum)
## on a logarithmic scale:
fb=cobra$df$Best
fb[1:(which(cobra$df$feasible)[1]-1)] <- NA  # invalidate iterates before the 1st feasible point
plot(fb-(-15),log="y",type="l",ylab="error",xlab="iteration")

