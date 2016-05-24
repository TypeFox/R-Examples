## ox - interactive example

require(graphics)

fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
    x1 <- x[1]
    x2 <- x[2]
    c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
       200 *      (x2 - x1 * x1))
}
ans1<-optimx(c(-1.2,1), fr)
ans1
tmpin<-readline("cont?")
str(ans1)
tmpin<-readline("cont?")
print(attr(ans1,"details"))
tmpin<-readline("cont?")
coef(ans1)
print(dim(coef(ans1)))
tmpin<-readline("cont?")
summary(ans1)
tmpin<-readline("cont?")
best1<-summary(ans1, sort.order=value)[1, ]
best1
dim(best1)

tmpin<-readline("cont?")

####
##JN Since this is a single method, the details could have wrong structure
ans2<-optimx(c(-1.2,1), fr, grr, method = "BFGS")
ans2
tmpin<-readline("cont?")
str(ans2)
tmpin<-readline("cont?")
det2<-attr(ans2,"details")
##JN cat("dim(det2)\n")
print(dim(det2))
tmpin<-readline("cont?")
coef(ans2)
print(dim(coef(ans2)))
tmpin<-readline("cont?")
ans2d<-as.data.frame(ans2)
print(ans2d)
attributes(ans2d)
tmpin<-readline("cont?")
ans2parval<-ans2[,1:(attr(ans2,"npar")+1)]
print(ans2parval)
class(ans2parval)

tmpin<-readline("cont?")
summary(ans2)
dim(summary(ans2))
tmpin<-readline("cont?")



## The next line will fail if executed because 'hessian = TRUE' no longer allowed
# ans3<-optimx(c(-1.2,1), fr, NULL, method = "BFGS", hessian = TRUE)
ans3<-optimx(c(-1.2,1), fr, NULL, method = "BFGS", hessian = TRUE)
ans3
tmpin<-readline("cont?")
str(ans3)
tmpin<-readline("cont?")
####q

ans4<-optimx(c(-1.2,1), fr, grr, method = "CG",control=list(trace=1))
ans4
tmpin<-readline("cont?")
str(ans4)
tmpin<-readline("cont?")
####

ans5<-optimx(c(-1.2,1), fr, grr, method = "CG", control=list(type=2))
ans5
tmpin<-readline("cont?")
str(ans5)
tmpin<-readline("cont?")
####

ans6<-optimx(c(-1.2,1), fr, grr, method = "L-BFGS-B")
ans6
tmpin<-readline("cont?")
str(ans6)
tmpin<-readline("cont?")
####

flb <- function(x)
    { p <- length(x); sum(c(1, rep(4, p-1)) * (x - c(1, x[-p])^2)^2) }
## 25-dimensional box constrained
ans7<- optimx(rep(3, 25), flb, NULL, method = "L-BFGS-B",
      lower=rep(2, 25), upper=rep(4, 25)) # par[24] is *not* at boundary
ans7
tmpin<-readline("cont?")
str(ans7)
tmpin<-readline("cont?")

## "wild" function , global minimum at about -15.81515
## fw <- function (x)
##    10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80
## plot(fw, -50, 50, n=1000, main = "optim() minimising 'wild function'")

## Suppressed for optimx() ans7 <- optimx(50, fw, method="SANN",
##             control=list(maxit=20000, temp=20, parscale=20))
## ans7a
## Now improve locally {typically only by a small bit}:
## newpar<-unlist(ans7$par) # NOTE: you need to unlist the parameters as optimx() has multiple outputs
##(r2 <- optimx(newpar, fw, method="BFGS"))
##points(r2$par, r2$value, pch = 8, col = "red", cex = 2)

## Show multiple outputs of optimx using all.methods
# genrose function code
genrose.f<- function(x, gs=NULL){ # objective function
## One generalization of the Rosenbrock banana valley function (n parameters)
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	fval<-1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
        return(fval)
}

genrose.g <- function(x, gs=NULL){
# vectorized gradient for genrose.f
# Ravi Varadhan 2009-04-03
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	gg <- as.vector(rep(0, n))
	tn <- 2:n
	tn1 <- tn - 1
	z1 <- x[tn] - x[tn1]^2
	z2 <- 1 - x[tn]
	gg[tn] <- 2 * (gs * z1 - z2)
	gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1
	return(gg)
}

genrose.h <- function(x, gs=NULL) { ## compute Hessian
   if(is.null(gs)) { gs=100.0 }
	n <- length(x)
	hh<-matrix(rep(0, n*n),n,n)
	for (i in 2:n) {
		z1<-x[i]-x[i-1]*x[i-1]
		z2<-1.0-x[i]
                hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
                hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
                hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
                hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
	}
        return(hh)
}

startx<-4*seq(1:10)/3.
ans8<-optimx(startx,fn=genrose.f,gr=genrose.g, hess=genrose.h, control=list(all.methods=TRUE, save.failures=TRUE, trace=0), gs=10)
ans8
tmpin<-readline("cont?")
str(ans8)
tmpin<-readline("cont?")
ans8[, "gevals"]
tmpin<-readline("cont?")
ans8["spg", ]
tmpin<-readline("cont?")
summary(ans8, par.select = 1:3)
tmpin<-readline("cont?")
summary(ans8, sort.order = value)[1, ] # show best value
tmpin<-readline("cont?")
head(summary(ans8, sort.order = value)) # best few
tmpin<-readline("cont?")
ans8.best<-summary(ans8, sort.order = value)[1,]
print(ans8.best)
attr(ans8.best,"details")
tmpin<-readline("cont?")

## order by value.  Within those values the same to 3 decimals order by fevals.
summary(ans8, sort.order = list(round(value, 3), fevals), par.select = FALSE)

summary(ans8, sort.order = rownames, par.select = FALSE) # order by method name
tmpin<-readline("cont?")

summary(ans8, par.select = FALSE) # same
tmpin<-readline("cont?")
## Test missing methods
ans8missmeth<-optimx(startx,fn=genrose.f,gr=genrose.g, method=c("spg", "Rcgmin"),
     control=list(save.failures=TRUE, trace=0), gs=10)
# There are just 2 methods here
ans8missmeth
tmpin<-readline("cont?")
# ans8missmeth["Rvmmin", ]
# When we try for a method that is not there in the structure
# the rowname (that is, the method) comes back "NA"
ans8missmeth["Rvmmin", ]
cat("Is method Rvmmin missing? ", (row.names(ans8missmeth["Rvmmin",])[[1]] == "NA"),"\n")

tmp<-readline("Getting items in attributes")
abest<-summary(ans8, order=value)[1,]
print(abest)
attributes(abest)
abestdetails<-attr(abest, "details")
abestdetails[2]
abestdetails[3]
abestdetails[,"ngatend"]
abestdetails[,"hev"]



tmpin<-readline("Polyalgorithm")

startx<-4*seq(1:10)/3.

####

## "Polyalgorithm with 200 steps NM followed by up to 75 of ucminf
ans9<-optimx(startx,fn=genrose.f,gr=genrose.g, hess=genrose.h, method=c("Nelder-Mead","ucminf"),
             itnmax=c(200,75), control=list(follow.on=TRUE, save.failures=TRUE,trace=0), gs=10)
ans9
tmpin<-readline("cont? Again, but with trace on")
ans9<-optimx(startx,fn=genrose.f,gr=genrose.g, hess=genrose.h, method=c("Nelder-Mead","ucminf"),
             itnmax=c(200,75), control=list(follow.on=TRUE, save.failures=TRUE,trace=1), gs=10)

####

startx<-4*seq(1:10)/3.
## 200 steps NM is not enough to terminate
ans10<-optimx(startx,fn=genrose.f,gr=genrose.g, method=c("Nelder-Mead"),
             itnmax=c(200), control=list(trace=0, save.failures=FALSE), gs=10)
## The answer should be NULL
ans10
tmpin<-readline("cont?")
str(ans10)
tmpin<-readline("cont?")

####

startx<-4*seq(1:10)/3.
## Try getting hessian but not kkt tests
ans11<-optimx(startx,fn=genrose.f,gr=genrose.g, hessian=TRUE,  
	control=list(all.methods=TRUE, trace=0, save.failures=FALSE, kkt=FALSE), gs=10)
ans11
tmpin<-readline("cont?")
str(ans11)
tmpin<-readline("cont?")

## Method 7 should be Rvmmin
attr(ans11,"details")[7]

####

startx<-4*seq(1:10)/3.
## Use analytic hessian and no KKT tests
ans12<-optimx(startx,fn=genrose.f,gr=genrose.g, hess=genrose.h, 
      control=list(trace=0, save.failures=FALSE, kkt=FALSE), gs=10)
ans12
attr(ans12,"details")
tmpin<-readline("cont?")
str(ans12)
tmpin<-readline("cont?")

#### 

## Maximization test

maxfn<-function(x) {
      	n<-length(x)
	ss<-seq(1,n)
	f<-10-(crossprod(x-ss))^2
	f<-as.numeric(f)
	return(f)
}

x0<-rep(pi,4)
ans.mx<-optimx(x0,maxfn,control=list(maximize=TRUE,all.methods=TRUE,save.failures=TRUE,trace=TRUE))
ans.mx
tmpin<-readline("cont?")
str(ans.mx)
tmpin<-readline("cont?")
best.mx<-summary(ans.mx, sort.order = value)[1, ]
best.mx
## Check that hessian eigenvalues are negative
attr(best.mx,"details")["hev"]
tmpin<-readline("cont?")
attr(best.mx,"details")$nhatend
tmpin<-readline("cont?")

cat("test names on parameters\n")

## Bates model: y = Asym/(1+exp((xmid-t)/scal))
##      = b1 / (1 + exp(xmid/scal)*exp(-t/scal))
##      = b1 / (1 + b2 * exp(-b3*t))
##  where b2 = exp(xmid/scal)
##        b3 = 1/scal        b1 = Asym
##    So Asym = b1,  scal = 1/b3, xmid = log(b2) * scal = log(b2)/b3



hbates.res<-function(x) { 
    if(length(x) != 3) stop("hobbs.res -- parameter vector n!=3")
    y<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
         38.558, 50.156, 62.948, 75.995, 91.972)
    t<-1:12
#    if(abs(12*x[3])>50) {
#       res<-rep(Inf,12)
#    } else {
       res<-x[1]/(1+exp((x[2]-t)/x[3])) - y
#    }
}

hbates.f<- function(x){ # # Hobbs weeds problem -- function
#    if (abs(12*x[3]) > 500) { # check computability
#       fbad<-.Machine$double.xmax
#       return(fbad)
#    }
    res<-hbates.res(x)
    f<-sum(res*res)
}

hbates.jac<-function(x){ # Jacobian of Hobbs weeds problem
   jj<-matrix(0.0, 12, 3)
   t<-1:12
    yy<-exp((x[2]-t)/x[3])
    zz<-1.0/(1+yy)
     jj[ ,1] <- zz
     jj[ ,2] <- -x[1]*zz*zz*yy/x[3]
     jj[ ,3] <- x[1]*zz*zz*yy*(x[2]-t)/(x[3]*x[3])
   return(jj)
}

hbates.g<-function(x){ # gradient of Hobbs weeds problem
    # NOT EFFICIENT TO CALL AGAIN
    jj<-hbates.jac(x)
    res<-hbates.res(x)
    gg<-as.vector(2.*t(jj) %*% res)
    return(gg)
}

start<-c(200,10,10)
names(start)<-NULL # just to be safe
# no names
anoname<-optimx(start, hbates.f, hbates.g, control=list(all.methods=TRUE))
print(summary(anoname))

names(start)<-c("Asym", "xmid", "scal")

anames<-optimx(start, hbates.f, hbates.g, control=list(all.methods=TRUE))
print(summary(anames))

tmpin<-readline("cont?")

cat("DONE!")



