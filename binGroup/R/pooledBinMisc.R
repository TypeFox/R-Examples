bracket.root <-
function(f, lower=0, upper=1, ..., scaling = 1.6, direction = c("both", "up", 
"down"), max.iter = 50)
{
direction <- match.arg(direction)
if(lower == upper)
stop("lower cannot equal upper")
f1 <- f(lower, ...)
f2 <- f(upper, ...)
done <- FALSE
num.iter <- 0
while(!done) {
#cprint(lower)
#cprint(upper)
num.iter <- num.iter + 1
if(num.iter >= max.iter)
stop(paste("root not bracketed in", max.iter, 
"iterations\n"))
if(f1 * f2 < 0) {
done <- TRUE
ans <- c(lower, upper)
}
switch(direction,
both = if(abs(f1) < abs(f2)) f1 <- f(lower <- lower +
scaling * (lower - upper), ...) else f2 <-
f(upper <- upper + scaling * (upper -
lower), ...),
down = f1 <- f(lower <- lower + scaling * (lower - 
upper), ...),
up = f2 <- f(upper <- upper + scaling * (upper - lower),
...))
}
sort(ans)
}



bracket.bounded.root <-
function(f,lower=0,upper=1,lbnd=0,ubnd=1,...,slicing=2,max.iter=50)
{
if(lower == upper) stop("lower cannot equal upper")
f1 <- f(lower,...)
f2 <- f(upper,...)
done <- FALSE
num.iter <- 0
while(!done){
num.iter <- num.iter + 1
if(num.iter >= max.iter) stop(paste("root not bracketed in",max.iter,"iterations\n"))
if(f1*f2 < 0) {
done <- TRUE
ans <- c(lower,upper)
} 
# march toward the bound in slicing steps
if(abs(f1) < abs(f2))
f1 <- f(lower <- (lower+lbnd)/slicing,...) 
else 
f2 <- f(upper <- (upper+ubnd)/slicing,...)
}
sort(ans)
}


printc <-
function(x,digits=options()$digits)
{
namex <- deparse(substitute(x))
cat(paste(namex,"=",round(x,digits=digits),"\n"))
}

score2 <-
function(s,d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
{
0.5*(score.p((s+d)/2,x1,m1,n1) + score.p((s-d)/2,x2,m2,n2))
}


score2.vec <-
function(s,d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
{
ns <- length(s)
s2 <- vector(length=ns)
for(i in 1:ns) s2[i] <- score2(s[i],d,x1,m1,x2,m2,n1,n2)
s2
}



score.p <-
function(p, x, m, n = rep(1,length(m)))
{
if(sum(x) == 0 & p >= 0) return(-sum(m*n)/(1-p))
if(p <= 0 | p >= 1) return(0)
sum(m*x/(1-(1-p)^m) - m*n)/(1-p)
}


Shatd <-
function(d, x1, m1, x2, m2, n1, n2)
{
#if(sum(x1)==0 & sum(x2)==0){
#N1 <- sum(m1*n1)
#N2 <- sum(m2*n2)
#return(2 - abs(N1-N2)*abs(d)/(N1+N2))
#}
if(sum(x1) == 0){
if(d < 0){
# g is value of dScore(s,d)/ds on the edge s = -d
# must be vectorized for uniroot()
g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
{
nd <- length(d)
gv <- vector(length=nd)
for(i in 1:nd) gv[i] <-  score2(-d[i],d[i],x1,m1,x2,m2,n1,n2)
gv
}
brkt <- bracket.bounded.root(g,lower=-0.5,lbnd=-1,upper=-0.001,ubnd=0,
x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
dstar <- uniroot(g,lower = brkt[1],upper = brkt[2], tol = .Machine$double.eps^0.75,
x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
#dstar <- uniroot(g,lower = -1 + .Machine$double.eps^0.5,upper = 0 -.Machine$double.eps^0.5,
#tol = .Machine$double.eps^0.75,
#x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
if(d <= dstar){ 
return(-d)
} else {

brkt <- bracket.bounded.root(score2,lower= -0.5,lbnd=-d,upper=0,ubnd=2+d,
d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
ans <- uniroot(score2,lower = brkt[1],upper = brkt[2],
tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root

#ans <- uniroot(score2,lower = -d + .Machine$double.eps^0.5,upper = 2 + d - .Machine$double.eps^0.5,
#tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
return(ans)
}
} 
else {

brkt <- bracket.bounded.root(score2,lower= d+0.001,lbnd=d,upper=2-d-0.001,ubnd=2-d,
d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
ans <- uniroot(score2,lower= brkt[1],upper = brkt[2],
tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root

#ans <- uniroot(score2,lower= d + .Machine$double.eps^0.5,upper = 2 - d -.Machine$double.eps^0.5,
#tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
    return(ans)
}
}
if(sum(x2) == 0){
if(d > 0){
# g is value of dScore(s,d)/ds on the edge s = d
# must be vectorized for uniroot()
g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
{
nd <- length(d)
gv <- vector(length=nd)
for(i in 1:nd) gv[i] <- score2(d[i],d[i],x1,m1,x2,m2,n1,n2)
gv
}
brkt <- bracket.bounded.root(g,lower=0.001,lbnd=0,upper=0.5,ubnd=1,
x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
dstar <- uniroot(g,lower = brkt[1],upper = brkt[2],
tol = .Machine$double.eps^0.75,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
#dstar <- uniroot(g,lower = 0 + .Machine$double.eps^0.5,upper = 1 -.Machine$double.eps^0.5,
#tol = .Machine$double.eps^0.75,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
if(d >= dstar) 
return(d)
else {

brkt <- bracket.bounded.root(score2,lower=d+0.001,lbnd=d,upper=2-d-0.001,ubnd=2-d,
d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
ans <- uniroot(score2,lower=brkt[1],upper = brkt[2],
tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root

#ans <- uniroot(score2,lower= d + .Machine$double.eps,upper = 2 - d - .Machine$double.eps,
#tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
return(ans)
}
} 
else {
brkt <- bracket.bounded.root(score2,lower=-d+0.001,lbnd=-d,upper=2+d-0.001,ubnd=2+d,
d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
ans <- uniroot(score2,lower= brkt[1],upper=brkt[2],
tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
#ans <- uniroot(score2,lower= -d+.Machine$double.eps^0.5,upper=2+d-.Machine$double.eps^0.5,
#tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
    return(ans)
}
}
if(d < 0){
brkt <- bracket.bounded.root(score2,lower=-d+0.001,lbnd=-d,upper=2+d-0.001,ubnd=2+d,
d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
ans <- uniroot(score2, lower =  brkt[1],
upper = brkt[2], 
tol = .Machine$double.eps^0.75, d = d, x1 = x1,
m1 = m1, x2 = x2, m2 = m2, n1 = n1,
n2 = n2)$root

#ans <- uniroot(score2, lower =  - d + .Machine$double.eps^0.5,
#upper = 2 + d - .Machine$double.eps^0.5, 
#tol = .Machine$double.eps^0.75, d = d, x1 = x1,
#m1 = m1, x2 = x2, m2 = m2, n1 = n1,
#n2 = n2)$root
}
else {
brkt <- bracket.bounded.root(score2,lower=d+0.025,lbnd=d,upper=2-d-0.05,ubnd=2-d,
d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
ans <- uniroot(score2, lower = brkt[1], upper = brkt[2],
tol = .Machine$double.eps^0.75,d = d, x1 = x1,m1 = m1, x2 = x2, m2 = m2, n1 = n1,n2 = n2)$root
#ans <- uniroot(score2, lower = d + .Machine$double.eps^0.5, upper = 2 - d - .Machine$double.eps^0.5,
#tol = .Machine$double.eps^0.75,d = d, x1 = x1,m1 = m1, x2 = x2, m2 = m2, n1 = n1,n2 = n2)$root
}
ans
}



Shatd.vec <-
function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
{
nd <- length(d)
shat <- vector(length=nd)
dhat <- pooledbinom.diff.mle(x1,m1,x2,m2,n1,n2)
for(i in 1:nd){
 shat[i] <- Shatd(d[i],x1,m1,x2,m2,n1,n2)
 }
shat
}

Z <-
function(d, s, x1, m1, x2, m2, n1=rep(1,length(x1)), n2=rep(1,length(x2)))
 {
0.5*(score.p((d + s)/2, x1, m1, n1) - score.p((s -
d)/2, x2, m2, n2)) * sqrt((pooledbinom.mle.var(
(d + s)/2, m1, n1) + pooledbinom.mle.var((
s - d)/2, m2, n2)))
 }



checkpooledBin<-function(x,m,n,alpha)
{

if (length(x)!=length(n)||length(x)!=length(m))
 {stop("Vectors x, n, m must be vectors of the same length")}

if ( any(x>n) )
 { xgn <- which(x>n)
 stop(paste("Element(s)",xgn, "of the vector x (number of positive pools) is (are) greater than the corresponding element(s) in vector n (number of pools)"))
 }

if(any(abs(round(n)-n) > 1e-07) | any(n<1))
 {stop("Number of pools n must be a vector of integer values > 0")}

if(any(abs(round(m)-m) > 1e-07) | any(m<1))
 {stop("Pool sizes m must be a vector of integer values > 0")}

if(any(abs(round(x)-x) > 1e-07) | any(x<0))
 {stop("Number of positive pools x must be a vector of non-negative integer values > 0")}

if(alpha<=0 | alpha>=1)
 {stop("Argument alpha must be a numeric value between 0 and 1")}

}


checkpooledBinDiff<-function(x1, x2, m1, m2, n1, n2, alpha)
{

if (length(x1)!=length(n1)||length(x1)!=length(m1)) {stop("Vectors x1, n1, m1 must be vectors of the same length")}
if (length(x2)!=length(n2)||length(x2)!=length(m2)) {stop("Vectors x2, n2, m2 must be vectors of the same length")}

if ( any(x1>n1) )
 { xgn1 <- which(x1>n1); stop(paste("Element(s)",xgn1, "of the vector x1 (number of positive pools in population 1) is (are) greater than the corresponding element(s) in vector n1 (number of pools in population 1)")) }
if ( any(x2>n2) )
 { xgn2 <- which(x2>n2); stop(paste("Element(s)",xgn2, "of the vector x2 (number of positive pools in population 2) is (are) greater than the corresponding element(s) in vector n2 (number of pools in population 2)")) }

if(any(abs(round(n1)-n1) > 1e-07) | any(n1<1)) {stop("Number of pools n1 must be a vector of integer values > 0")}
if(any(abs(round(n2)-n2) > 1e-07) | any(n2<1)) {stop("Number of pools n2 must be a vector of integer values > 0")}


if(any(abs(round(m1)-m1) > 1e-07) | any(m1<1)) {stop("Pool sizes m1 must be a vector of integer values > 0")}
if(any(abs(round(m2)-m2) > 1e-07) | any(m2<1)) {stop("Pool sizes m2 must be a vector of integer values > 0")}

if(any(abs(round(x1)-x1) > 1e-07) | any(x1<0)) {stop("Number of positive pools x must be a vector of non-negative integer values > 0")}
if(any(abs(round(x2)-x2) > 1e-07) | any(x2<0)) {stop("Number of positive pools x must be a vector of non-negative integer values > 0")}

if(alpha<=0 | alpha>=1) {stop("Argument alpha must be a numeric value between 0 and 1")}

}

