pooledBin <-
function(x,m,n=rep(1,length(x)), 
 pt.method = c("bc-mle","mle","mir"),
 ci.method = c("skew-score","bc-skew-score","score","lrt","wald","mir"),
 scale=1, alpha=0.05, tol=.Machine$double.eps^0.5)
{

checkpooledBin(x=x, m=m, n=n, alpha=alpha)

pt.method <- match.arg(pt.method)
ci.method <- match.arg(ci.method)
if(ci.method=="mir" | pt.method=="mir"){
ci.method <- "mir"
pt.method <- "mir"
}
switch(pt.method,
"mle" = { p <- pooledbinom.mle(x,m,n,tol)},
"bc-mle" ={ p <- pooledbinom.cmle(x,m,n,tol)},
"mir" = {p <- pooledbinom.mir(x,m,n)}
)
if(p < 0 & pt.method=="bc-mle"){
pt.method <- "mle"
warning("Bias-correction results in negative point estimate; using MLE\n")
p <- pooledbinom.mle(x,m,n,tol)
}
switch(ci.method,
"skew-score" = { ci.p <- pooledbinom.cscore.ci(x,m,n,tol,alpha)[2:3]},
"bc-skew-score" ={ ci.p <- pooledbinom.bcscore.ci(x,m,n,tol,alpha)[2:3]},
"score" = {ci.p <- pooledbinom.score.ci(x,m,n,tol,alpha)[2:3]},
"lrt" = {ci.p <- pooledbinom.lrt.ci(x,m,n,tol,alpha)[2:3]},
"wald" = {ci.p <- pooledbinom.wald.ci(x,m,n,tol,alpha)[2:3]},
"mir" = {ci.p <- pooledbinom.mir.ci(x,m,n,tol,alpha)[2:3]}
)
structure(list(p=p,lcl=ci.p[1],ucl=ci.p[2],pt.method=pt.method,ci.method=ci.method,alpha=alpha,x=x,m=m,n=n,scale=scale),class="poolbin")
}

pooledbinom.bcscore.ci <-
function(x, m, n = rep(1., length(x)), tol = 1e-008, alpha = 0.05)
{
f <- function(p, x, mm, n, alpha)
{
gamm <- function(p0, mm, n = rep(1, length(mm)))
{
pooledbinom.mu3(p0, mm, n) * pooledbinom.mle.var(p0,
mm, n)^(3/2)
}
# NOTE: I use the square of the score statistic and chisq
(score.p(p, x, mm, n) * sqrt(pooledbinom.mle.var(p, mm, n)) -
pooledbinom.bias(p, mm, n) - (gamm(p, mm, n) * (qchisq(
1 - alpha, 1) - 1))/6)^2 - qchisq(1. - alpha, 1)
}
# The MLE is just used for a sensible cut-value for the search ...
# it's not needed in the computation
p.hat <- pooledbinom.cmle(x, m, n, tol)
p.hat <- pooledbinom.cmle(x, m, n, tol)
if(sum(x) == 0)
return(pooledbinom.score.ci(x, m, n, tol, alpha))
# Effectively "else"
root.brak <- bracket.bounded.root(f, lower = p.hat/10, lbnd = 0., upper
 = p.hat/2, ubnd = p.hat, x = x, mm = m, n = n, alpha = alpha)
lower.limit <- uniroot(f, lower = root.brak[1], upper = root.brak[
2], tol = .Machine$double.eps^0.75, x = x, mm = m, n = n, alpha
 = alpha)$root
# use mm because m is confused with maxiter
root.brak <- bracket.bounded.root(f, lower = (p.hat + 1)/10, lbnd = 
p.hat, upper = (p.hat + 1)/2, ubnd = 1., x = x, mm = m, n = n,
alpha = alpha)
upper.limit <- uniroot(f, lower = root.brak[1], upper = root.brak[
2], tol = .Machine$double.eps^0.75, x = x, mm = m, n = n, alpha
 = alpha)$root
c(p = p.hat, lower = lower.limit, upper = upper.limit, alpha = alpha
)
}

pooledbinom.bias <-
function(p,m,n=rep(1,length(m)))
{
if(p > 0)
sum((m-1)*m^2*n*(1-p)^(m-3)/(1-(1-p)^m)) * 
pooledbinom.mle.var(p,m,n)^2 / 2
else 0
}

pooledbinom.cmle <-
function(x, m, n = rep(1, length(m)), tol= 1e-8)
{
phat <- pooledbinom.mle(x,m,n)
bias <- pooledbinom.bias(phat,m,n)
phat - bias
}

pooledbinom.cscore.ci <-
function(x, m, n = rep(1., length(x)), tol = 1e-008, alpha = 0.05)
{
# use gamm not gam b/c of the function gam()
f <- function(p, x, mm, n, alpha)
{
gamm <- function(p0, mm, n = rep(1, length(mm)))
{
pooledbinom.mu3(p0, mm, n) * pooledbinom.mle.var(p0,
mm, n)^(3/2)
}
# NOTE: I use the square of the score statistic and chisq
(score.p(p, x, mm, n) * sqrt(pooledbinom.mle.var(p, mm, n)) -
(gamm(p, mm, n) * (qchisq(1 - alpha, 1) - 1))/6)^2 -
qchisq(1. - alpha, 1)
}
# The MLE is just used for a sensible cut-value for the search ...
# it's not needed in the computation
p.hat <- pooledbinom.mle(x, m, n, tol) # used to be cmle
if(sum(x) == 0)
return(pooledbinom.score.ci(x, m, n, tol, alpha))
# Effectively "else"
root.brak <- bracket.bounded.root(f, lower = p.hat/10, lbnd = 0., upper
 = p.hat/2, ubnd = p.hat, x = x, mm = m, n = n, alpha = alpha)
lower.limit <- uniroot(f, lower = root.brak[1], upper = root.brak[
2], tol = .Machine$double.eps^0.75, x = x, mm = m, n = n, alpha
 = alpha)$root
# use mm because m is confused with maxiter
root.brak <- bracket.bounded.root(f, lower = (p.hat + 1)/10, lbnd = 
p.hat, upper = (p.hat + 1)/2, ubnd = 1., x = x, mm = m, n = n,
alpha = alpha)
upper.limit <- uniroot(f, lower = root.brak[1], upper = root.brak[
2], tol = .Machine$double.eps^0.75, x = x, mm = m, n = n, alpha
 = alpha)$root
c(p = p.hat, lower = lower.limit, upper = upper.limit, alpha = alpha
)
}

pooledbinom.I <-
function(p,m,n=rep(1,length(m)))
{
if(p > 0)
sum((m^2 * n * (1 - p)^(m - 2))/(1 - (1 - p)^m))
else 0 # Not sure whether to set this to 0 or not
}

pooledbinom.loglike <-
function(p,x,m,n=rep(1,length(m)))
{
if(p > 0)
sum(x*log((1-(1-p)^m))) + log(1-p)*sum(m*(n-x))
else 0
}

pooledbinom.loglike.vec <-
function(p,x,m,n=rep(1,length(m)))
{
np <- length(p)
lik <- vector(length=np)
for(i in 1:np) lik[i] <- pooledbinom.loglike(p[i],x,m,n)
lik
}

pooledbinom.lrt.ci <-
function(x, m, n = rep(1., length(x)), tol = 1e-008, alpha = 0.05)
{
p.hat <- pooledbinom.mle(x, m, n, tol)
f <- function(p, x, mm, n, f.tol, alpha, p.hat)
{
if(p.hat == 0.)
return( - Inf)
else -2. * (pooledbinom.loglike(p, x, mm, n) - 
pooledbinom.loglike(p.hat, x, mm, n)) - qchisq(
1. - alpha, 1.)
}
if(sum(x) == 0)
return(pooledbinom.score.ci(x, m, n, tol, alpha))
# Effectively "else"
root.brak <- bracket.bounded.root(f, lower = p.hat/10, lbnd = 0., upper
 = p.hat/2, ubnd = p.hat, x = x, mm = m, n = n, f.tol = tol,
alpha = alpha, p.hat = p.hat)
lower.limit <- uniroot(f, lower = root.brak[1], upper = root.brak[
2], tol = .Machine$double.eps^0.75, x = x, mm = m, n = n, f.tol
 = tol, alpha = alpha, p.hat = p.hat)$root
# use mm because m is confused with maxiter
#print(lower.limit)
root.brak <- bracket.bounded.root(f, lower = (p.hat + 1)/10, lbnd = 
p.hat, upper = (p.hat + 1)/2, ubnd = 1., x = x, mm = m, n = n,
f.tol = tol, alpha = alpha, p.hat = p.hat)
upper.limit <- uniroot(f, lower = root.brak[1], upper = root.brak[
2], tol = .Machine$double.eps^0.75, x = x, mm = m, n = n, f.tol
 = tol, alpha = alpha, p.hat = p.hat)$root
c(p = p.hat, lower = lower.limit, upper = upper.limit, alpha = alpha)
}


pooledbinom.mir.ci <-
function(x, m, n = rep(1, length(x)), tol = 1e-008, alpha = 0.05)
{
N <- sum(m*n)
mir <- sum(x)/N
mir.stderr <- sqrt(mir*(1-mir)/N)
z <- qnorm(1-alpha/2)
c(p=mir,lower=max(0, mir - z * mir.stderr), upper = min(1, mir + z * mir.stderr),alpha=alpha)
}

pooledbinom.mir <-
function(x,m,n=rep(1,length(x)))
{
sum(x)/sum(m*n)
}

pooledbinom.mle <-
function(x, m, n = rep(1., length(x)), tol = 1e-008)
{
#
# This is the implementation using Newton-Raphson, as given
# in the Walter, Hildreth, Beaty paper, Am. J. Epi., 1980
#
if(length(m) == 1.) m <- rep(m, length(x)) else if(length(m) != length(x))
stop("\n ... x and m must have same length if length(m) > 1")
if(any(x > n))
stop("x elements must be <= n elements")
if(all(x == 0.))
return(0.)
if(sum(x) == sum(n)) return(1)
p.new <- 1 - (1 - sum(x)/sum(n))^(1/mean(m)) # starting value
done <- 0
N <- sum(n * m)
while(!done) {
p.old <- p.new
p.new <- p.old - (N - sum((m * x)/(1 - (1 - p.old)^m)))/
sum((m^2 * x * (1 - p.old)^(m - 1))/(1 - (1 - 
p.old)^m)^2)
if(abs(p.new - p.old) < tol)
done <- 1
}
p.new
}

pooledbinom.mle.var <-
function(p, m, n = rep(1, length(m)))
{
if(p > 0 & p < 1)
1/sum((m^2 * n * (1 - p)^(m - 2))/(1 - (1 - p)^m))
else 0 #1/sum(m*n)
}

pooledbinom.mu3 <-
function(p, m, n = rep(1,length(m)))
{
if(p > 0)
sum((m^3 * n * (1-p)^(m-3) * (2 * (1-p)^m - 1))/(1-(1-p)^m)^2)
else 0
}

pooledbinom.score.ci <-
function(x, m, n = rep(1., length(x)), tol = .Machine$double.eps^0.75, alpha = 
0.05)
{
FUN <- function(p0, x, mm, n, alpha)
{
# NOTE: I use the square of the score statistic and chisq
score.p(p0, x, mm, n)^2 * pooledbinom.mle.var(p0, mm, n) - 
qchisq(1 - alpha, 1)
}
# The MLE is just used for a sensible cut-value for the search ...
# it's not needed in the computation
p.hat <- pooledbinom.mle(x, m, n, tol)
if(sum(x) == 0) {
lower.limit <- 0
}
else {
root.brak <- bracket.bounded.root(f=FUN, lower = p.hat/10, lbnd = 
0., upper = p.hat/2, ubnd = p.hat, x = x, mm = m, n = n,
alpha = alpha)
lower.limit <- uniroot(f=FUN, lower = root.brak[1], upper = 
root.brak[2], tol = .Machine$double.eps^0.75, x = x,
mm = m, n = n, alpha = alpha)$root
}
if(sum(x) == sum(n)) {
upper.limit <- 1
}
else {
root.brak <- bracket.bounded.root(f=FUN, lower = (p.hat + 1)/10,
lbnd = p.hat, upper = (p.hat + 1)/2, ubnd = 1., x = x,
mm = m, n = n, alpha = alpha)
upper.limit <- uniroot(f=FUN, lower = root.brak[1], upper = 
root.brak[2], tol = .Machine$double.eps^0.75, x = x,
mm = m, n = n, alpha = alpha)$root
}
c(p = p.hat, lower = lower.limit, upper = upper.limit, alpha = alpha
)
}


pooledbinom.wald.ci <-
function(x, m, n = rep(1, length(x)), tol = 1e-008, alpha = 0.05)
{
if(length(m) == 1)
m <- rep(m, length(x))
p.hat <- pooledbinom.mle(x, m, n, tol)
p.stderr <- sqrt(pooledbinom.mle.var(p.hat, m, n))
z <- qnorm(1 - alpha/2)
c(p = p.hat, lower = max(0, p.hat - z * p.stderr), upper = min(1, p.hat +
z * p.stderr), alpha = alpha)
}




























