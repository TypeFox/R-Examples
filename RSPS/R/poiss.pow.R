poiss.pow <-
function(n,lambda1,k,alpha = 0.05,seed = 20,numsim=2000,monitor=TRUE,sig=3)
{
#------------------------------------------
# Error check
#------------------------------------------
n.integer <- function(x, tol = .Machine$double.eps) 
{ abs(x - round(x)) < tol }

if(alpha >= 1 || alpha <= 0) stop("Error: alpha must be between 0 and 1")
if( any(k <= 0) ) stop("Error: k <= 0")
if( mean(n.integer(n)) !=1 || n<= 0) stop("n must be positive integer")

#------------------------------------------
l <- list(n=n,lambda1=lambda1,k=k,alpha=alpha,seed = seed,numsim=numsim)
store <- expand.grid(l)

inner.fcn <- function(n,lambda1,k,alpha,seed,numsim)
{
set.seed(seed)
T <- NULL

for (i in 1:numsim)
{
x <- rpois(n,lambda = lambda1)
y <- rpois(n,lambda = lambda1)

se <- sqrt(mean(x)/n + mean(y)/n) 
t <- ( mean(x) - mean(y) ) / se
T[i] <- t
}

T_crit <- c(quantile(T,alpha/2),quantile(T,1-alpha/2))

# Computing the same test statistic under alternate hypothesis.
T_alt <- NULL

for (i in 1:numsim)
{
x <- rpois(n,lambda = lambda1)
y <- rpois(n,lambda = k*lambda1)

se <- sqrt(mean(x)/n + mean(y)/n) 
t <- ( mean(x) - mean(y) ) / se
T_alt[i] <- t
}

# Computing the Power.
val1 <- length(T_alt[T_alt > T_crit[2]])
val2 <- length(T_alt[T_alt < T_crit[1]]) 
power <- sum(val1,val2)/numsim
}
#------------------------------------------
# Organizing Output
#------------------------------------------
if (monitor == TRUE){
pb <- txtProgressBar(min = 0, max = dim(store)[1], style = 3)
for (i in 1:dim(store)[1])
{
power <- mapply(inner.fcn,store[,1],store[,2],store[,3],store[,4],store[,5],store[,6])
Std.Err <- sqrt(power*(1-power)/store[,6])
out <- cbind(store[,1:4],round(power,sig),round(Std.Err,sig))
colnames(out) <- c("N","Mean.Null","Effect.Size","Type.I.Error","Power","Std.Err")
setTxtProgressBar(pb, i)
}
close(pb)
return(out)}

if (monitor == FALSE)
{
power <- mapply(inner.fcn,store[,1],store[,2],store[,3],store[,4],store[,5],store[,6])
Std.Err <- sqrt(power*(1-power)/store[,6])
out <- cbind(store[,1:4],round(power,sig),round(Std.Err,sig))
colnames(out) <- c("N","Mean.Null","Effect.Size","Type.I.Error","Power","Std.Err")
}
return(out)
}
