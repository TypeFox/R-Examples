negbin.pow <-
function(n,lambda1,k,disp=1.5,alpha,seed = 200, numsim = 1000,monitor=TRUE,sig=3)
{
#------------------------------------------
# Error check
#------------------------------------------
n.integer <- function(x, tol = .Machine$double.eps) 
{ abs(x - round(x)) < tol }

if(alpha >= 1 || alpha <= 0) stop("Error: alpha must be between 0 and 1")
if( any(k <= 0) ) stop("Error: k <= 0")
if( any(disp <= 0) ) stop("Error: Dispersion <= 0. See documentation for underdispersion specification")
if( mean(n.integer(n)) !=1 || n<= 0) stop("n must be positive integer")


#------------------------------------------
# Random vals from Overdispersed poisson
negbin.disp<-function (n, lambda,disp) 
{
if (disp==1) {rpois(n, lambda)}
else {rnbinom(n, size=(1/(disp-1)), mu=lambda)}
}


#------------------------------------------
l <- list(n=n,lambda1=lambda1,k=k,disp=disp,alpha=alpha,seed = seed,numsim=numsim)
store <- expand.grid(l)

inner.fcn <- function(n,lambda1,k,disp,alpha,seed,numsim) 
{
T_nb <- NULL
set.seed(seed)
for (i in 1:numsim)
{
x <- negbin.disp(n, lambda1,disp)
y <- negbin.disp(n, lambda1,disp)
a <- mean(x); b <- mean(y)
se <- sqrt((a+(disp-1)*a^2)/n + (b+(disp-1)*b^2)/n) 
t <- ( mean(x) - mean(y) ) / se
T_nb[i] <- t
}
T_nb <- na.omit(T_nb)
T_crit <- c(quantile(T_nb,alpha/2),quantile(T_nb,1-alpha/2))

# Computing the same test statistic under alternate hypothesis.
T_nb.alt <- NULL

for (i in 1:numsim)
{
x <- negbin.disp(n, lambda1,disp)
y <- negbin.disp(n, k*lambda1, disp)
a <- mean(x); b <- mean(y)
se <- sqrt((a+(disp-1)*a^2)/n + (b+(disp-1)*b^2)/n) 
t <- ( mean(x) - mean(y) ) / se
T_nb.alt[i] <- t
}

summary(T_nb.alt)
# Computing the Power.
val1 <- length(T_nb.alt[T_nb.alt > T_crit[2]])
val2 <- length(T_nb.alt[T_nb.alt < T_crit[1]]) 
power <- sum(val1,val2)/numsim
}
#------------------------------------------
# Organizing the output
#------------------------------------------
if (monitor == TRUE)
{
pb <- txtProgressBar(min = 0, max = dim(store)[1], style = 3)
for (i in 1:dim(store)[1])
{
power <- mapply(inner.fcn,store[,1],store[,2],store[,3],store[,4],store[,5],store[,6],store[,7])
std.err <- sqrt(power*(1-power)/store[,7])
out <- cbind(store[,1:5],round(power,sig),round(std.err,3))
colnames(out) <- c("N","Mean.Null","Effect.Size","Disp.Par","Type.I.Error","Power","Std.Err")
setTxtProgressBar(pb, i)
}
close(pb)
return(out)
}

if (monitor == FALSE)
{
power <- mapply(inner.fcn,store[,1],store[,2],store[,3],store[,4],store[,5],store[,6],store[,7])
std.err <- sqrt(power*(1-power)/store[,7])
out <- cbind(store[,1:5],round(power,sig),round(std.err,3))
colnames(out) <- c("N","Mean.Null","Effect.Size","Disp.Par","Type.I.Error","Power","Std.Err")
return(out)
}

}
