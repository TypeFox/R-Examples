diffDetect <- function(N,sigma,alpha=0.05, power=0.8, two.tailed=TRUE)
{
#Error Checking
if ((alpha >= 1) || (alpha <= 0) || (power <= 0) || (power >= 1))
        stop("Sorry, the alpha and power must lie within (0,1)")

for (i in 1:length(N))
{
if (N[i] <= 0)
{
        stop("Sorry, the specified values of N must be strictly positive...")
}
}

for (i in 1:length(sigma))
{
if (sigma[i] <=0) 
        stop("Sorry, the specified value of sigma must be strictly positive...")
}

#Initialize Parameters
r <- NULL;

r$n <- N; r$sigma <- sigma; r$alpha <- alpha; r$power <- power; r$two.tailed <- two.tailed;

r$delta <- matrix(0,nrow=length(N), ncol=length(sigma));
#Label rows and columns
colnames(r$delta) <- sigma;
rownames(r$delta) <- N;

#Compute delta for one/two-sided test
if (two.tailed)
{
for (i in 1:length(N))
{

for (j in 1:length(sigma))

{
r$delta[i,j] <- (qnorm(1 - alpha/2) + qnorm(power))*sqrt((2*sigma[j]^2)/N[i]);
}

}

}

if (!two.tailed)
{
for (i in 1:length(N))
{

for (j in 1:length(sigma))

{
r$delta[i,j] <- (qnorm(1 - alpha) + qnorm(power))*sqrt((2*sigma[j]^2)/N[i]);
}

}
}

class(r) <- "diffDetect";
return(r);
}

#Print Method
print.diffDetect <- function(x, ...)
{
cat("The Minimum Detectable Difference between two populations for fixed N (in rows) and sigma (in columns) is: \n \n");
print(x$delta)
}

#Summary Method
summary.diffDetect <- function(object, ...)
{
cat("The Minimum Detectable Difference between two populations for fixed N (in rows) and sigma (in columns) is: \n \n");
print(object$delta)
cat("\n This assumes: \n")
cat("Type I Error Rate (alpha) = ", object$alpha, " (with two.tailed=", object$two.tailed, ") and Power = ", object$power, "\n \n",sep="")
cat("Note: Original N and Sigma Vectors are available as $n and $sigma \n")
}