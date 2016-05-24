n4propsEB <- function(ICC, varICC=0, from=0, to,  pe,pc,m, iter=1000, alpha=0.05, power=0.8, two.tailed=TRUE, digits=3, plot=TRUE)
{

if ((alpha >= 1) || (alpha <= 0) || (power <= 0) || (power >= 1))
        stop("Sorry, the alpha and power must lie within (0,1)")

if ((pe >= 1) || (pe <= 0) || (pc <= 0) || (pc >= 1))
        stop("Sorry, the prior proportions must lie within (0,1)")

for (i in 1:length(ICC))
{
if (ICC[i] <= 0)
        stop("Sorry, the ICC must lie within (0,1)")
}

if (m <=1) 
        stop("Sorry, the (average) cluster size, m, should be greater than one...")

if (to < from) 
        stop("From and To form the range of the estimated density for the ICC...")

#If m is a decimal, round up to generate a more conservative sample size.
m <- ceiling(m);

#Initialize Parameters
r <- NULL
r$pe <- pe; r$pc <- pc; r$ICC <- ICC; r$varICC <- varICC; r$m <- m; 
r$alpha <- alpha; r$power <- power; r$two.tailed <- two.tailed; 
r$digits <- digits; r$from <- from; r$to <- to;


#One or two-tailed tests
if (two.tailed)
{
ZA <- -qnorm(alpha/2);
ZB <- -qnorm(1 - power);
}

else {
ZA <- -qnorm(alpha);
ZB <- -qnorm(1 - power);
}

#Initialization of Results Vectors
r$ResRho <- NULL;
r$ResK <- NULL;

#Compute Density Function; n must be much larger than iter.


if (sum(varICC) != 0)
{
dens <- density(ICC, n=2^16, from=from, to=to, weights = 1/(varICC + var(varICC))/sum(1/(varICC+var(varICC))))
}
else
{
dens <- density(ICC, n=2^16, from=from, to=to)
}

rhoVector <- sample(dens$x, size=iter, prob=dens$y)

#Computational Loop

for (i in 1:iter)

{
rho <- rhoVector[i]

k <- (((ZA + ZB)^2)*(pe*(1-pe) + pc*(1-pc))*(1 + (m - 1)*ICC))/(m*(pe - pc)^2);

r$ResRho <- append(r$ResRho, rho);
r$ResK <- append(r$ResK, k);
}

if (plot)
{
hist(ICC, freq = FALSE, main="Histogram of Values of ICC and
Empirical Density", ylab="Density", xlab="ICC Estimates",
xlim=c(from, to), ylim=c(0,25));
par(new=TRUE);
plot(dens, xlim=c(from,to), ylim=c(0,25), main="", xlab="");
}


class(r) <- "n4propsEB";
return(r);

}

#Print Method
print.n4propsEB <- function(x, ...)
{
cat("Simulation of the Empirical Density suggests that appropriate quantiles \n")
cat("for the number of clusters to be randomized in each group are: \n")
print(round(quantile(x$ResK, probs=c(0,0.25,0.5,0.75, 1.0)),digits=x$digits))
cat("With ICC quantiles: \n ")
print(round(quantile(x$ResRho, probs=c(0,0.25,0.5,0.75, 1.0)),digits=x$digits))
}

#Summary Method
summary.n4propsEB <- function(object, ...)
{
cat("Simulation Based Sample Size Estimation (Empirical Density) to Compare Means of Two Populations", "\n \n")
cat("Assuming:", "\n")
cat("Treatment Rate = ", object$pe, "\n")
cat("Control Rate = ", object$pc, "\n");
cat("Cluster Size (average) = ", object$m, "\n");
cat("ICCs = ", object$ICC, "\n");
cat("Variance of ICCs = ", object$varICC, "\n");
cat("Type I Error Rate (alpha) = ", object$alpha, " and Power = ", object$power, "\n \n",sep="")

cat("Simulation of the Empirical Density suggests that appropriate quantiles \n")
cat("for the number of clusters to be randomized in each group are: \n")
print(round(quantile(object$ResK, probs=c(0,0.1, 0.25,0.5,0.75, 0.9, 1.0)),digits=object$digits))
cat("With ICC quantiles: \n")
print(round(quantile(object$ResRho, probs=c(0,0.1, 0.25,0.5,0.75, 0.9, 1.0)),digits=object$digits))


}