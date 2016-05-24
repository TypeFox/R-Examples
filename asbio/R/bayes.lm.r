bayes.lm <- function(Y, X, model = "anova", length = 1000, cred = .95)
{
n <- length(Y); p <- ncol(X)
beta.hat <- solve(qr(X, LAPACK=TRUE), Y)
V.beta <- solve(t(X)%*%X)
MSE<-(t(Y)%*%(Y)-t(beta.hat)%*%t(X)%*%Y)/(n - p)
sigma.sq <- rinvchisq(length, n - p, MSE)

beta.f.x <- matrix(ncol = length(beta.hat), nrow = length)
for(i in 1 : length){
    beta.f.x[i,] <- rmvnorm(1, beta.hat, V.beta *sigma.sq[i])
    }

alpha <- 1 - cred
hi <- 1 - (alpha/2)
low <- 1 - hi
 
posterior <- data.frame(cbind(sigma.sq, beta.f.x))
lower <- apply(posterior, 2, function(x) quantile(x, low))
upper <- apply(posterior, 2, function(x) quantile(x, hi))
med <- apply(posterior, 2, median)

if(model == "reg") rnames <- c("sigma.sq", paste("beta", seq(0, (ncol(X) -1)), sep = ""))
if(model == "anova") rnames <- c("sigma.sq", "mu", paste("alpha", seq(1, (ncol(X) -1)), sep = ""))
cnames <- c(paste(low * 100, "%", sep =""), "Median", paste(hi * 100, "%", sep =""))
head <- c("Bayesian linear model with standard uniform priors")

summary <- data.frame(cbind(lower, med, upper)) 
names(summary) <- cnames
rownames(summary) <- rnames

res <- list(summary = summary, head = head, posterior = posterior)
class(res) <- "blm"
res
}  

print.blm <- function(x, ...){
cat("\n")
cat(x$head, "\n\n")
print(x$summary)
cat("\n")
invisible(x)
}   