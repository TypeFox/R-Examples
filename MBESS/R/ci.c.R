ci.c <- function(means=NULL, s.anova=NULL, c.weights=NULL, n=NULL, N=NULL,
Psi=NULL, conf.level=.95, alpha.lower=NULL, alpha.upper=NULL, df.error=NULL, ...)
{


if(length(n)==1) 
{
n <- rep(n, length(means))
}

if(length(n)!=length(c.weights)) stop("The lengths of 'n' and 'c.weights' differ, which should not be the case.")

if(sum(c.weights[c.weights>0])>1) stop("Please use fractions to specify the contrast weights")
  
if(!(sum(c.weights)==0)) stop("The sum of the contrast weights ('c.weights') should equal zero.")

part.of.se <- sqrt(sum((c.weights^2)/n))


if(!is.null(Psi))
{
if(!is.null(means)) stop("Since the contrast effect ('Psi') was specified, you should not specify the vector of means ('means').")
if(is.null(s.anova)) stop("You must specify the standard deviation of the errors (i.e., the square root of the error variance).")
if(is.null(n)) stop("You must specify the vector per group/level sample size ('n').")
if(is.null(c.weights)) stop("You must specify the vector of contrast weights ('c.weights').")
}


if(!is.null(means))
{
Psi <- sum(c.weights*means)
}

if(is.null(alpha.lower) & is.null(alpha.upper))
{
alpha.lower <- (1-conf.level)/2
alpha.upper <- (1-conf.level)/2
}

if(is.null(N)&& is.null(n)) stop("You must specify the either total sample size ('N'), or sample sizes per group('n').")
if(is.null(N) && !is.null(n)) N <- sum(n)
if(is.null(df.error)) df.2 <- N - length(means)

CV.up <- qt(1-alpha.upper, df=df.2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
CV.low <- qt(alpha.lower, df=df.2, ncp = 0, lower.tail = TRUE, log.p = FALSE) 

Result <- list(Lower.Conf.Limit.Contrast = Psi + CV.low*part.of.se*s.anova,
Contrast = Psi, Upper.Conf.Limit.Contrast = Psi + CV.up*part.of.se*s.anova)
        
print(paste("The", 1 - (alpha.lower + alpha.upper), "confidence limits for the contrast are given as:"))
return(Result)
}
