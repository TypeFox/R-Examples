"conf.limits.ncf" <-
function(F.value=NULL, conf.level=.95, df.1=NULL, df.2=NULL, alpha.lower=NULL, alpha.upper=NULL, tol=1e-9, Jumping.Prop=.10)
{
if(Jumping.Prop <=0 | Jumping.Prop >= 1) stop("The Jumping Proportion (\'Jumping.Prop\') must be between zero and one.")
if(is.null(F.value)) stop("Your \'F.value\' is not correctly specified.")
if(F.value < 0) stop("Your \'F.value\' is not correctly specified.")
if(is.null(df.1) | is.null(df.2)) stop("You must specify the degrees of freedom (\'df.1\' and \'df.2\').")
if(is.null(alpha.lower) & is.null(alpha.upper) & is.null(conf.level)) stop("You need to specify the confidence interval parameters.")
if((!is.null(alpha.lower) | !is.null(alpha.upper)) & !is.null(conf.level)) stop("You must specify only one method of defining the confidence limits.")

if(!is.null(conf.level))
{
if(conf.level >=1 | conf.level <= 0) stop("Your confidence level (\'conf.level\') must be between 0 and 1.")
alpha.lower <- alpha.upper <- (1-conf.level)/2
}

if(alpha.lower==0) alpha.lower <- NULL
if(alpha.upper==0) alpha.upper <- NULL

# Critical value for lower tail.
################################################################################################
FAILED <- NULL
if(!is.null(alpha.lower))
{
LL.0 <- qf(p=alpha.lower*.0005, df1=df.1, df2=df.2) # Obtain a lower value by using the central F distribution
Diff <- pf(q=F.value, df1=df.1, df2=df.2, ncp=LL.0) - (1-alpha.lower)

if(pf(q=F.value, df1=df.1, df2=df.2, ncp=LL.0) < (1-alpha.lower))
{
FAILED <- if(pf(q=F.value, df1=df.1, df2=df.2, ncp=0) < 1-alpha.lower) 
LL.0 <- .00000001
if(pf(q=F.value, df1=df.1, df2=df.2, ncp=LL.0) < 1-alpha.lower) FAILED <- TRUE
#if(FAILED==TRUE) warning("The size of the effect combined with the degrees of freedom is too small to determine a lower confidence limit for the \'alpha.lower\' (or the (1/2)(1-\'conf.level\') symmetric) value specified (set to zero).", call. = FALSE)
}

if(is.null(FAILED))
{
LL.1 <- LL.2 <- LL.0 # Define both in case there is no need for the while loop (LL.2 is overwritten later if the while loop is used).

while(Diff > tol) # Find a value that is too small and one that is too big.
{
LL.2 <- LL.1*(1+Jumping.Prop) 
Diff <- pf(q=F.value, df1=df.1, df2=df.2, ncp=LL.2) - (1-alpha.lower)
LL.1 <- LL.2
}
LL.1 <- LL.2/(1+Jumping.Prop) # Produces the value directly before failure (a Lambda value that is too small.)

LL.Bounds <- c(LL.1, (LL.1+LL.2)/2, LL.2) # The middle value is in the middle.

Diff <- pf(q=F.value, df1=df.1, df2=df.2, ncp=LL.Bounds[2])-(1-alpha.lower)
while(abs(Diff) > tol) # Run the while loop to home in on the value satisfying the conditions (i.e., the lower limit).
{
Diff.1 <- pf(q=F.value, df1=df.1, df2=df.2, ncp=LL.Bounds[1])-(1-alpha.lower) > tol
Diff.2 <- pf(q=F.value, df1=df.1, df2=df.2, ncp=LL.Bounds[2])-(1-alpha.lower) > tol
Diff.3 <- pf(q=F.value, df1=df.1, df2=df.2, ncp=LL.Bounds[3])-(1-alpha.lower) > tol

if(Diff.1==TRUE & Diff.2==TRUE & Diff.3==FALSE)
{
LL.Bounds <- c(LL.Bounds[2], (LL.Bounds[2]+LL.Bounds[3])/2, LL.Bounds[3])
}

if(Diff.1==TRUE & Diff.2==FALSE & Diff.3==FALSE)
{
LL.Bounds <- c(LL.Bounds[1], (LL.Bounds[1]+LL.Bounds[2])/2, LL.Bounds[2])
}

Diff <- pf(q=F.value, df1=df.1, df2=df.2, ncp=LL.Bounds[2])-(1-alpha.lower)

}
LL <- LL.Bounds[2] # Confidence limit.
}
}
if(!is.null(FAILED)) LL <- NA
################################################################################################

# Critical value for upper tail.
################################################################################################
if(!is.null(alpha.upper))
{
FAILED.Up <- NULL
UL.0 <- qf(p=1-alpha.upper*.0005, df1=df.1, df2=df.2)
Diff <- pf(q=F.value, df1=df.1, df2=df.2, ncp=UL.0)-alpha.upper

if(Diff < 0) UL.0 <- .00000001

Diff <- pf(q=F.value, df1=df.1, df2=df.2, ncp=UL.0)-alpha.upper
if(Diff < 0) 
{
FAILED.Up <- TRUE
# warning("The size of the effect combined with the degrees of freedom is too small to determine an upper confidence limit for the \'alpha.upper\' (or (1/2)(1-\'conf.level\') symmetric) value specified.", call. = FALSE)
}

if(is.null(FAILED.Up))
{
UL.1 <- UL.2 <- UL.0
while(Diff > tol)
{
UL.2 <- UL.1*(1+Jumping.Prop) 
Diff <-  pf(q=F.value, df1=df.1, df2=df.2, ncp=UL.2) - alpha.upper
UL.1 <- UL.2
}
UL.1 <- UL.2/(1+Jumping.Prop) 

UL.Bounds <- c(UL.1, (UL.1+UL.2)/2, UL.2)

Diff <- pf(q=F.value, df1=df.1, df2=df.2, ncp=UL.Bounds[2])-alpha.upper
while(abs(Diff) > tol)
{
Diff.1 <- pf(q=F.value, df1=df.1, df2=df.2, ncp=UL.Bounds[1])-alpha.upper > tol
Diff.2 <- pf(q=F.value, df1=df.1, df2=df.2, ncp=UL.Bounds[2])-alpha.upper > tol
Diff.3 <- pf(q=F.value, df1=df.1, df2=df.2, ncp=UL.Bounds[3])-alpha.upper > tol

if(Diff.1==TRUE & Diff.2==TRUE & Diff.3==FALSE)
{
UL.Bounds <- c(UL.Bounds[2], (UL.Bounds[2]+UL.Bounds[3])/2, UL.Bounds[3])
}

if(Diff.1==TRUE & Diff.2==FALSE & Diff.3==FALSE)
{
UL.Bounds <- c(UL.Bounds[1], (UL.Bounds[1]+UL.Bounds[2])/2, UL.Bounds[2])
}

Diff <- pf(q=F.value, df1=df.1, df2=df.2, ncp=UL.Bounds[2])-alpha.upper

}
UL <- UL.Bounds[2] # Confidence limit.
}
if(!is.null(FAILED.Up)) UL <- NA
}
################################################################################################
if(!is.null(alpha.lower) & !is.null(alpha.upper)) return(list(Lower.Limit=LL, Prob.Less.Lower=1-pf(q=F.value, df1=df.1, df2=df.2, ncp=LL), Upper.Limit=UL, Prob.Greater.Upper=pf(q=F.value, df1=df.1, df2=df.2, ncp=UL)))
if(is.null(alpha.lower) & !is.null(alpha.upper)) return(list(Upper.Limit=UL, Prob.Greater.Upper=pf(q=F.value, df1=df.1, df2=df.2, ncp=UL)))
if(!is.null(alpha.lower) & is.null(alpha.upper)) return(list(Lower.Limit=LL, Prob.Less.Lower=1-pf(q=F.value, df1=df.1, df2=df.2, ncp=LL)))
}
