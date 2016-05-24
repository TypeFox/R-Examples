ss.power.reg.coef <- function(Rho2.Y_X=NULL, Rho2.Y_X.without.j=NULL, p=NULL, desired.power=.85, alpha.level=.05, Directional=FALSE, beta.j=NULL, sigma.X=NULL, sigma.Y=NULL, Rho2.j_X.without.j=NULL, RHO.XX=NULL, Rho.YX=NULL, which.predictor=NULL, Cohen.f2=NULL, Specified.N=NULL, Print.Progress=FALSE)
{

# Method 1 to define the effect size.
if(!is.null(Rho2.Y_X) & !is.null(Rho2.Y_X.without.j))
{
if(!is.null(RHO.XX) | !is.null(Rho.YX)) stop("Since \'Rho2.Y_X,\' and \'Rho2.Y_X.without.j\', have been specified, do not specify \'RHO.XX\' or \'Rho.YX\'.")
if(!is.null(which.predictor)) stop("Since you have specified \'Rho2.Y_X\' and \'Rho2.Y_X.without.j\', do not specifiy \'which.predictor\'.")
f <- sqrt((Rho2.Y_X - Rho2.Y_X.without.j)/(1-Rho2.Y_X))
}


# Method 2 to define the effect size. 
if(is.null(Rho2.Y_X) & is.null(Rho2.j_X.without.j) & is.null(Rho2.Y_X.without.j) & !is.null(Rho.YX) & !is.null(RHO.XX))
{
if(is.null(which.predictor)) stop("Since you have specified \'RHO.XX\' and \'Rho.YX\', you must also specifiy the \'which.predictor\'.")
if(!is.null(p)) stop("There is no need to specify \'p\' in this situation.")
if((dim(RHO.XX)[1]!=dim(RHO.XX)[2]) | (dim(RHO.XX)[1]!=length(Rho.YX))) stop("There is a problem with \'RHO.XX\' and/or \'Rho.YX\'")
p <- dim(RHO.XX)[1]
Rho2.Y_X <- (Rho.YX%*%solve(RHO.XX)%*%Rho.YX)
Rho2.Y_X.without.j <- (Rho.YX[-which.predictor]%*%solve(RHO.XX[-which.predictor,-which.predictor])%*%Rho.YX[-which.predictor])
f <- sqrt((Rho2.Y_X - Rho2.Y_X.without.j)/(1-Rho2.Y_X))
}

# Method 3 to define the effect size.
if(!is.null(Rho2.Y_X) & !is.null(beta.j) & !is.null(Rho2.j_X.without.j))
{
if(!is.null(which.predictor)) stop("Since you have specified \'Rho2.Y_X\', \'beta.j\', and \'Rho2.j_X.without.j\', do not specifiy \'which.predictor\'.")
if(is.null(sigma.X) | is.null(sigma.Y)) stop("Since you have specified \'Rho2.Y_X\', \'beta.j\', and \'Rho2.j_X.without.j\', you must also specify \'sigma.X\' and \'sigma.Y\'.")
f <- beta.j * sqrt((1-Rho2.j_X.without.j)/(1-Rho2.Y_X)) * (sigma.X/sigma.Y)
}

# Method 4, input the effect size, a la Cohen, directly.
if(!is.null(Cohen.f2))
{
if(!is.null(which.predictor)) stop("You do not need to specify \'which.predictor\'.")
if(!is.null(sigma.X) | !is.null(sigma.Y)) stop("You do not need to specify \'sigma.X\' or \'sigma.Y\'.")
if(is.null(p)) stop("You need to specify \'p\'.")
if(!is.null(RHO.XX) | !is.null(Rho.YX)) stop("Since the effect size, \'Cohen.f2\' was specified directly, you do not need to specify \'RHO.XX\' or \'Rho.YX\'.")
if(!is.null(Rho2.Y_X) | !is.null(Rho2.j_X.without.j) | !is.null(Rho2.Y_X.without.j)) stop("Since the effect size, \'Cohen.f2\' was specified directly, you do not need to specify \'Rho2.Y_X\', \'Rho2.j_X.without.j\', or \'Rho2.Y_X.without.j\'.")
if(Cohen.f2<0) stop("The effect size, \'Cohen.f2\', cannot be negative.")
f <- sqrt(Cohen.f2)
}

# Now that the noncentral t parameter has been defined, contine to sample size determination.
#################################################

if(is.null(Specified.N))
{
# Minimum sample size to start.
N.i <- p+1+1

Dif <- 1
while(Dif > 0)
{
N.i <- N.i + 1

# Given the critical value, which changes with sample size updates, find the noncentrality parameter that leads to
# the critical value having Power*100% of the alternative distribution beyond the critical value.
CV.for.test.of.Null <- ifelse(Directional==FALSE, qt((1-alpha.level/2), df=(N.i-p-1), lower.tail = TRUE, log.p = FALSE), qt((1-alpha.level), df=(N.i-p-1), lower.tail = TRUE, log.p = FALSE))

Actual.Power <- 1 - pt(CV.for.test.of.Null, df=(N.i-p-1), ncp=(sqrt(N.i)*abs(f)), lower.tail = TRUE, log.p = FALSE)

if(Print.Progress==TRUE) cat(c(Current.Power=Actual.Power, Current.NC.t.Parm=(sqrt(N.i)*abs(f)), Current.N=N.i), "\n")

# Algorithm stops once Actual.Power is greater than the desired power. 
Dif <- desired.power - Actual.Power
}
VALUE <- list(Necessary.Sample.Size=N.i, Actual.Power=Actual.Power, Noncentral.t.Parm=(sqrt(N.i)*abs(f)), Effect.Size.NC.t=f)
}

if(!is.null(Specified.N))
{
CV.Specified.N <- ifelse(Directional==FALSE, qt((1-alpha.level/2), df=(Specified.N-p-1), lower.tail = TRUE, log.p = FALSE), qt((1-alpha.level), df=(Specified.N-p-1), lower.tail = TRUE, log.p = FALSE))
Actual.Power.Specified.N <- 1 - pt(CV.Specified.N, df=(Specified.N-p-1), ncp=(sqrt(Specified.N)*abs(f)), lower.tail = TRUE, log.p = FALSE)
VALUE <- list(Sample.Size=Specified.N, Actual.Power=Actual.Power.Specified.N, Noncentral.t.Parm=(sqrt(Specified.N)*abs(f)), Effect.Size.NC.t=f)
}


return(VALUE)
}
