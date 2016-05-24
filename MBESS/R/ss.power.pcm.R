ss.power.pcm <- function(beta, tau, level.1.variance, frequency, duration, desired.power=NULL, N=NULL, alpha.level=.05, standardized=TRUE, directional=FALSE)
{
if(is.null(desired.power) & is.null(N)) stop("You must specify either 'desired.power' or 'N'.")
if(!is.null(desired.power) & !is.null(N)) stop("You must specify either 'desired.power' or 'N', but not both.")

sum.c2.pm <- function(M, p=1, K.p=1/12)
{
K.p*(factorial(M+p)/factorial(M-p-1))
}


M <- frequency*duration + 1

if(is.null(N))
{
N.i <- 4 # Minimum (total) number of sample size possible to start algorithm.

Dif <- 1
while(Dif > 0)
{
N.i <- N.i + 2 # Plus 2 here because N.i is the total sample size

V <- level.1.variance/sum.c2.pm(M=M)

var.beta <- 4*(tau+V)/N.i

reliability <- tau/(tau+V)

if(standardized==FALSE)
{
Lambda <- (beta^2)/var.beta
}
if(standardized==TRUE)
{
delta <- beta
Lambda <- (N.i*delta^2*reliability)/4
}

# Given the critical value, which changes with sample size updates, find the noncentrality parameter that leads to
# the critical value having Power*100% of the alternative distribution beyond the critical value.
if(directional==FALSE) CV.for.test.of.Null <- qt((1-alpha.level/2), df=(N.i-2), lower.tail = TRUE, log.p = FALSE)
if(directional==TRUE)  CV.for.test.of.Null <- qt((1-alpha.level), df=(N.i-2), lower.tail = TRUE, log.p = FALSE)

Actual.Power <- 1 - pt(CV.for.test.of.Null, df=(N.i-2), ncp=sqrt(Lambda), lower.tail = TRUE, log.p = FALSE)
Dif <- desired.power - Actual.Power
}
}

if(is.null(desired.power))
{
V <- level.1.variance/sum.c2.pm(M=M)

var.beta <- 4*(tau+V)/N

reliability <- tau/(tau+V)

if(standardized==FALSE)
{
Lambda <- (beta^2)/var.beta
}

if(standardized==TRUE)
{
delta <- beta
Lambda <- (N*delta^2*reliability)/4
}
# Standardized

# Given the critical value, which changes with sample size updates, find the noncentrality parameter that leads to
# the critical value having Power*100% of the alternative distribution beyond the critical value.
if(directional==FALSE) CV.for.test.of.Null <- qt((1-alpha.level/2), df=(N-2), lower.tail = TRUE, log.p = FALSE)
if(directional==TRUE)  CV.for.test.of.Null <- qt((1-alpha.level), df=(N-2), lower.tail = TRUE, log.p = FALSE)

Actual.Power <- 1 - pt(CV.for.test.of.Null, df=(N-2), ncp=sqrt(Lambda), lower.tail = TRUE, log.p = FALSE)
N.i <- N # To make the output simpler to print.
}
return(list(
Design.features=list(
Necessary.SS.Control=N.i/2, 
Necessary.SS.Treatment=N.i/2, 
Total.SS=N.i, 
Actual.Power=Actual.Power,
Frequency=frequency,
Duration=duration,
Total.Measurement.Occasions=M),
Parameters=list(
Regression.Coefficient=beta,
Standardized.Regression.Coefficient=beta/sqrt(tau),
Level.1.error.variance=level.1.variance,
true.variance.of.slopes=tau,
error.variance.of.slopes=V,
Reliability=reliability,
Noncentral.t.parameter=sqrt(Lambda))))
}

# Example from Raudenbush and Liu (2001)
ss.power.pcm(beta=-.4, tau=.003, level.1.variance=.0262, frequency=2, duration=2, desired.power=.80, alpha.level=.05, standardized=TRUE, directional=FALSE)
ss.power.pcm(beta=-.4, tau=.003, level.1.variance=.0262, frequency=2, duration=2, N=238, alpha.level=.05, standardized=TRUE, directional=FALSE)

####
# The standardized effect size is obtained as beta/sqrt(tau): -.4/sqrt(.003) = -.0219.
ss.power.pcm(beta=-.0219, tau=.003, level.1.variance=.0262, frequency=2, duration=2, desired.power=.80, alpha.level=.05, standardized=FALSE, directional=FALSE)
ss.power.pcm(beta=-.0219, tau=.003, level.1.variance=.0262, frequency=2, duration=2, N=238, alpha.level=.05, standardized=FALSE, directional=FALSE)

#

