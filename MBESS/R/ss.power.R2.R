ss.power.R2 <- function(Population.R2=NULL, alpha.level=.05, desired.power=.85, p, Specified.N=NULL, Cohen.f2=NULL, Null.R2=0, Print.Progress=FALSE, ...)
{
if(alpha.level > 1 | alpha.level < 0) stop("\'alpha.level\' has been specified incorrectly.")
if(Null.R2 > 1 | Null.R2 < 0) stop("\'Null.R2\' has been specified incorrectly.")   

if(!is.null(Population.R2))
{
if(!is.null(Cohen.f2)) stop("Since you specified \'Population.R2\', do not specify \'Cohen.f2\'.")
if(Population.R2 > 1 | Population.R2 < 0) stop("\'Population.R2\' has been specified incorrectly.")
if(desired.power > 1 | desired.power < 0) stop("\'desired.power\' has been specified incorrectly.")
f.Squared <- (Population.R2-Null.R2)/(1-Population.R2)
}

if(!is.null(Cohen.f2)) f.Squared <- Cohen.f2

if(is.null(Specified.N))
{
# Starting value for the iterative N.
N.i <- p + 1 + 1

Dif <- 1
while(Dif > 0)
{
N.i <- N.i + 1

# Given the critical value, which changes with sample size updates, find the noncentrality parameter that leads to
# the critical value having Power*100% of the alternative distribution beyond the critical value.
CV.for.test.of.Null <- qf((1-alpha.level), df1=p, df2=(N.i-p-1), lower.tail = TRUE, log.p = FALSE)

Actual.Power <- 1 - pf(CV.for.test.of.Null, df1=p, df2=(N.i-p-1), ncp=(N.i*f.Squared), lower.tail = TRUE, log.p = FALSE)

if(Print.Progress==TRUE) cat(c(Current.Power=Actual.Power, Current.NC.F.Parm=(N.i*f.Squared), Current.N=N.i), "\n")

# Algorithm stops once Actual.Power is greater than the desired power. 
Dif <- desired.power - Actual.Power
}
VALUE <- list(Necessary.Sample.Size=N.i, Actual.Power=Actual.Power, Noncentral.F.Parm=(N.i*f.Squared), Effect.Size=f.Squared)
}

if(!is.null(Specified.N))
{
CV.Specified.N <- qf((1-alpha.level), df1=p, df2=(Specified.N-p-1), lower.tail = TRUE, log.p = FALSE)
Actual.Power.Specified.N <- 1 - pf(CV.Specified.N, df1=p, df2=(Specified.N-p-1), ncp=(Specified.N*f.Squared), lower.tail = TRUE, log.p = FALSE)
VALUE <- list(Specified.Sample.Size=Specified.N, Actual.Power=Actual.Power.Specified.N, Noncentral.F.Parm=(Specified.N*f.Squared), Effect.Size=f.Squared)
}
return(VALUE)
}
