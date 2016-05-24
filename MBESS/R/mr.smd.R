mr.smd <- function(A, structural.cost, epsilon, d, n, sampling.cost, pilot=FALSE, m0=4, gamma=.49) 
{  
# d current estimate of the standardized mean difference. 
# n current sample size per group.
# A the cost that the researcher would be willing to pay to have a maximum allowable difference of the population standardized mean difference and its estimate of $\epsilon$.
# sampling.cost the cost of collecting an observation; thus the cost is 2n for increasing sample size by 1 in each group. 
# gamma
# m pilot sample size.
# Here we check to see if the current sample size is sufficiently large
  
if(!missing(A))
{
if(!missing(structural.cost)) stop("Because you specified \'A\' directly, you should not also specify \'structural.cost\'.")
if(!missing(epsilon)) stop("Because you specified \'A\' directly, you should not also specify \'epsilon\'.")
if(A<=0) stop("A should be a non-zero positive value")
}  

if(missing(A))
{
if(missing(structural.cost)) stop("Because you did not specificy \'A\' directly, you must specify \'structural.cost\'.")
if(missing(epsilon)) stop("Because you did not specificy \'A\' directly, you must specify \'epsilon\'.")
A <- structural.cost/(epsilon^2) 
if(A<=0) stop("A should be a non-zero positive value")
}    
  
if(pilot==FALSE)
{
Stop <- FALSE

Criterion <- sqrt(A/(2*sampling.cost))*(sqrt(2 + d^2/4)+n^(-gamma))

if(n >= Criterion) Stop <- TRUE
if(n < Criterion) Stop <- FALSE

##### add risk here?
Rk <- A*((1/n + 1/n) + (d^2)/(2*(n + n))) + sampling.cost*(n + n)
Outcome <- rbind(list(Risk=Rk, n1=n, n2=n, d=d, "Is.Satisfied?"=Stop))

}  


if(pilot==TRUE)
{
if(m0 < 4) stop("The value of 'm0' must be 4 or greater.")
Outcome <-c(Pilot.SS=max(m0, ceiling((A/(2*sampling.cost))^(1/(2+2*gamma)))))
}

return(Outcome)
}
