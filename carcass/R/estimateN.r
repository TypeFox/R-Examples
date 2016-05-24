# Function to estimate the number of killed animals from carcass counts 
# when the detection probability is estimated with uncertainty (95% CI)
# Author: Fraenzi Korner-Nievergelt
# This function wraps the code given in Appendix II of the following paper where you find detailed comments in each step:

# Korner-Nievergelt F, Korner-Nievergelt P, Behr O, Niermann I, Brinkmann R & Hellriegel B (2011)
# A new method to determine bird and bat fatality at wind energy turbines from carcass searches. Wildl. Biol. 17: 1-14
#-------------------------------------------------------------------------------



# R-Code written and tested for R. 3.1.0
#-------------------------------------------------------------------------------

# track of changes
# 25.11.2015, fk: calculation of the shapeparameters of a moved outside the is.finit(f) loop; 
#                 this was a severe bug! Thanks Sergio Cabrera-Cruz for reporting inconsistencies!
# 17.7.2014, fk: integration of persistence probabilities over the time intervals 
#                 to account for continuous arrival of carcasses, argument "arrival" added.
# 8. 3. 2016, fk: Persistence probabilities are constrained so that they cannot 
#                 become exactly 1: sr[sr>0.9999] <- 0.9999. Bug reported by Greg Forcey.

#-------------------------------------------------------------------------------


#Estimating the number of fatalities with a credible interval
# Describe the uncertainty in the estimates for searcher efficiency and carcass persistence probability by a beta-distribution, i.e. transform the 95 % CI into the shape parameters of the beta-distributions.

estimateN <- function(count, p=NA, p.lower=NA, p.upper=NA, f=NA, f.lower=NA, f.upper=NA, 
                        s=NA, s.lower=NA, s.upper=NA, arrival="discrete", a=1, a.lower=1, a.upper=1,
                        pform="korner", d=1, n=NA, J=NA, maxn=1000, nsim=1000, 
                        plot=TRUE, postdist=FALSE, k=1, x=c(1:10)){

# count   number of carcasses found
# p       estimate for detection probability
# p.lower lower limit of 95% CI of detection probability
# p.upper upper limit of 95% CI of detection probability

# either p or f and s have to be provided!!

# f       estimate for searcher efficiency
# f.lower lower limit of 95% CI of searcher efficiency
# f.upper upper limit of 95% CI of searcher efficiency
# s       estimate for persistence probability
# s.lower lower limit of 95% CI of persistence probability
# s.upper upper limit of 95% CI of persistence probability
# arrival assumption of the distribution of arrival times: 
#         if "uniform" it is assumed that carcasses arrive and are removed at a constant rate 
#         during the whole time interval.
#         if "discrete" it is assumed that all carcasses arrive at the same time point and 
#         and persistence probability given for the first day corresponds to the time from
#         the discrete arrival to the first day.
  
# a       estimate for the proportion of killed animals that falled into a searched area
# a.lower lower limit of 95% CI of proportion of killed animals in searched area
# a.upper upper limit of 95% CI of proportion of killed animals in searched area
# d       search interval, the number days between two searches
# n       number of searches
# J       vector with the lengths of the search intervals (if search intervals are irregular)
# pform   formula used to estimate p, one of "korner", "huso", "erickson", "etterson"

# maxn   maximal possible number of animals killed for which the posterior probability is estimated (should not be too high but also not be too low!)
# nsim   number of Monte Carlo simulations
# x      the posterior probability P(true number of dead animals > x) is given as result

  
# Warnings
# Say which parameters are ignored if too many are specified
if (is.finite(p) & is.finite(f))  {
    warning('f, f.lower, f.upper, s, s.lower, s.upper, and d will be ignored')
  }
 
   
if(a.lower<a.upper){
    a.a <- shapeparameter(a, a.lower, a.upper)$a
    a.b <- shapeparameter(a, a.lower, a.upper)$b
}
  
if(is.finite(f)[1]){
 if (f.lower>=f.upper) stop("Something is wrong with the CI of f. If f is known without uncertainty, use posteriorN instead of estimateN")
 if (s.lower>=s.upper) stop("Something is wrong with the CI of s. If s is known without uncertainty, use posteriorN instead of estimateN")
 if (a.lower>a.upper) stop("Something is wrong with the CI of a.")
 if (length(s)==1 & s==1) stop("s=1: no removal. This is very unlikely. The formulas implemented here are not made for this case.")

f.a <- shapeparameter(f, f.lower, f.upper)$a
f.b <- shapeparameter(f, f.lower, f.upper)$b
s.a <- shapeparameter(s, s.lower, s.upper)$a
s.b <- shapeparameter(s, s.lower, s.upper)$b

Npostdist <- numeric(maxn+1) 


# messages concerning s: declaration of assumptions concerning time of carcass arrival and persistence probability
if(length(s)==1&arrival=="discrete"){ 
cat("You assume that all carcasses arrive simultaneously 
once a day right after the time the search takes place (if there is one) 
and that the probability to persist from death to the first day 
(the time the searches normally are performed) is equal to the probability
to persist from one day to the next. To change that, 
use a different persistence probability for the first day 
or choose arrival='uniform' to allow carcasses to arrive continuously.\n")
}

if(arrival=="uniform" & length(s)>1)   cat("You assume constant (continuous) carcass arrival. To change that, use arrival='discrete'")
if(arrival=="uniform" & length(s)==1)  cat("You assume that persistence probability is independent of age and that carcasses arrive constantly (continuously).
To change that, use different persistence probabilities for every day after 
death and/or arrival='discrete'\n")
if(arrival=="discrete" & length(s)>1){
cat("You assume that all carcasses arrive simultaneously and that the first 
persistence probability is the probability to persist from death to 
the first day at the time when searches take place (if there is a search at that day).\n")
}
  





for(i in 1:nsim){

fr <- rbeta(length(f), f.a, f.b)

if(length(s)>1){
  sr <- numeric(length(s))
  sr[1] <- rbeta(1, s.a[1], s.b[1])
  psr <- pbeta(sr[1], s.a[1], s.b[1])
  srelevant <- s>0.001  # to avoid warnings in the qbeta code below, do not use quantiles for s close to zero
  sr[srelevant] <- qbeta(psr, s.a[srelevant], s.b[srelevant])
  sr[!srelevant] <- 0
}

if(length(s)==1) sr <- rbeta(1, s.a, s.b)
sr[sr>0.9999] <- 0.9999

# integration of persistence probability over the interval when constant arrival
if(arrival=="uniform") sr <- integrate.persistence(sr, n=n, d=d)
  
ar <- ifelse(a.lower<a.upper, rbeta(1, a.a, a.b), a)

if(pform=="korner") pr <- pkorner(s=sr, f=fr, d=d, n=n, k=k, search.efficiency.constant=ifelse(k==1, TRUE, FALSE))
if(pform=="huso") pr <- phuso(s=sr, f=fr, d=d)
if(pform=="erickson") pr <- perickson(t.bar=-1/log(sr), f=fr, d=d)  # 
if(pform=="etterson") {
  if(length(J)==1) if(is.na(J)) stop("Please, provide J")
  if(length(sr)==1) sr <- rep(sr, sum(J))
  if(length(fr)==1) fr <- rep(fr, length(J))
  pr <- ettersonEq14v2(s=sr, f=fr, J=J)
}

postNtemp <- posteriorN(nf=count, p=pr*ar, maxN=maxn,  plot=FALSE, dist=TRUE)

Npostdist <- Npostdist + postNtemp$pN    #Sum the posterior densities over all nsim simulations.

} # close loop i

if(pform=="korner") pm <- pkorner(s=s, f=f, d=d, n=n, k=k, search.efficiency.constant=ifelse(k==1, TRUE, FALSE))
if(pform=="huso") pm <- phuso(s=s, f=f, d=d)
if(pform=="erickson") pm <- perickson(t.bar=-1/log(s), f=f, d=d)  # 
if(pform=="etterson") {
  if(length(J)==1) if(is.na(J)) stop("Please, provide J")
  if(length(s)==1) s <- rep(s, sum(J))
  if(length(f)==1) f <- rep(f, length(J))
  pm <- ettersonEq14v2(s=s, f=f, J=J)
}


HT.estimate <- count/(pm*a) # Horwitz-Thompson estimate
} # close sim based on f and s


if(is.finite(p)){
if(p.lower>=p.upper) stop("Something is wrong with the CI of p. If p is known without error, use posteriorN instead of estimateN")
p.a <- shapeparameter(p, p.lower, p.upper)$a
p.b <- shapeparameter(p, p.lower, p.upper)$b
HT.estimate <- count/(p*a) # Horwitz-Thompson estimate

Npostdist <- numeric(maxn+1) 


for(i in 1:nsim){
pr <- rbeta(1, p.a, p.b)
ar <- ifelse(a.lower<a.upper, rbeta(1, a.a, a.b),a)
postNtemp <- posteriorN(nf=count, p=pr*ar, maxN=maxn,  plot=FALSE, dist=TRUE)
Npostdist <- Npostdist + postNtemp$pN    #Sum the posterior densities over all nsim simulations.

} # close loop i
} # close sim based on p


Npostdist.sc <- Npostdist/nsim
indexLower <- cumsum(Npostdist.sc) < 0.025
indexMedian <- cumsum(Npostdist.sc) < 0.5 
indexUpper <- cumsum(Npostdist.sc) < 0.975
lower <- min(c(0:maxn)[!indexLower])  
estimate.median <- min(c(0:maxn)[!indexMedian])  
#estimate.mean <- sum(c(0:maxn)*Npostdist.sc)
upper <- min(c(0:maxn)[!indexUpper])

# post probability that true mortality is larger x
plarger <- 1-cumsum(Npostdist.sc)[is.element(c(0:maxn),x)]
names(plarger) <- paste0("x=", x)

if(plot) plot(0:maxn, Npostdist.sc , type="h", lwd=5, lend="butt", xlab="Number of fatalities", ylab="Posterior density", xlim=c(0,50))
if(!postdist) {
result <- list(estimate=estimate.median, lower=lower, upper=upper, 
               HT.estimate=HT.estimate,
               P.true.larger.x=plarger)
}
if(postdist) {
result <- list(estimate=estimate.median, lower=lower, upper=upper, HT.estimate=HT.estimate, 
               postdist=Npostdist.sc, P.true.larger.x=plarger)
}

return(result)
}






