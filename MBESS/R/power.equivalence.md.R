`power.equivalence.md` <-
function(alpha,logscale,ltheta1,ltheta2,ldiff,sigma,n,nu){

# power of the two one-sided tests procedure in bioquivalence on either the regular or the logarithmic scale

# Regular scale:
# K.F. Phillips
# Power of the Two One-Sided Tests Procedure in Bioquivalence
# Journal of Pharmacokinetics and Biopharmaceutics, Vol. 18, No. 2, 1990

# Logarithmic scale:
# E. Diletti, D. Hauschke, V.W. Steinijans
# Sample size determination of bioequivalence assessment by means of confidence intervals
# International Journal of Clinical Pharmacology, Therapy and Toxicology
# Vol. 29, No, 1 - 1991, p1-8


# parameters on the ordinary scale are:
# alpha:       alpha level for the 2 t-tests (usually 0.05).  Confidence interval for full test is at level 1 - 2 alpha 
# logscale:    FALSE=use regular scale; TRUE=use logarithmic scale
# ltheta1:     lower limit of equivalence interval
# ltheta2:     upper limit of equivalence interval
# ldiff:       true difference or ratio (log scale) in treatment means
# sigma:       sqrt(error variance); ie root MSE from ANOVA (eg, 2 for 2-period crossover)
# n:           number of subjects per treatment
# nu:          degrees of freedom for sigma

theta1 = ltheta1
theta2 = ltheta2
diff  =  ldiff

if (logscale){
theta1 = log(ltheta1)
theta2 = log(ltheta2)
diff  =  log(ldiff)}

power_t <- qt(1.0-alpha,nu,lower.tail=TRUE, log.p = FALSE)                                                                                                             
power_l <- (theta2 - theta1)/(2.0*power_t*sqrt(2.0/n)) -.000001  

#  At sigma=.15, diff=1.1275, n=40, nu=38 the power was underestimated by about 3%
#  possibly due to a problem in function  integrate.
#  Decreasing the upper limit very slightly (.000001) seems to fix the problem

power<-integrate(power.density.equivalence.md, lower=0,upper=power_l,alpha=alpha,theta1= theta1,theta2= theta2,diff=diff,sigma=sigma,n=n,nu=nu,subdivisions=10000,
      rel.tol = .Machine$double.eps^0.25, abs.tol = .Machine$double.eps^0.25, stop.on.error = TRUE, keep.xy = FALSE, aux = NULL)[[1]]

return(power)}
