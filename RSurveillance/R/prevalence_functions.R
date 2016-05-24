#############################################################
# functions for estimating true and apparent prevalence
##############################################################

#################################################################
# apparent prevalence
#################################################################
##' Sample size for apparent prevalence
##' @description Calculates sample size for estimating apparent 
##' prevalence (simple proportion)
##' @param p expected proportion, scalar or vector of values
##' @param precision absolute precision, +/- proportion equivalent to
##' half the width of the desired confidence interval, scalar or vector of values,
##' note: at least one of p and precision must be a scalar
##' @param conf level of confidence required, default = 0.95 (scalar)
##' @return a vector of sample sizes
##' @keywords methods
##' @export
##' @examples
##' # examples of n.ap
##' n.ap(0.5, 0.1)
##' n.ap(0.5, 0.1, conf=0.99)
##' n.ap(seq(0.1, 0.5, by = 0.1), 0.05)
##' n.ap(0.2, c(0.01, 0.02, 0.05, 0.1))
n.ap<- function(p, precision, conf=0.95) {
  tails<- 2
  z.conf<- qnorm(1 - (1 - conf)/tails, 0, 1)
  # calculate n
  n<- ceiling(z.conf^2*p*(1-p)/(precision^2))
    return(n)
}


##' Agresti-Coull confidence limits
##' @description Calculates Agresti-Coull confidence limits for 
##' a simple proportion (apparent prevalence)
##' @param x number of positives in sample
##' @param n sample size, note: either x or n can be a vector, 
##' but at least one must be scalar
##' @param conf level of confidence required, default 0.95 (scalar)
##' @return a dataframe with 6 columns, x, n, proportion, lower confidence limit,
##' upper confidence limit, confidence level and CI method
##' @keywords methods
##' @export
##' @examples 
##' # test binom.agresti
##' binom.agresti(25, 200)
##' binom.agresti(seq(10, 100, 10), 200)
##' binom.agresti(50, seq(100, 1000, 100))
binom.agresti<- function(x, n, conf=0.95) {
  # agresti-coull
  tails<- 2
  z.conf<- qnorm(1 - (1 - conf)/tails, 0, 1)
  n.ac<- n + z.conf^2
  x.ac<- x + z.conf^2/2
  p.ac<- x.ac/n.ac
  q.ac<- 1 - p.ac
  lc<- p.ac - z.conf*(p.ac*q.ac)^0.5 * n.ac^-0.5
  uc<- p.ac + z.conf*(p.ac*q.ac)^0.5 * n.ac^-0.5
  return(data.frame(x=x, n=n, proportion=p.ac, lower=lc,
               upper=uc, conf.level=conf, method="agresti-coull"))
}


##' Jeffreys confidence limits
##' @description Calculates Jeffreys confidence limits for 
##' a simple proportion (apparent prevalence)
##' @param x number of positives in sample
##' @param n sample size, note: either x or n can be a vector, 
##' but at least one must be scalar
##' @param conf level of confidence required, default = 0.95 (scalar)
##' @return a dataframe with 6 columns, x, n, proportion, lower confidence limit,
##' upper confidence limit, confidence level and CI method
##' @keywords methods
##' @export
##' @examples 
##' # test binom.jeffreys
##' binom.jeffreys(25, 200)
##' binom.jeffreys(seq(10, 100, 10), 200)
##' binom.jeffreys(50, seq(100, 1000, 100))
binom.jeffreys<- function(x, n, conf=0.95) {
  # jeffreys interval
  tails<- 2
  lc<- qbeta((1 - conf)/tails, x+0.5, n-x+0.5)
  uc<- qbeta(1 - (1 - conf)/tails, x+0.5, n-x+0.5)
  p<- x/n
  return(data.frame(x=x, n=n, proportion=p,
               lower=lc, upper=uc, conf.level=conf, method="jeffreys"))
}


##' Clopper-Pearson exact confidence limits
##' @description Calculates Clopper-Pearson exact binomial confidence limits for 
##' a simple proportion (apparent prevalence)
##' @param x number of positives in sample
##' @param n sample size, note: either x or n can be a vector, 
##' but at least one must be scalar
##' @param conf level of confidence required, default = 0.95 (scalar)
##' @return a dataframe with 6 columns, x, n, proportion, lower confidence limit,
##' upper confidence limit, confidence level and CI method
##' @keywords methods
##' @export
##' @examples 
##' # test binom.cp
##' binom.cp(25, 200)
##' binom.cp(seq(10, 100, 10), 200)
##' binom.cp(50, seq(100, 1000, 100))
binom.cp<- function(x, n, conf=0.95) {
  # clopper-pearson exact interval
  tails<- 2
  lc<- qbeta((1 - conf)/tails, x, n-x+1)
  uc<- qbeta(1 - (1 - conf)/tails, x+1, n-x)
  p<- x/n
  return(data.frame(x=x, n=n, proportion=p, lower=lc,
               upper=uc, conf.level=conf, method="clopper-pearson"))
}


##' Apparent prevalence
##' @description Estimates apparent prevalence and confidence limits for
##' given sample size and result, assuming representative sampling
##' @param x number of positives in sample
##' @param n sample size, note: either x or n can be a vector, 
##' but at least one must be scalar
##' @param type method for estimating CI, one of c("normal", "exact", "wilson", "jeffreys", "agresti-coull", "all"),
##' default = "wilson"
##' @param conf level of confidence required, default = 0.95 (scalar)
##' @return either 1) if type = "all", a list with 5 elements, each element
##' a matrix with 6 columns, x, n, proportion, lower confidence limit,
##' upper confidence limit, confidence level and CI method; or
##' 2) a matrix of results for the chosen method
##' @keywords methods
##' @export
##' @examples 
##' # examples for ap function
##' n<- 200
##' x<- 25
##' conf<- 0.95
##' ap(x, n)
##' ap(seq(10, 100, 10), 200, type = "agresti")
##' ap(seq(10, 100, 10), 200, type = "all")
ap<- function(x, n, type = "wilson", conf = 0.95) {
  types<- c("normal", "clopper-pearson", "wilson", "agresti-coull", "jeffreys", "all")
#  require(epitools)
  ap.exact<- binom.cp(x, n, conf)
  ap.normal<- epitools::binom.approx(x, n, conf)
  ap.wilson<- epitools::binom.wilson(x, n, conf)
  ap.jeffreys<- binom.jeffreys(x, n, conf)
  ap.agresti<- binom.agresti(x, n, conf)
  ap.all<- list(normal=ap.normal, "exact"=ap.exact, "wilson"=ap.wilson, "jeffreys"=ap.jeffreys, "agresti" = ap.agresti)
  if (type == "all") {
    return(ap.all)
  } else {
    return(ap.all[[type]])
  }
}

############################################################
# estimated true prevalence
#############################################################
##' Sample size for true prevalence
##' @description Calculates sample size for estimating true prevalence 
##' using normal approximation
##' @param p estimated true prevalence (scalar or vector)
##' @param se test sensitivity (scalar or vector)
##' @param sp test specificity (scalar or vector) 
##' @param precision absolute precision, +/- proportion equal to
##' half the width of the desired confidence interval (scalar or vector)
##' @param conf desired level of confidence for CI, default = 0.95 (scalar or vector)
##' @return a vector of sample sizes
##' @keywords methods
##' @export
##' @examples 
##' # examples for n.tp
##' n.tp(0.1, 0.9, 0.99, 0.05)
##' n.tp(0.1, 0.9, 0.99, 0.05, conf = 0.99)
##' n.tp(c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5), 0.9, 0.99, 0.05)
##' n.tp(0.5, 0.9, 0.99, c(0.01, 0.02, 0.05, 0.1, 0.2))
n.tp<- function(p, se, sp, precision, conf=0.95) {
  tails<- 2
  z.conf<- qnorm(1 - (1 - conf)/tails, 0, 1)
  n<- ceiling((z.conf/precision)^2*(se*p + (1 - sp)*(1 - p))*(1 - se*p - (1 - sp)*(1 - p))/(se + sp - 1)^2)
  return(n)
}

##' Standard deviation of true prevalence estimate
##' @description Calculates the standard deviation of true prevalence estimate
##' assuming se and sp known exactly, used to calculate normal approximation CI for estimate
##' @param x number of positive results in sample (scalar or vector)
##' @param n sample size (scalar or vector)
##' @param se test sensitivity (scalar or vector)
##' @param sp test specificity (scalar or vector)
##' @return vector of standard deviation values for true prevalence estimates
##' @keywords methods
##' @export
##' @examples 
##' # example of sd.tp
##' sd.tp(1:10, 20, 0.9, 0.99)
sd.tp<- function(x, n, se, sp) {
  ap<- x/n
  var.tp<- ap*(1-ap)/(n*(se + sp - 1)^2)
  return(sqrt(var.tp))
} # end of sd.tp function


##' Normal approximation confidence limits for true prevalence
##' @description Estimates true prevalence and confidence limits for 
##' estimates based on normal approximation
##' @param x number of positive results in sample (scalar or vector)
##' @param n sample size (scalar or vector)
##' @param se test unit sensitivity (scalar or vector)
##' @param sp test unit specificity (scalar or vector)
##' @param conf desired level of confidence for CI, default = 0.95 (scalar or vector)
##' @return list with 2 elements, a matrix of apparent prevalence and wilson lower and upper confidence limits
##' and a matrix of true prevalence and normal approximation lower and upper confidence limits
##' @keywords methods
##' @export
##' @examples 
##' # examples for tp.normal
##' tp.normal(25, 120, 0.9, 0.99)
##' tp.normal(seq(5, 25, by=5), 120, 0.9, 0.99)
tp.normal<- function(x, n, se, sp, conf=0.95) {
#  require(epitools)
  tails<- 2
  z.conf<- qnorm(1 - (1 - conf)/tails, 0, 1)
  wilson.ci<- epitools::binom.wilson(x, n, conf)
  ap<- wilson.ci$proportion
  tp<- (ap + sp - 1)/(se + sp - 1)
  tp.ci<- array(0, dim = c(length(tp), 2))
  for (i in 1:length(tp)) tp.ci[i,]<- tp[i] + c(-1, 1)* z.conf*sd.tp(ap[i], n, se, sp)
  ap<- cbind(est=ap, lower=wilson.ci$lower, upper=wilson.ci$upper)
  tp<- cbind(est=tp, lower=tp.ci[,1], upper=tp.ci[,2])
  return(list(ap=ap, tp=tp))
}


##' True prevalence
##' @description Estimates true prevalence and confidence limits for
##' given sample size and result, according to specified method
##' @param x number of positive units (scalar)
##' @param n sample size (no. units sampled) (scalar)
##' @param se test sensitivity (scalar)
##' @param sp test specificity (scalar)
##' @param type method for estimating CI, one of c("normal", "c-p", "sterne", "blaker", "wilson", "all")
##' @param conf desired level of confidence for CI, default = 0.95 (scalar)
##' @return list with 2 elements, a matrix of apparent prevalence and 
##'   lower and upper confidence limits
##'   and a matrix of true prevalence and lower and upper 
##'   confidence limits using the chosen method(s)
##' @keywords methods
##' @export
##' @examples 
##' # examples for tp
##' x<- 20
##' n<- 120
##' se<- 0.9
##' sp<- 0.99
##' conf<- 0.95
##' tp(x, n, se, sp, "all")
##' tp(x, n, se, sp, "c-p")
##' tp(x, n, 0.95, 0.9, "c-p")
tp<- function(x, n, se, sp, type = "blaker", conf=0.95) {
#  require(epiR)
  tp.cp<-epiR::epi.prev(x, n, se, sp, "c-p", conf)
  tp.wilson<-epiR::epi.prev(x, n, se, sp, "wilson", conf)
  tp.blaker<-epiR::epi.prev(x, n, se, sp, "blaker", conf)
  tp.sterne<-epiR::epi.prev(x, n, se, sp, "sterne", conf)
  tp.norm<- tp.normal(x, n, se, sp, conf)
  tp<- rbind(normal=tp.norm$tp, blaker=tp.blaker$tp, "c-p"=tp.cp$tp, "wilson"=tp.wilson$tp, sterne=tp.sterne$tp)
  ap<- rbind(normal=tp.norm$ap, blaker=tp.blaker$ap, "c-p"=tp.cp$ap, "wilson"=tp.wilson$ap, sterne=tp.sterne$ap)
  result<- list(normal=tp.norm, blaker=tp.blaker, "c-p"=tp.cp, wilson=tp.wilson, sterne=tp.sterne, all=list(ap=ap, tp=tp))
  return(result[[type]])
}
