#' Power calculations for proportion tests (two-sided)
#'
#' @param n           Total number of observations
#' @param h           Effect size, Cohen's h
#' @param power       Power of test
#' @param sig.level   Significance level
#' @param type        Type of proportion tests, must be one of "one","two" (default), or "unequal". "one" means one sample proportion test. "two"/"unequal" means two sample (equal size/unequal size) proportion test.
#' @param ratio       The ratio of sample size 1 to sample size 2. Only will be used when 'type' is "unequal".
#' @export
#' @examples
#' ## one sample
#' power.proportions(n=600,h=0.3,type="one")
#' 
#' ## two sample with same sample size
#' power.proportions(h=0.2,n=600)
#' 
#' ## two sample with different sample size
#' power.proportions(h=0.2,n=1200,type="unequal",ratio=2)
#' 
power.proportions <- function (n = NULL, h = NULL, power = NULL, sig.level = 0.05,
                               type=c("two","one","unequal"),ratio=1){
    if (sum(sapply(list(n, h, power, sig.level), is.null)) != 1){
      stop("exactly one of 'n', 'h', 'power', and 'sig.level' must be NULL")
    }
    if (!is.null(n) && n < 2){
      stop("number of observations in each group must be at least 2")  
    }
    if (!is.null(power) && !is.numeric(power) || any(0>power | power>1)){
      stop("'power' must be numeric in [0, 1]")
    }
    if (!is.null(sig.level) && !is.numeric(sig.level) || any(0>sig.level | sig.level>1)){
      stop("'sig.level' must be numeric in [0, 1]")
    }
    type <- match.arg(type)
    
    if(type == "unequal" && any(ratio==1 | ratio <0)){
      stop("when 'type' is 'unequal', ratio must be not equal to 1 and great than 0")
    }
    
    nsample <- switch(type, one=1, two=2, unequal=(1+ratio)^2/ratio)
    h <- abs(h)
    
    p.body <- quote({
      pnorm(qnorm(sig.level/2, lower = FALSE) - h * sqrt(n/nsample), lower = FALSE) +
             pnorm(qnorm(sig.level/2, lower = TRUE) - h * sqrt(n/nsample), lower = TRUE)                      
    })
    
    if (is.null(power)){
      power <- eval(p.body)
    }else if (is.null(n)){
      n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+07))$root
    }else if (is.null(h)){
      h <- uniroot(function(h) eval(p.body) - power, c(1e-10, 5))$root     
    }else if (is.null(sig.level)){
      sig.level <- uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root             
    }
    
    if(type=="two"){
      NOTE <- "same sample sizes, n is the total sample size"
      METHOD <- "Difference of proportion power calculation for binomial distribution (arcsine transformation)"
      results <- structure(list(n = n, power = power, h = h, sig.level = sig.level,  
                                note=NOTE, method = METHOD), class = "power.htest")
    }else if(type=="one"){
      NOTE <- "n is the number of observation"
      METHOD <- "proportion power calculation for binomial distribution (arcsine transformation)"
      results <- structure(list(n = n, power = power, h = h, sig.level = sig.level,  
                                note=NOTE, method = METHOD), class = "power.htest")
    }else if(type=="unequal"){
      NOTE <- "different sample sizes, n is the total sample size"
      METHOD <- "difference of proportion power calculation for binomial distribution (arcsine transformation)"
      n1=n*ratio/(1+ratio)
      n2=n-n1
      results <- structure(list(n = n, n1=n1,n2=n2, power = power, h = h, sig.level = sig.level,  
                                note=NOTE, method = METHOD), class = "power.htest")
    }
    return(results)
}
