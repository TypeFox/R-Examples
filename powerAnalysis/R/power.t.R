#' Power calculations for t-test
#'
#' @param es              effect size. 
#' @param n               total number of observations/pairs
#' @param power           power of study
#' @param sig.level       significance level
#' @param ratio           the ratio of sample size 1 to sample size 2. Only will be used when 'type' is "unequal".
#' @param type            type of t test, must be one of "one","two" (default), "paired", or "unequal". "one" means one sample t test, which test whether the population mean is equal to a specified value. "two"/"unequal" means two sample (equal size/unequal size) t test, which is used to ascertain how likely an observed mean difference between two groups would be to occur by chance alone. "paired" means paired t-test (also called the correlated t-test and the t-test for dependent means), which is used to ascertain how likely the difference between two means that contain the same (or matched) observations is to occur by chance alone.
#' @param alternative     One- or two-sided test, must be one of "two.sided" (default), "left", "right"
#' @seealso               \code{\link{ES.t.one}}
#' @seealso               \code{\link{ES.t.two}}
#' @seealso               \code{\link{ES.t.paired}}
#' @export
#' @examples
#' ## one sample two sided test, calculate power
#' power.t(es=0.2,n=60,sig.level=0.10,type="one",alternative="two.sided")
#' 
#' ## one sample one sided (left tail) test, calculate power
#' power.t(es=0.2,n=60,sig.level=0.10,type="one",alternative="left")
#' 
#' ## one sample one sided (right tail) test, calculate power
#' power.t(es=0.2,n=60,sig.level=0.10,type="one",alternative="right")
#' 
#' ## one sample two sided test, calculate sampe size
#' power.t(es=0.2,power=0.8,sig.level=0.05,type="one",alternative="two.sided")
#' 
#' ## one sample two sided test, calculate effect size
#' power.t(n=200,power=0.8,sig.level=0.05,type="one",alternative="two.sided")
#' 
#' ## one sample two sided test, calculate sig.level
#' power.t(es=0.2,n=200,power=0.8,type="one",alternative="two.sided")
#' 
#' ## paired sample two sided test, calculate power
#' power.t(es=0.559,n=40,sig.level=0.05,type="paired",alternative="two.sided")
#' 
#' ## paired sample two sided test, calculate sample size
#' power.t(es=0.15,power=0.8,sig.level=0.05,type="paired",alternative="two.sided")
#' 
#' ## paired sample two sided test, calculate effect size
#' power.t(n=200,power=0.8,sig.level=0.05,type="paired",alternative="two.sided")
#' 
#' ## two sample two sided test, calculate power
#' power.t(es=0.15,n=300,sig.level=0.05,type="two",alternative="two.sided")
#' 
#' ## two sample two sided test, calculate sample size
#' power.t(es=0.15,power=0.8,sig.level=0.05,type="two",alternative="two.sided")
#' 
#' ## two sample two sided test, calculate effect size
#' power.t(n=300,power=0.8,sig.level=0.05,type="two",alternative="two.sided")
#' 
#' ## two sample (unequal size), calculate sample size
#' power.t(es=0.15,power=0.8,sig.level=0.05,type="unequal",ratio=2,alternative="two.sided")
#' 
## two sample (unequal size), calculate power
#' power.t(es=0.1,n=3000,sig.level=0.05,type="unequal",ratio=2,alternative="two.sided")
power.t <- function (es=NULL, n=NULL, power = NULL, sig.level = NULL, ratio=1, 
                     type=c("two","paired","one","unequal"), 
                     alternative = c("two.sided", "left","right")){
  if (sum(sapply(list(es, n, power, sig.level), is.null)) != 1){
    stop("exactly one of 'es', 'n', 'power', and 'sig.level' must be NULL")
  }
  
  if (!is.null(n) && is.numeric(n) && n < 1){
    stop("'n' must be numeric larger than 1")
  }
  
  if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 1)){
    stop("'power' must be numeric in [0, 1]")
  }
    
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > sig.level | sig.level > 1)){
    stop("'sig.level' must be numeric in [0, 1]")
  }
   
  alternative <- match.arg(alternative)
  ttside<-switch(alternative, left = 1, two.sided = 2, right=3)
  tside <- switch(alternative, left = 1, two.sided = 2, right =1)
  
  type <- match.arg(type)
  tsample <- switch(type, one = 1, two = 2, unequal=2, paired = 1)
  
  if(tside == 2 && !is.null(es)){
    es <- abs(es)
  }
  n1=NULL
  n2=NULL
  if(type!="unequal"){
    if(ttside==1){
      p.body <- quote({
        nu <- (n - 1) * tsample
        pt(qt(sig.level/tside, nu, lower.tail = TRUE), nu, ncp = sqrt(n/tsample) * es, 
           lower.tail = TRUE)
      })
    }else if(ttside==2){
      p.body <- quote({
        nu <- (n - 1) * tsample
        pt(qt(sig.level/tside, nu, lower.tail = FALSE), nu, ncp = sqrt(n/tsample) * es, 
           lower.tail = FALSE) + pt(-qt(sig.level/tside, nu, lower.tail = FALSE), nu, 
                                    ncp = sqrt(n/tsample) * es, lower.tail = TRUE)
      })
    }else if(ttside==3){
      p.body <- quote({
        nu <- (n - 1) * tsample
        pt(qt(sig.level/tside, nu, lower.tail = FALSE), nu, ncp = sqrt(n/tsample) * es, 
           lower.tail = FALSE) 
      })
    }
  }else{
    if(ratio==1 || ratio <= 0){
      stop("when 'type' is 'unequal', 'ratio' must be greater than 0 and does not equal to 1")
    }
    
    if(ttside==1){
      p.body <- quote({
        nu <- n-2
        pt(qt(sig.level/tside, nu, lower = TRUE), nu, 
           ncp = es*(1/sqrt((1+ratio)/ratio/n + (1+ratio)/n)),lower = TRUE)
      })
    }else if(ttside==2){
      p.body <- quote({
        nu <- n-2
        qu <- qt(sig.level/tside, nu, lower = FALSE)
        pt(qu, nu, ncp = es*(1/sqrt((1+ratio)/ratio/n + (1+ratio)/n)), 
           lower = FALSE) + pt(-qu, nu,ncp = es*(1/sqrt((1+ratio)/ratio/n + (1+ratio)/n)), 
                               lower = TRUE)
        })
    }else if(ttside==3){
      p.body <- quote({
        nu <- n-2
        pt(qt(sig.level/tside, nu, lower = FALSE), nu, 
           ncp = es*(1/sqrt((1+ratio)/ratio/n + (1+ratio)/n)), lower = FALSE)
      })
    }
  }

  
  if(is.null(es)){
    es <- uniroot(function(es) eval(p.body) - power, c(1e-08, 20))$root                                                                                                                                 
  }else if(is.null(n)){ 
    n <- uniroot(function(n) eval(p.body) - power, c(2+1e-02, 1e+10))$root
  }else if(is.null(power)){
    power <- eval(p.body)
  }else if(is.null(sig.level)){
    sig.level <- uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root                    
  }
  
  if(type == "unequal"){
    n1=ratio*n/(1+ratio)
    n2=n/(1+ratio)
  }
  NOTE <- "n is number in *each* group"
  METHOD <- "Independent two-sample t-test power calculation"
  
  NOTE <- switch(type, paired = "n is number of *pairs*", two = "n is number in *each* group", 
                 one="n is the number of observations", unequal="n is the total number of observations")
  METHOD <- paste(switch(type, one = "One-sample", two = "Two-sample (equal size)", 
                         paired = "Paired", unequal="Two-sample (unequal size)"), 
                  "t test power calculation")
  results=NULL
  if(type!="unequal"){
    results=structure(list(es = es, n=n, power = power, sig.level = sig.level, 
                           alternative = alternative, note = NOTE, method = METHOD), 
                           class = "power.htest")              
  }else{
    results=structure(list(es = es, power = power, n=n, n1=n1, n2=n2, ratio=ratio, 
                           sig.level = sig.level, alternative = alternative,
                            note = NOTE, method = METHOD), class = "power.htest") 
  }
  return(results)
}
