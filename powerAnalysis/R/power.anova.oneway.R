#' Power calculations for balanced one-way analysis of variance tests
#'
#' @param groups          Number of groups
#' @param n               Number of observations (per group)
#' @param f               Effect size, Cohen's f
#' @param power           power of study
#' @param sig.level       significance level
#' @seealso               \code{\link{ES.anova.oneway}}
#' @export
#' @examples
#' power.anova.oneway(groups=4,n=20,f=0.28)
power.anova.oneway<-function (groups=NULL, n=NULL, f=NULL, power=NULL, sig.level=0.05)
{
  if (sum(sapply(list(groups, n, f, power, sig.level), is.null)) != 1){
    stop("exactly one of 'groups', 'n', 'f', 'power', and 'sig.level' must be NULL")
  } 
  if (!is.null(groups) && groups < 2){
    stop("number of groups must be at least 2")
  }
  if (!is.null(n) && n < 2){
    stop("number of observations in each group must be at least 2")
  }
  if (!is.null(f) && f < 0){
    stop("f must be positive")
  }
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 >sig.level | sig.level>1)){
    stop("'sig.level' must be numeric in [0, 1]")
  }

  p.body <- quote({
    lambda <- groups * n * f^2
    pf(qf(sig.level, groups - 1, (n - 1) * groups, lower.tail = FALSE), 
       groups - 1, (n - 1) * groups, lambda, lower.tail = FALSE)
  })
  
  if (is.null(power)){
    power <- eval(p.body)
  }else if (is.null(groups)){
    groups <- uniroot(function(groups) eval(p.body) - power,  c(2, 100))$root               
  }else if (is.null(n)){
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+08))$root
  }else if (is.null(f)){
    f <- uniroot(function(f) eval(p.body) - power, c(1e-07, 1e+07))$root                                              
  }else if (is.null(sig.level)){
    sig.level <- uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root                     
  }
  NOTE <- "n is number in each group"
  METHOD <- "Balanced one-way analysis of variance power calculation"
  structure(list(groups=groups, n=n, f=f, power=power, sig.level=sig.level, 
                 note=NOTE, method=METHOD), class="power.htest")
}
