#' Linear model p-value
#' 
#' Extract F-test p-value from a linear model object. Can also use \code{broom::glance(fit)}. Originally described at \url{http://www.gettinggeneticsdone.com/2011/01/rstats-function-for-extracting-f-test-p.html}.
#'
#' @param modelobject A model object of class \code{lm}. 
#'
#' @return The p-value on the f-test of a linear model object testing the null hypothesis that R^2==0.
#' 
#' @importFrom stats pf
#'
#' @examples
#' # simulate some (e.g. SNP genotype) data
#' set.seed(42)
#' n=20
#' d=data.frame(x1=rbinom(n,2,.5), x2=rbinom(n,2,.5))
#' d=transform(d, y=x1+x2+rnorm(n))
#' #fit the linear model
#' fit=lm(y ~ x1 + x2, data=d)
#' #shows that the F-test is 0.006641
#' summary(fit)
#' #can't access that p-value using this
#' names(summary(fit)) 
#' # this doesn't work either
#' names(fit)
#  # use the lmp() function:
#' lmp(fit)
#' 
#' @export
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm'.")
    f <- summary(modelobject)$fstatistic
    # f[1]=value, f[2]=numeratorDF, f[3]=denominatorDF
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}
