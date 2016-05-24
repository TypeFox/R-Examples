
.Normal <- function()
    list(p = pnorm, d = dnorm, q = qnorm, 
         ### see also MiscTools::ddnorm
         dd = function(x) -dnorm(x = x) * x,
         ddd = function(x) dnorm(x = x) * (x^2 - 1), 
         name = "normal")

.Logistic <- function()
    list(p = plogis, d = dlogis, q = qlogis,
         dd = function(x) 
             (exp(x) - exp(2*x)) / (1 + exp(x))^3,
         ddd = function(x) 
             (exp(x) - 4*(exp(2*x)) + exp(3*x)) / (1 + exp(x))^4, 
         name = "logistic")
#         (2 * exp(-x)^2 / (1 + exp(-x))^3) - 
#          (exp(-x) / (1 + exp(-x))^2))

.MinExtrVal <- function()
    list(p = function(x) 1 - exp(-exp(x)),
         q = function(p) log(-log(1 - p)),
         d = function(x, log = FALSE) {
             ret <- x - exp(x)
             if (!log) return(exp(ret))
             ret
         },
         dd = function(x)
             (exp(x) - exp(2*x)) / exp(exp(x)),
         ddd = function(x)
             (exp(x) - 3*exp(2*x) + exp(3*x)) / exp(exp(x)),
         name = "minimum extreme value")

.distr <- function(which = c("Normal", "Logistic", 
                             "MinExtrVal")) {
    which <- match.arg(which)
    do.call(paste(".", which, sep = ""), list())
}
