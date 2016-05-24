#===============================================================================
# Name           : cv
# Author(s)      : José Cláudio Faria/UESC/DCET
#                  Walmes Marques Zeviani UFPR/DE
# Date           : 2014/08/16 - 09:24:11
# Version        : v7
# Aim            : Calculate coefficient variation from lm, aov and
#                  aovlist objects
#===============================================================================

# Arguments:
# x    : an lm, aov or aovlist object
# round: rounds the values of cv to the specified number of decimal places (default 2)

cv <- function(x, 
               round=2)
{
  if(is.null(x))
    stop('An object of class lm, aov or aovlist must be informed!')

  if(inherits(x,
              'lm')) {
    qmee <- with(x,
                 sum(residuals^2) / df.residual)
    m <- mean(x$fitted.values)
    cv   <-  100 * sqrt(qmee) / m
    names(cv) <- 'cv'
    return(round(cv,
                 round))
  }
  if(inherits(x,
              'aovlist')) {
    sm.x <- summary(x)
    qmee <- sapply(sm.x,
                   function(i){
                     i[[1]]["Residuals", 3]
                   })
    qmee <- qmee[!is.na(qmee)]
    m <- x$'(Intercept)'$coefficients[[1]]
    cv <- 100 * sqrt(qmee) / m
    names(cv) <- paste('cv(',
                       letters[1:length(qmee)],
                       ')',
                       sep='')
    return(round(cv,
                 round))
  }
}
