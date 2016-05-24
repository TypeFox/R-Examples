##
## INPUT:
## method: character string describing optimization method used
##
## RETURN:
## vector of integer values corresponding to convergence codes indicating
## success (which differ by method in maxLik)
## 
convergenceCriterion <- function(method)
{
    switch(tolower(method),
           `newton-raphson` = c(1L, 2L),
           nr = c(1L, 2L),
           bfgs = 0L,
           bfgsr = c(1L, 2L),
           `bfgs-r` = c(1L, 2L),
           bhhh = c(1L, 2L),
           `nelder-mead` = 0L,
           nm = 0L,
           sann = 0L)
}
