### $Id: aaa.R 151 2011-01-02 15:26:27Z bhm $

###
### Start-up setup
###
#.onLoad <- function(...){
    ## To avoid warnings in R CMD check:
    #if(!exists(".baseline.current"))
    #    assign(".baseline.current", list(), .GlobalEnv)
#}

.baselineEnv <- new.env(parent=emptyenv())
baselineEnv <- function() .baselineEnv
putBaselineEnv <- function(x, value) assign(x, value, envir=baselineEnv())
getBaselineEnv <- function(x, mode="any") get(x, envir=baselineEnv(), mode=mode, inherits=FALSE)

putBaselineEnv("baseline.result", list())
putBaselineEnv("baseline.current", list())
