# non-exported utility functions

# fr_checkstart
# Provides packageswide checking of start values
fr_checkstart <- function(start, start_name){
    if(length(start)==0){
        stop(paste0("You didn't provide starting values.\n",
                    "   It's impossible to fit anything without knowing what to optimise!"))
    }
    if(!is.list(start) | is.null(names(start))){
        stop(paste0(start_name, " must be a list containing single, named numeric values."))
    }
    if(any(lapply(start, length)!=1)){
        stop(paste0("The items in ", start_name, " must be single, named numeric values."))
    }
    if(!(all(is.numeric(unlist(start))))){
        stop(paste0("The items in ", start_name, " must be single, named numeric values."))
    }
}


# fr_setupout
# Utility function to help clean up the code in the various statistic functions
# NB: This setups a single vector, as used by the various dual-duty functional response functions (e.g. rogersII_fit) NOT the output of frair_fit
# Returns a vector of (correctly) named NAs, based on the start and fixed input to the statitic function
fr_setupout <- function(start, fixed, samp){
    out <- c(rep(NA, times=length(names(start))*2), rep(NA, times=length(names(fixed))), samp)
    outnames <- NULL
    # Setup optimised variable output
    if(!is.null(start)){
        for (i in 1:length(names(start))){
            outnames <- c(outnames, names(start)[i], paste(names(start)[i], 'var', sep=''))
        }
    }
    # Setup fixed variable output
    if(!is.null(fixed)){
        for (i in 1:length(names(fixed))){
            outnames <- c(outnames, names(fixed)[i])
        }
    }
    names(out) <- c(outnames, rep('', times=length(samp)))
    return(out)
}

# fr_setpara parallel
# Utility to setup parallel processing in Windows
fr_setpara <- function(boot, windows){
        if(boot && windows){
        is_load <- require(frair, warn.conflicts=FALSE, quietly=TRUE)
        if(is_load==FALSE){
            stop('Error establishing workspace for parallel computing in Windows.')
        }
    }
}

# From http://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
fr_catchlist <- function(expr) {
    val <- NULL
    myWarnings <- NULL
    wHandler <- function(w) {
        myWarnings <<- c(myWarnings, paste0('Warning: ', w$message))
        invokeRestart("muffleWarning")
    }
    myError <- NULL
    eHandler <- function(e) {
        myError <<- paste0('Error: ', e$message)
        NULL
    }
    val <- tryCatch(withCallingHandlers(expr, warning = wHandler), error = eHandler)
    list(value = val, warnings = myWarnings, error=myError)
}

## The startup method
# .onAttach <- function(lib, pkg)  {
#     packageStartupMessage('This is FRAIR (v. ', utils::packageDescription("frair", field="Version"), ')', appendLF = TRUE)
# }