.sjdbcenv <- new.env(parent=emptyenv())

"sjdbcOptions"<-
function(...)
{
    # This function is used to get/set sjdbc options, 
    # similar to the S-PLUS options() function. Since
    # there is not frame 0 in R/TERR, it uses an
    # environment, .sjdbcenv, in the sjdbc package to
    #  store values durring a session.
    if(exists(".sjdbcOptions", envir = .sjdbcenv)) {
        current <- .sjdbcenv$.sjdbcOptions
    }
    else {
        current <- list(driverClass = "", con = "", user = "", password = "", keepAlive = FALSE, batchSize = 1000)
    }
    if(nargs() == 0) {
        return(current)
    }
    temp <- list(...)
    if(length(temp) == 1 && is.null(names(temp))) {
        arg <- temp[[1]]
        switch(mode(arg),
            list = temp <- arg,
            character = return(current[arg]),
            stop(paste("invalid argument:", arg)))
    }
    if(length(temp) == 0)
        return()
    n <- names(temp)
    if(is.null(n))
        stop("options must be given by name")
    changed <- current[n]
    names(changed) <- n
    current[n] <- temp
    assign(".sjdbcOptions", current, envir=.sjdbcenv)
    invisible(changed)
}
