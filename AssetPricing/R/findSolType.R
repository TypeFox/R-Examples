findSolType <- function(S,prices) {
#
# Examine the nature of the price sensitivity function and check
# for the existence of "prices", and accordingly return the resulting
# "solution type" ("cont", "disc", or "pwl").
    if(is.function(S)) {
        if(!is.null(prices)) {
                return("disc")
        } else {
            funtype <- attr(S,"funtype")
            if(!is.null(funtype) && funtype == "pwl") return("pwl")
            stop(paste("Argument \"S\" is a function but not of type \"pwl\"\n",
                       "and \"prices\" is not specified.  Something is wrong.\n"))
        }
    }
# Expression.
    if(is.expression(S)) {
        return("cont")
    }
# List.
    if(is.list(S)) {
# List of functions --- not allowed in the "pwl" setting.
        if(all(sapply(S,is.function))) {
	    if(!is.null(prices)) {
                return("disc")
            } else {
                stop(paste("All entries of \"S\" are functions but \"prices\"\n",
                           "is not specified.  Something is wrong.\n"))
            }
        }
# List of expressions.
        if(all(sapply(S,is.expression))) {
                return("cont")
        } else {
            stop(paste("At least one entry of \"S\" is neither a function\n",
                       "nor an expression.\n"))
        }
    }
    stop(paste("Argument \"S\" must be a function, an expression\n",
                       "or a list of such.\n"))
}
