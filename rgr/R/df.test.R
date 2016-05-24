df.test <-
function(dfname)
{
# Function to check for the existence of a data frame and whether it
# is already attached.  If the data frame is 'legitimate' the variable
# names are displayed with its dimensions.
#
    dfn <- deparse(substitute(dfname))
    if (exists(dfn)) {
        if (match(dfn, search(), nomatch = 0)){ 
            cat(" ", dfn, "is already attached\n")
        }
        else cat(" ", dfn, "is not attached\n")
        cat(paste("  Names of variables in ", dfn, ":\n ", sep = ""), 
            names(dfname), "\n")
        cat(paste("  Dimensions of ", dfn, ": ", sep = ""), dim(dfname), "\n")
    }
    else cat(" ", dfn, "is not in the workspace\n")
    invisible()
}
