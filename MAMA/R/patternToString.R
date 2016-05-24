########################################################################
# find pattern names for 0-1-2 matrix as X.discret
##########
patternToString <- function(X.discret)
{
 library(xtable)

    W <- ncol(X.discret)
#    words <- apply(X.discret,2,prettyNum)
    words <- apply(X.discret,2,as.character)
    l <- apply(words,1,paste,collapse="")
    names(l) <- rownames(X.discret)
    return(l)
   }

