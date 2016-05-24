##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snijders/siena
## *
## * File: effectsMethods.r
## *
## * Description: This file contains the print and edit methods for the
## * class sienaEffects
## *
## ****************************************************************************/
##@print.sienaEffects Methods
print.sienaEffects <- function(x, fileName=NULL, includeOnly=TRUE,
                               expandDummies=FALSE, ...)
{
    if (!inherits(x, "sienaEffects"))
        stop("not a legitimate Siena effects object")

    if (typeof(fileName)=="character")
    {
        sink(fileName, split=TRUE)
    }

    interactions <- x[x$shortName %in% c("unspInt", "behUnspInt") & x$include &
                            x$effect1 > 0, ]
    if (expandDummies)
    {
        if(includeOnly && !all(x[x$include, "timeDummy"] == ",")
                          || !all(x[, "timeDummy"] == "," ))
        {
            x <- sienaTimeFix(x)$effects
            x <- fixUpEffectNames(x)
        }
        else
        {
            if (nrow(interactions) > 0)
            {
                x <- fixUpEffectNames(x)
            }
        }
    }
    if (nrow(x) > 0)
    {
        nDependents <- length(unique(x$name))
        userSpecifieds <- x$shortName[x$include] %in% c("unspInt", "behUnspInt")
        endowments <- !x$type[x$include] %in% c("rate", "eval")
        timeDummies <- !x$timeDummy[x$include] == ","
        specs <- x[, c("name", "effectName", "include", "fix", "test",
                       "initialValue", "parm")]
        if (includeOnly)
        {
            specs <- specs[x$include, ]
        }
        if (nDependents == 1)
        {
            specs <- specs[, -1]
        }
        if (any(endowments))
        {
            specs <- cbind(specs, type=x[x$include, "type"])
        }
        if (any(timeDummies))
        {
            specs <- cbind(specs, timeDummy=x[x$include, "timeDummy"])
        }
        if (any(userSpecifieds))
        {
            specs <- cbind(specs, x[x$include, c("effect1", "effect2")])
            if (any (x$effect3[x$include] > 0))
            {
                specs <- cbind(specs, effect3=x[x$include, "effect3"])
            }
        }
        specs[, "initialValue"] <- format(round(specs$initialValue,digits=5),
                                          width=10)
        if (nrow(specs) > 0)
        {
            row.names(specs) <- 1:nrow(specs)
            print(as.matrix(specs), quote=FALSE)
        }
        else
        {
            print.data.frame(specs)
        }
    }
    else
    {
        print.data.frame(x)
    }
    if (typeof(fileName)=="character")
    {
        sink()
    }

    invisible(x)
}

##@summary.sienaEffects Methods
summary.sienaEffects <- function(object, fileName=NULL, includeOnly=TRUE,
                                 expandDummies=FALSE, ...)
{
    if (!inherits(object, "sienaEffects"))
        stop("not a legitimate Siena effects object")
    if (expandDummies && (includeOnly && !all(object[object$include,
                                                     "timeDummy"] == ",")
        || !all(object[, "timeDummy"] == ",")))
    {
        object <- sienaTimeFix(object)$effects
    }
    if (includeOnly)
    {
        object <- object[object$include, ]
    }
    class(object) <- c("summary.sienaEffects", class(object))
    object
}

##@print.summary.sienaEffects Methods
print.summary.sienaEffects <- function(x, fileName=NULL, ...)
{
    if (!inherits(x, "summary.sienaEffects"))
        stop("not a legitimate summary of a Siena effects object")
    ## find out if any columns need removing because they will not print

    problem <- sapply(x, function(x)
                  {
                      inherits(try(unlist(x[[1]]), silent=TRUE), "try-error")
                  }
                      )
    if (any(problem))
    {
        x1 <- x[, problem]
        x <- x[, !problem]
    }
    if (typeof(fileName)=="character")
    {
        sink(fileName, split=TRUE)
    }
    print.data.frame(x)
    if (any(problem))
    {
        lapply(x1, print)
    }
    if (typeof(fileName)=="character")
    {
        sink()
    }
    invisible(x)
}

edit.sienaEffects <- function(name, ...)
{
	if (!interactive())
	{
		return(name)
	}
    ## store the original column order
    originalNames <- names(name)
    ## move function name and other things out of the way
    priorityColumns <- c("name", "effectName", "type", "include", "fix",
                         "test", "initialValue", "parm", "shortName")
    priorityX <- name[, priorityColumns]
    notPriorityX <- name[, -c(match(priorityColumns, names(name)))]
    name <- cbind(priorityX, notPriorityX)

    ## get edit.data.frame to do the actual edit
    tmp <- NextMethod(, , edit.row.names=FALSE)

    ## re-sort the columns
    tmp <- tmp[, match(originalNames, names(tmp))]
    class(tmp) <- c("sienaEffects", class(tmp))
    tmp
}
