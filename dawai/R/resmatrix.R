resmatrix <-
function(text, numgroups, dimension)
{
    if(class(text) != "character" || nchar(text) <= 2 || !(tolower(substr(text, 1, 1))
       %in% c("s", "t")) || !(substr(text, 2, 2) %in% c(">", "<")))
    {
        cat("Invalid constraints text.\n\n")
        return(NULL)
    }
    variables <- unlist(strsplit(substr(text, 3, nchar(text)), ","))
    if(suppressWarnings(sum(is.na(as.numeric(variables))) > 0))
    {
        cat("Variables must be denoted by numbers varying from 1 to dimension.\n\n")
        return(NULL)
    }
    variables <- as.numeric(variables)
    if(sum(variables > dimension) > 0)
    {
        cat("Variables must be denoted by numbers varying from 1 to dimension.\n\n")
        return(NULL)
    }
    ll <- length(variables)
    conmatrix <- array(0, c({numgroups - 1} * ll, numgroups * dimension))
    if(tolower(substr(text, 1, 1)) == "s")
    {
        for(group in 2:numgroups)
        {
            conmatrix[cbind(1:ll + {group - 2} * ll, variables + {group - 2} * dimension)] <- -1
            conmatrix[cbind(1:ll + {group - 2} * ll, variables + {group - 1} * dimension)] <- 1
        }
    }else
        for(group in 2:numgroups)
        {
            conmatrix[cbind(1:ll + {group - 2} * ll, variables)] <- -1
            conmatrix[cbind(1:ll + {group - 2} * ll, variables + {group - 1} * dimension)] <- 1
        }
    if(substr(text, 2, 2) == "<")
        conmatrix <- {-1} * conmatrix
    return(conmatrix)
}
