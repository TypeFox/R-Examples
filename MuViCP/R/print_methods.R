print.bpa <- function(x, verbose = FALSE, ...)
    {
        cat('A bpa object.', '\n')
        print(x$get.m())
        if(verbose == TRUE)
            {
                cat('The following functions are available:\n')
                print(names(x))
            }
    }

print.bpamat <- function(x, ...)
    {
        cat('A bpa mat object with information on',
            length(x$get.setlist()) - 1  , 'classes and',
            length(x$get.pointlist()) - 1, 'points', '\n')
        cat('The following functions are available:\n')
        print(names(x))        
    }
