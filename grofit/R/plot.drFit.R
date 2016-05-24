plot.drFit <-
function(x, ...)
{
# x an object of class drFit

n <- length(x$drFittedSplines)

# /// plot all drFitSpline objects
for (i in 1:n){
#x11()
try(plot(x$drFittedSplines[[i]]))
title(as.character(x$drFittedSplines[[i]]$drID))
}

}

