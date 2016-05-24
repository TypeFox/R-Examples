#############################################################
#                                                           #
#	mle.inversegaussian function                        #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: April, 18, 2011                               #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2011 Claudio Agostinelli              #
#                                                           #
#############################################################

mle.inversegaussian <- function(x) {
  x <- as.vector(x)
  result <- list()
  mu <- mean(x)
  lambda <- mu^2*length(x)/sum(((x-mu)^2/x))
  result$mu <- mu
  result$lambda <- lambda
  result$call <- match.call()
  class(result) <- "mle.inversegaussian"
  return(result)
}

#############################################################
#                                                           #
#	print.mle.inversegaussian function                  #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: April, 18, 2011                               #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2011 Claudio Agostinelli              #
#                                                           #
#############################################################

print.mle.inversegaussian <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("Location:\n")
    print.default(format(x$mu, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("Precision:\n")
    print.default(format(x$lambda, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    invisible(x)
}



