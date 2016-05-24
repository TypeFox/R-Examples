#############################################################
#                                                           #
#	mme.gamma function                                 #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: May 18, 2007                                  #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2007 Claudio Agostinelli              #
#                                                           #
#############################################################

mme.gamma <- function(x) {

x <- as.vector(x)
result <- list()
  
#### central moments
m1 <- function(x) mean(x)
m2 <- function(x, m1) c(sum((x-m1)^2)/length(x))
m3 <- function(x, m1) c(sum((x-m1)^3)/length(x))

#### wmme solutions
alphap <- function(m2, m3) 4*m2^3/m3^2          #shape
betap <- function(m2, m3) m3/(2*m2)             #scale
gammap <- function(m1, m2, m3) (m1 - 2*m2^2/m3) #location

    m1x <- m1(x)
    m2x <- m2(x, m1=m1x)
    m3x <- m3(x, m1=m1x)
    aa <- alphap(m2x, m3x)
    bb <- betap(m2x, m3x)
    cc <- gammap(m1x, m2x, m3x)
   result$scale <- c(bb)
   result$rate <- c(1/bb)  
   result$shape <- c(aa)
   result$location <- c(cc)
   result$call <- match.call()
   class(result) <- "mme.gamma"
   return(result)
}

#############################################################
#                                                           #
#	print.mme.gamma function                            #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: May, 18, 2007                                 #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2007 Claudio Agostinelli              #
#                                                           #
#############################################################

print.mme.gamma <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("Scale:\n")
    print.default(format(x$scale, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("Rate:\n")
    print.default(format(x$rate, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")  
    cat("Shape:\n")
    print.default(format(x$shape, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")  
    cat("Location:\n")
    print.default(format(x$location, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    invisible(x)
}






