#############################################################
#                                                           #
#	wle.smooth function                                 #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: October, 10, 2000                             #
#	Version: 0.3                                        #
#                                                           #
#	Copyright (C) 2000 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.smooth <- function(weight=0.31,costant=3,level=0.2,dimension=1,raf="HD",interval=c(0.00001,0.5),tol=10^-6,max.iter=1000) {

raf <- switch(raf,
	HD = 1,
	NED = 2,
	SCHI2 = 3,
	-1)

if (raf==-1) stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

delta <- function(smooth,costant,level,dimension){
level*(((smooth+1)/smooth)^(dimension/2)*exp(costant^2/(2*(dimension+1)))-1)}

if (raf==3) {w <- function(smooth,costant,level,dimension,weight){(1-((delta(smooth,costant,level,dimension)**2)/((delta(smooth,costant,level,dimension)**2) + 2)))-weight}
} else {
if (raf==2) {
adelta <- function(d) {2-(2+d)*exp(-d)} 
} else {
adelta <- function(d) {2*(sqrt(d+1)-1)}
}
w <- function(smooth,costant,level,dimension,weight){
(adelta(delta(smooth,costant,level,dimension))+1)/(delta(smooth,costant,level,dimension)+1)-weight
}
}

result <- uniroot(w,interval=interval,costant=costant,level=level,dimension=dimension,weight=weight,maxiter=max.iter,tol=tol)

result$call <- match.call()

class(result) <- "wle.smooth"

return(result)
}

#############################################################
#                                                           #
#	print.wle.smooth function                           #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: October, 10, 2000                             #
#	Version: 0.3                                        #
#                                                           #
#	Copyright (C) 2000 Claudio Agostinelli              #
#                                                           #
#############################################################

print.wle.smooth <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("\nBandwidth: ",format(x$root, digits=digits))
    cat("\n")
    invisible(x)
}



