#    Function: mrm main routine for fitting mixed Rasch models
#    Copyright (C) 2011  David Preinerstorfer
#    david.preinerstorfer@univie.ac.at
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details. A copy may be obtained at
#    http://www.r-project.org/Licenses/

mrm <-function(data.matrix, cl, max.it = 1000, conv.crit = .0001){

#check and prepare data 

if(is.data.frame(data.matrix) == F & is.matrix(data.matrix) == F){stop("Inappropriate data object!")}
data.matrix <- na.omit(as.matrix(data.matrix))	
csum <- apply(data.matrix, 2, sum)
keep.items <- (csum < dim(data.matrix)[1] & csum > 0)
which.removed <- which(keep.items == FALSE)
if(sum((keep.items)==FALSE)>0){cat("Removed items: ", which.removed, "\n"); warning("Uninformative Item(s) have been removed!")}

data.matrix <- data.matrix[,keep.items]
dim.data <- dim(data.matrix)
N.whole <- dim.data[1]	
k <- dim.data[2]

#get response patterns 
resp.pat <- resp.patterns(data.matrix)
resp.pat <- resp.pat[ifelse(resp.pat[, k + 2] > 0 & resp.pat[, k + 2] < k, TRUE, FALSE), ]

N <- dim(resp.pat)[1]

#get starting values for EM algorithm; 
starting.vals <- starting.values(resp.pat, cl)

beta <- starting.vals$beta
prelim.beta <- beta * 0
pi.r.c <- starting.vals$cond.pat.freq
class.size <- starting.vals$pi.c

#call "em.cpp"

fit <-	.C("em", as.integer(k), as.integer(N.whole), as.integer(dim(resp.pat)[1]), as.integer(k-1), 
	         as.integer(cl), as.integer(max.it), as.double(conv.crit), 
		 as.double(as.vector(beta)), as.double(as.vector(pi.r.c)), as.double(as.vector(resp.pat)), as.double(as.vector(class.size)), 
		 as.double(as.vector(rep(0, 3))), as.integer(0), as.integer(0), as.integer(0), PACKAGE = "mRm")

B <- matrix(as.vector(unlist(fit[[8]])), nrow = k, ncol = cl)
P <- matrix(as.vector(unlist(fit[[9]])), nrow = dim(pi.r.c)[1], ncol = cl)
C <- matrix(as.vector(unlist(fit[[11]])), nrow = 1, ncol = cl)
L <- as.vector(unlist(fit[[12]]))
I <- unlist(fit[[13]])
n.para <- unlist(fit[[14]])
boundary <- unlist(fit[[15]])

obj <- list("beta" = B, "pi.r.c" = P, "class.size" = C, "logLik" = L[1], "AIC" = L[2], "BIC" = L[3], "number.of.iterations" = I,
"number.of.parameters" = n.para, "conv.to.bound" = boundary)

class(obj) <- "mrm"

return(obj)

return(); 
}

