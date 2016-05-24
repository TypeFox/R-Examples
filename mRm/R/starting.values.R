#    Function: Provides starting values for the EM algorithm (Rost 1991).
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

starting.values <- function(resp.pat, cl){
k <- dim(resp.pat)[2] -2
N <- sum(resp.pat[, k+1])

#generates matrix of Sum = 0 standardized item parameters

item.para <- matrix(runif(k*cl, -2, 2), nrow = k, ncol = cl)
fun.stand.0 <- function(x){x - mean(x)}
item.para <- apply(item.para, 2, fun.stand.0)

#class size parameters with restriction sum(pi.c) = 1 - pi.0 - pi.k

pi.0 <- sum(resp.pat[, k+1][resp.pat[, k+2] == 0])/N
pi.k <- sum(resp.pat[, k+1][resp.pat[, k+2] == k])/N

vec <- runif(cl, 0, 1)
pi.c <- vec*(1-pi.0 - pi.k)/(sum(vec))

#latent score probabilities with restriction sum_r=1^(k-1) pi.r.c = 1 for all cl \in {1,...,C}

cond.pat.freq <- matrix(runif(cl*(k-1), 0, 1), ncol = cl, nrow = k-1)
fun.stand.1 <- function(x){x/sum(x)}
cond.pat.freq <- apply(cond.pat.freq, 2, fun.stand.1)
return(list("beta" = item.para, "pi.c" = pi.c, "cond.pat.freq" = cond.pat.freq))}

