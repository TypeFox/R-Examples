#    Function: Simulates data matrices conforming to mixed Rasch models.
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

sim.mrm <- function(N.sample, N.items, cl.prob, item.para = NULL, pers.para = NULL, seed = NULL){

if (!is.null(seed)) set.seed(seed)
N.classes <- length(cl.prob)
memb <- rmultinom(1, N.sample, cl.prob)

if (is.null(item.para)){
	item.para <- apply(matrix(seq(-2, 2, length = N.items), byrow = T, ncol = N.items, nrow = N.classes), 1, sample)
	}

Q <- apply(item.para, 1, rep.int, times = memb)

if (is.null(pers.para)){
	pers.para <- rnorm(N.sample)
}

pers.para <- matrix(pers.para, byrow = F, nrow = N.sample, ncol = N.items)

probs <- exp(pers.para + Q)/(1 + exp(pers.para + Q))
M <- matrix(runif(N.items * N.sample), ncol = N.items, nrow = N.sample)

data.matrix <- (probs>=M)*1

return(list("data.matrix" = data.matrix, "beta" = item.para, "emp.probs" = memb/N.sample, "xi" = pers.para[,1]))
}
