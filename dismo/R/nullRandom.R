# Author: Robert J. Hijmans
# Date : April 2010
# Version 1
# Licence GPL v3


nullRandom <- function(x, model, n=25, rep=25, pa=FALSE) {
	e <- list()
	stopifnot( n < nrow(x) )
	for (r in 1:rep) {
		#select n presence records
		i <- sample(nrow(x), n)
		pres <- x[i, ]
		absc <- x[-i, ]
		if (pa) {
			d <- data.frame(rbind(pres, absc))
			v <- c(rep(1, nrow(pres)), rep(0, nrow(absc)))
			m <- model(d, v)
		} else {
			m <- model(pres)
		}
		e[[r]] <- evaluate(pres, absc, m)
		message('-', appendLF = FALSE)
		if (r%%50 == 0) message(" ", r)
		flush.console()
	}
	if (r%%50 != 0) { message(" ", r) } else { message("") }
	e
}


.nullRandom2 <- function(p, a, model, n=25, rep=25, pa=FALSE) {
	e <- list()
	stopifnot(n < nrow(p))
	for (r in 1:rep) {
		#select n presence records
		i <- sample(nrow(p), n)
		trainpres <- p[i, ]
		testpres <- p[-i, ]
		j <- sample(nrow(a), nrow(testpres) * 2)
		testabs <- a[j, ]
		trainabs <- a[-j, ]

		if (pa) {
			x <- rbind(trainpres, trainabs)
			v <- c(rep(1, nrow(trainpres)), rep(0, nrow(trainabs)))
			m <- model(x, v)
		} else {
			m <- model(trainpres)
		}
		e[[r]] <- evaluate(testpres, testabs, m)
		message('-', appendLF=FALSE)
		if (r%%50 == 0) cat(" ", r)
		flush.console()
	}
	if (r%%50 != 0) { message(" ", r) } else { message("") }
	e
}

