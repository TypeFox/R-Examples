"mklevel" <-

function(x, y, z, p = c(0.10000000000000001, 0.25, 0.5, 0.75, 



	0.90000000000000002))



{



	p <- 1. - p



	n <- length(x)



	qx <- c(x[2.] - x[1.], x[3.:n] - x[2.:(n - 1.)], x[n] - x[



		n - 1.])/2.



	n <- length(y)



	qy <- c(y[2.] - y[1.], y[3.:n] - y[2.:(n - 1.)], y[n] - y[



		n - 1.])/2.



	qxy <- outer(qx, qy)



	it <- order(z)



	temp <- z[it] * qxy[it]



	temp <- cumsum(temp)/sum(temp)



	z <- z[it]



	levels <- rep(NA, length(p))



	for(k in 1.:length(p)) {



		levels[k] <- min(z[(temp > p[k])])



	}



	list(p = 1. - p, levels = levels)



}

