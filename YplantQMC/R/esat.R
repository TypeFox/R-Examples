	
esat <- function(TdegC, Pa=101){	
	a <- 611.21
    b <- 17.502
    c <- 240.97
    f <- 1.0007 + 3.46 * 10^-8 * Pa * 1000
    esatval <- f * a * (exp(b * TdegC/(c + TdegC)))
	return(esatval)
}

