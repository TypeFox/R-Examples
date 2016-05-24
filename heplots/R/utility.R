# critical value of Roy's largest root test

#lambda.crit <- function(alpha, p, dfh, dfe){
#    d <- max(p, dfh)
#    nu <- dfe - d + dfh
#    (d/nu) * qf(alpha, d, nu, lower.tail=FALSE)
#    }

lambda.crit <- function(alpha, p, dfh, dfe, test.statistic=c("Roy", "HLT", "Hotelling-Lawley")){
	test.statistic <- match.arg(test.statistic)
	switch(test.statistic,
		Roy = Roy.crit(alpha, p, dfh, dfe),
		HLT = HLT.crit(alpha, p, dfh, dfe),
		"Hotelling-Lawley" = HLT.crit(alpha, p, dfh, dfe)
		)
}
# see: http://wiki.math.yorku.ca/index.php/Statistics:_Ellipses
## Critical value for \lambda_1 in Roy test
Roy.crit <- function(alpha, p, dfh, dfe){
    df1 <- max(p, dfh)
    df2 <- dfe - df1 + dfh
    (df1/df2) * qf(alpha, df1, df2, lower.tail=FALSE)
}

## Critical value for \bar{\lambda_i} in HLT test
HLT.crit <- function ( alpha, p, dfh, dfe) {
	s <- min(p, dfh)
	m <- (abs(p-dfh)-1)/2
	n <- (dfe-p-1)/2
	df1 <- 2*m + s + 1
	df2 <- 2*(s*n +1)
  s * (df1/df2) * qf(alpha, df1, df2, lower.tail=FALSE)	
}


# extend HE parmeters for given number of terms
#   return vector in the form H1, H2, ..., E
he.rep <- function (x, n) {
    if (length(x) < 2) x <- rep(x, 2)
    x <- c(rep(x[-1], n)[1:n], x[1])
    return(x)
	}

last <- function(x) {x[length(x)]}

# copied from stats::: to avoid using :::
Pillai <- function (eig, q, df.res) 
{
	test <- sum(eig/(1 + eig))
	p <- length(eig)
	s <- min(p, q)
	n <- 0.5 * (df.res - p - 1)
	m <- 0.5 * (abs(p - q) - 1)
	tmp1 <- 2 * m + s + 1
	tmp2 <- 2 * n + s + 1
	c(test, (tmp2/tmp1 * test)/(s - test), s * tmp1, s * tmp2)
}

Wilks <- function (eig, q, df.res) 
{
	test <- prod(1/(1 + eig))
	p <- length(eig)
	tmp1 <- df.res - 0.5 * (p - q + 1)
	tmp2 <- (p * q - 2)/4
	tmp3 <- p^2 + q^2 - 5
	tmp3 <- if (tmp3 > 0) 
				sqrt(((p * q)^2 - 4)/tmp3)
			else 1
	c(test, ((test^(-1/tmp3) - 1) * (tmp1 * tmp3 - 2 * tmp2))/p/q, 
			p * q, tmp1 * tmp3 - 2 * tmp2)
}

HL <- function (eig, q, df.res) 
{
	test <- sum(eig)
	p <- length(eig)
	m <- 0.5 * (abs(p - q) - 1)
	n <- 0.5 * (df.res - p - 1)
	s <- min(p, q)
	tmp1 <- 2 * m + s + 1
	tmp2 <- 2 * (s * n + 1)
	c(test, (tmp2 * test)/s/s/tmp1, s * tmp1, tmp2)
}

Roy <- function (eig, q, df.res) 
{
	p <- length(eig)
	test <- max(eig)
	tmp1 <- max(p, q)
	tmp2 <- df.res - tmp1 + q
	c(test, (tmp2 * test)/tmp1, tmp1, tmp2)
}
