local_slope <-
function(p, k) {
	n <- length(p)
	o <- ordered_values(n)
	z <- -log(p[order(p, decreasing = TRUE)])
	sum1 <- 0
	for (j in 2:k) {
		sum1 <- sum1 + z[n - j]
		}
	eta <- (1/k) * ((2 * z[n - 1]) + sum1 - ((k + 1) * z[n - k + 1]))
	test_statistic <- (z[n] - z[n - 1])/eta
	Fvalue <- qf(0.95, df1 = 2, df2 = 2*k, lower.tail = TRUE)
	pvalue <- pf(test_statistic, df1 = 2, df2 = 2*k, lower.tail = FALSE)
	return(list("local_slope" = eta, "test_statistic" = test_statistic, "Fvalue" = Fvalue, "pvalue" = pvalue))
	}
