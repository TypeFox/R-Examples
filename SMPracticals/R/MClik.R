"MClik" <-
function(d)
{
	max.ord <- 3
	n <- length(d)
	i <- c(1:(n - max.ord))	# third-order model 
	four <- table(d[i], d[i + 1], d[i + 2], d[i + max.ord])
	ts <- apply(four, c(1, 2, 3), sum)
	ts <- sweep(four, c(1, 2, 3), ts, "/")
	L3 <- four * log(ts)
	L3 <- sum(ifelse(is.na(L3), 0, L3))	# second-order model 
	three <- apply(four, c(2, 3, 4), sum)
	ts <- apply(three, c(1, 2), sum)
	ts <- sweep(three, c(1, 2), ts, "/")
	L2 <- three * log(ts)
	L2 <- sum(ifelse(is.na(L2), 0, L2))	# first-order model
	two <- apply(three, c(2, 3), sum)
	ts <- apply(two, 1, sum)
	ts <- sweep(two, 1, ts, "/")
	L1 <- two * log(ts)
	L1 <- sum(ifelse(is.na(L1), 0, L1))	# independence model
	one <- apply(two, 2, sum)
	ts <- sum(one)
	L0 <- sum(one * log(one/ts))	# end calcuations
	S <- length(one)
	df <- (S - 1) * S^c(0:max.ord)
	L <- c(L0, L1, L2, L3)
	AIC <- -2 * (L - df)
	list(order = 0:3, df = df, L = L, AIC = AIC, one = one, two = two, three = three, 
		four = four)
}

