#
#	Watson-Wheeler test for homogeneity
#
#	Allows to compare the distribution of angles in two or more samples.
#
#	Based on
#		Circular statistics in biology, Batschelet, E (1981)
#		§6.3, p. 104
#		Biostatistical analysis, Zar, J H (1999)
#		§27.5, p. 640
#
# (c) Copyright 2010-211 Jean-Olivier Irisson
#     GNU General Public License, v.3
#
#------------------------------------------------------------

# Generic function
watson.wheeler.test <- function(x, ...) {
	UseMethod("watson.wheeler.test", x)
}

# Default method, for an angle vector and a grouping vector
watson.wheeler.test.default <- function(x, group, ...) {

	# get data name
	data.name <- paste(deparse(substitute(x)), "by", deparse(substitute(group)))

	# check arguments
	ok <- complete.cases(x, group)
	x <- x[ok]
	group <- group[ok]
	if (length(x)==0 | length(table(group)) < 2) {
		stop("No observations or no groups (at least after removing missing values)")
	}

	# remove circular attributes, if any
	if (is.circular(x)) {
		attr(x, "class") <- attr(x, "circularp") <- NULL
	}
	# NB: Since we will only work on ranks the units have no influence as long as they are consistent across groups. We do not need to store the angles attributes

	# check for ties
	nbTies <- sum(duplicated(x))
	if (nbTies > 0) {
		mess = ifelse(nbTies == 1, "There is 1 tie", paste("There are", nbTies, "ties"))
		warning(mess, " in the data.\n  Ties will be broken appart randomly and may influence the result.\n  Re-run the test several times to check the influence of ties.")
	}

	# check sample size per group
	ns <- as.numeric(table(group))
	if (!all(ns >= 10)) {
		warning("Some groups have less than 10 elements : ", paste(ns[ns < 10], collapse=", "),".\n  The Chi-squared approximation of the p-value is incorrect.")
	}

	result <- WatsonWheelerTestRad(x, group)
	result$data.name <- data.name

	return(result)
}

# Method for a list
watson.wheeler.test.list <- function(x, ...) {
	# fecth or fill list names
	k <- length(x)
	if (is.null(names(x))) {
		names(x) <- 1:k
	}
	# get data name
	data.name <- paste(names(x), collapse=" and ")

	# convert into x and group
	ns <- lapply(x, length)
	group <- rep(names(x), times=ns)
	x <- do.call("c", x)
	# NB: unlist() removes the circular attributes here

	# call default method
	result <- watson.wheeler.test.default(x, group)
	result$data.name <- data.name

	return(result)
}

# Method for a formula
watson.wheeler.test.formula <- function(formula, data, ...) {
	# convert into x and group
	d <- model.frame(as.formula(formula), data)

	# get data name
	data.name <- paste(names(d), collapse=" by ")

	# call default method
	result <- watson.wheeler.test.default(d[,1], d[,2])
	result$data.name <- data.name

	return(result)
}


# Computation in the usual trigonometric space
WatsonWheelerTestRad <- function(x, group) {

	# number of groups
	group <- as.factor(group)
	k <- nlevels(group)
	# total sample size
	n <- length(x)
	# sample size per group
	ns <- as.numeric(table(group))

	# ranks
	r <- rank(x, ties.method="random")
	# circular rank (or uniform score)
	cr <- r * 2*pi / n

	# compute
	C <- tapply(cos(cr), group, sum)
	S <- tapply(sin(cr), group, sum)

	if (k == 2) {
		W <- 2 * (n-1) * (C[1]^2 + S[1]^2) / prod(ns)
		names(W) <- NULL
		df <- 2
	} else {
		W <- 2 * sum( (C^2 + S^2) / ns)
		df <- 2*(k-1)
	}

	p.value <- pchisq(W, df=df, lower.tail=FALSE)

	# return result
	result <- list(
		method = "Watson-Wheeler test for homogeneity of angles",
		parameter = c(df=df),
		statistic = c(W=W),
		p.value = p.value
	)
	class(result) <- "htest"

	return(result)
}


# Test data (from Zar)
# x1 <- c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120)
# x2 <- c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165)
#
# watson.wheeler.test(list(x1,x2))
