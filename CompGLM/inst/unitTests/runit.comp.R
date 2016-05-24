
test.dcomp <- function() {
	checkEqualsNumeric(dpois(-5:5, 2.5, FALSE), dcomp(-5:5, 2.5, 1, 100L, FALSE))
	checkEqualsNumeric(dpois(1:2, 2:5, FALSE), dcomp(1:2, 2:5, 1, 100L, FALSE))
	checkEqualsNumeric(dpois(2:5, 1:2, FALSE), dcomp(2:5, 1:2, 1, 100L, FALSE))

	checkEqualsNumeric(dpois(-5:5, 2.5, TRUE), dcomp(-5:5, 2.5, 1, 100L, TRUE))
	checkEqualsNumeric(dpois(1:2, 2:5, TRUE), dcomp(1:2, 2:5, 1, 100L, TRUE))
	checkEqualsNumeric(dpois(2:5, 1:2, TRUE), dcomp(2:5, 1:2, 1, 100L, TRUE))
	
	checkException(dcomp(1, -1, 2), silent = TRUE)
	checkException(dcomp(1, 1, -2), silent = TRUE)
}

test.pcomp <- function() {
	checkEqualsNumeric(ppois(-5:5, 2.5, TRUE, FALSE), pcomp(-5:5, 2.5, 1, 100L, TRUE, FALSE))
	checkEqualsNumeric(ppois(1:2, 2:5, TRUE, FALSE), pcomp(1:2, 2:5, 1, 100L, TRUE, FALSE))
	checkEqualsNumeric(ppois(2:5, 1:2, TRUE, FALSE), pcomp(2:5, 1:2, 1, 100L, TRUE, FALSE))

	checkEqualsNumeric(ppois(-5:5, 2.5, FALSE, TRUE), pcomp(-5:5, 2.5, 1, 100L, FALSE, TRUE))
	checkEqualsNumeric(ppois(1:2, 2:5, FALSE, FALSE), pcomp(1:2, 2:5, 1, 100L, FALSE, FALSE))
	checkEqualsNumeric(ppois(1:2, 2:5, TRUE, TRUE), pcomp(1:2, 2:5, 1, 100L, TRUE, TRUE))

	checkException(pcomp(1, -1, 2), silent = TRUE)
	checkException(pcomp(1, 1, -2), silent = TRUE)
}

test.glm.comp <- function() {
	set.seed(1)
	n <- 5000
	x1 <- rnorm(n, -1.0, 0.5)
	x2 <- rnorm(n, 1.0, 0.7)
	x3 <- rnorm(n, 2.0, 0.4)
	y <- rpois(n, exp(-0.5 + 0.3 * x1 + 0.8 * x2 + 0.2 * x3))
	
	data <- data.frame(y, x1, x2, x3)
	
	poissonModel <- glm(y ~ ., poisson, data)
	compModel <- glm.comp(y ~ ., data = data)
	
	checkEqualsNumeric(coef(poissonModel), coef(compModel)$beta, "", 0.1)
	checkEqualsNumeric(0.0, coef(compModel)$zeta, "", 0.1)
}
