"f.fb.neglog.deriv" <-
function(t, mu, tau, c1)
{
	.expr3 <- (sqrt((2 * pi))) * t
	.expr4 <- t^2
	.expr7 <- (.expr4 * (tau^2)) + t
	.expr9 <- .expr3 * (sqrt(.expr7))
	.expr10 <- c1/.expr9
	.expr12 <- c1 - (mu * t)
	.expr13 <- .expr12^2
	.expr15 <- 2 * .expr7
	.expr17 <- exp((( - .expr13)/.expr15))
	.expr18 <- .expr10 * .expr17
	.expr22 <- 2 * (t * .expr12)
	.expr23 <- .expr22/.expr15
	.expr24 <- .expr17 * .expr23
	.expr25 <- .expr10 * .expr24
	.expr37 <- .expr18^2
	.expr42 <- .expr4 * (2 * tau)
	.expr43 <- 2 * .expr42
	.expr44 <- .expr13 * .expr43
	.expr45 <- .expr15^2
	.expr46 <- .expr44/.expr45
	.expr47 <- .expr17 * .expr46
	.expr54 <- .expr7^-0.5
	.expr57 <- .expr3 * (0.5 * (.expr42 * .expr54))
	.expr58 <- c1 * .expr57
	.expr59 <- .expr9^2
	.expr60 <- .expr58/.expr59
	.expr66 <- (.expr10 * .expr47) - (.expr60 * .expr17)
	.expr70 <-  - ((((.expr10 * ((.expr47 * .expr23) - (.expr17 * ((.expr22 *
		.expr43)/.expr45)))) - (.expr60 * .expr24))/.expr18) - ((
		.expr25 * .expr66)/.expr37))
	.expr71 <- 1/.expr9
	.expr76 <- 2 * .expr12
	.expr77 <- .expr76/.expr15
	.expr78 <- .expr17 * .expr77
	.expr86 <- (.expr71 * .expr17) - (.expr10 * .expr78)
	.expr90 <-  - ((((.expr71 * .expr24) + (.expr10 * ((.expr17 * ((2 * t)/
		.expr15)) - (.expr78 * .expr23))))/.expr18) - ((.expr25 * 
		.expr86)/.expr37))
	.expr94 <- .expr4 * 2
	.expr107 <- .expr60 * .expr47
	.expr150 <-  - (((((.expr71 * .expr47) + (.expr10 * ((.expr17 * ((
		.expr76 * .expr43)/.expr45)) - (.expr78 * .expr46)))) - (((
		.expr57/.expr59) * .expr17) - (.expr60 * .expr78)))/.expr18) - (
		(.expr66 * .expr86)/.expr37))
	.expr153 <- .expr71 * .expr78
	.value <-  - (log(.expr18))
	.grad <- array(0, c(length(.value), 3), list(NULL, c("mu", "tau", "c1")
		))
	.hess <- array(0, c(length(.value), 3, 3), list(NULL, c("mu", "tau", 
		"c1"), c("mu", "tau", "c1")))
	.grad[, "mu"] <-  - (.expr25/.expr18)
	.grad[, "tau"] <-  - (.expr66/.expr18)
	.grad[, "c1"] <-  - (.expr86/.expr18)
	.hess[, "mu", "mu"] <-  - (((.expr10 * ((.expr24 * .expr23) - (.expr17 * (
		(2 * (t * t))/.expr15))))/.expr18) - ((.expr25 * .expr25)/
		.expr37))
	.hess[, "tau", "mu"] <- .expr70
	.hess[, "c1", "mu"] <- .expr90
	.hess[, "mu", "tau"] <- .expr70
	.hess[, "tau", "tau"] <-  - (((((.expr10 * ((.expr47 * .expr46) + (
		.expr17 * (((.expr13 * (2 * .expr94))/.expr45) - ((.expr44 * (2 *
		(.expr43 * .expr15)))/(.expr45^2)))))) - .expr107) - (((((c1 * (
		.expr3 * (0.5 * ((.expr94 * .expr54) + (.expr42 * (-0.5 * (
		.expr42 * (.expr7^-1.5))))))))/.expr59) - ((.expr58 * (2 * (
		.expr57 * .expr9)))/(.expr59^2))) * .expr17) + .expr107))/
		.expr18) - ((.expr66 * .expr66)/.expr37))
	.hess[, "c1", "tau"] <- .expr150
	.hess[, "mu", "c1"] <- .expr90
	.hess[, "tau", "c1"] <- .expr150
	.hess[, "c1", "c1"] <- ((.expr153 + (.expr153 + (.expr10 * ((.expr17 * (
		2/.expr15)) - (.expr78 * .expr77)))))/.expr18) + ((.expr86 * 
		.expr86)/.expr37)
	attr(.value, "gradient") <- .grad
	attr(.value, "hessian") <- .hess
	.value
}

