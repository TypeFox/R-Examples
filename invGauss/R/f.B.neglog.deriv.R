"f.B.neglog.deriv" <-
function(t, mu, tau, c1)
{
	.expr1 <- mu * t
	.expr2 <- c1 - .expr1
	.expr3 <- t^2
	.expr4 <- tau^2
	.expr6 <- (.expr3 * .expr4) + t
	.expr7 <- sqrt(.expr6)
	.expr8 <- .expr2/.expr7
	.expr10 <- 2 * c1
	.expr13 <- 2 * (c1^2)
	.expr16 <- exp(((.expr10 * mu) + (.expr13 * .expr4)))
	.expr18 <- .expr10 * t
	.expr21 <- (( - c1) - (.expr18 * .expr4)) - .expr1
	.expr22 <- .expr21/.expr7
	.expr23 <- pnorm(.expr22)
	.expr25 <- (pnorm(.expr8)) - (.expr16 * .expr23)
	.expr28 <- dnorm(.expr8)
	.expr29 <- t/.expr7
	.expr31 <- .expr16 * .expr10
	.expr33 <- dnorm(.expr22)
	.expr34 <- .expr33 * .expr29
	.expr37 <- (.expr28 * .expr29) + ((.expr31 * .expr23) - (.expr16 * 
		.expr34))
	.expr39 <- .expr8 * .expr28
	.expr44 <- .expr31 * .expr34
	.expr46 <- .expr22 * .expr33
	.expr55 <- .expr25^2
	.expr58 <- 2 * tau
	.expr59 <- .expr3 * .expr58
	.expr60 <- .expr6^-0.5
	.expr62 <- 0.5 * (.expr59 * .expr60)
	.expr63 <- .expr2 * .expr62
	.expr64 <- .expr7^2
	.expr65 <- .expr63/.expr64
	.expr66 <- .expr39 * .expr65
	.expr69 <- (t * .expr62)/.expr64
	.expr72 <- .expr13 * .expr58
	.expr73 <- .expr16 * .expr72
	.expr76 <- .expr18 * .expr58
	.expr78 <- .expr21 * .expr62
	.expr80 <- (.expr76/.expr7) + (.expr78/.expr64)
	.expr81 <- .expr33 * .expr80
	.expr85 <- .expr46 * .expr80
	.expr98 <- (.expr28 * .expr65) + ((.expr73 * .expr23) - (.expr16 * 
		.expr81))
	.expr101 <- ((((.expr66 * .expr29) - (.expr28 * .expr69)) + ((((.expr73 *
		.expr10) * .expr23) - (.expr31 * .expr81)) - ((.expr73 * 
		.expr34) + (.expr16 * ((.expr85 * .expr29) - (.expr33 * .expr69
		))))))/.expr25) + ((.expr37 * .expr98)/.expr55)
	.expr103 <- 2 * .expr10
	.expr105 <- (2 * mu) + (.expr103 * .expr4)
	.expr106 <- .expr16 * .expr105
	.expr111 <- 2 * t
	.expr113 <- 1 + (.expr111 * .expr4)
	.expr114 <- .expr113/.expr7
	.expr115 <- .expr33 * .expr114
	.expr119 <- .expr46 * .expr114
	.expr124 <- 1/.expr7
	.expr125 <- .expr39 * .expr124
	.expr133 <- (.expr28 * .expr124) - ((.expr106 * .expr23) - (.expr16 * 
		.expr115))
	.expr136 <- (((((((.expr106 * .expr10) + (.expr16 * 2)) * .expr23) - (
		.expr31 * .expr115)) - ((.expr106 * .expr34) + (.expr16 * (
		.expr119 * .expr29)))) - (.expr125 * .expr29))/.expr25) - ((
		.expr37 * .expr133)/.expr55)
	.expr146 <- 0.5 * (((.expr3 * 2) * .expr60) + (.expr59 * (-0.5 * (
		.expr59 * (.expr6^-1.5)))))
	.expr150 <- 2 * (.expr62 * .expr7)
	.expr152 <- .expr64^2
	.expr162 <- .expr73 * .expr81
	.expr167 <- .expr76 * .expr62
	.expr214 <- ((((.expr28 * (.expr62/.expr64)) - (.expr125 * .expr65)) + (
		((((.expr106 * .expr72) + (.expr16 * (.expr103 * .expr58))) * 
		.expr23) - (.expr73 * .expr115)) - ((.expr106 * .expr81) + (
		.expr16 * ((.expr119 * .expr80) + (.expr33 * (((.expr111 * 
		.expr58)/.expr7) - ((.expr113 * .expr62)/.expr64))))))))/
		.expr25) - ((.expr98 * .expr133)/.expr55)
	.expr223 <- .expr106 * .expr115
	.value <-  - (log(.expr25))
	.grad <- array(0, c(length(.value), 3), list(NULL, c("mu", "tau", "c1")
		))
	.hess <- array(0, c(length(.value), 3, 3), list(NULL, c("mu", "tau", 
		"c1"), c("mu", "tau", "c1")))
	.grad[, "mu"] <- .expr37/.expr25
	.grad[, "tau"] <- .expr98/.expr25
	.grad[, "c1"] <-  - (.expr133/.expr25)
	.hess[, "mu", "mu"] <- ((((.expr39 * .expr29) * .expr29) + ((((.expr31 * 
		.expr10) * .expr23) - .expr44) - (.expr44 + (.expr16 * ((
		.expr46 * .expr29) * .expr29)))))/.expr25) + ((.expr37 * 
		.expr37)/.expr55)
	.hess[, "tau", "mu"] <- .expr101
	.hess[, "c1", "mu"] <- .expr136
	.hess[, "mu", "tau"] <- .expr101
	.hess[, "tau", "tau"] <- ((((.expr66 * .expr65) + (.expr28 * (((.expr2 * 
		.expr146)/.expr64) - ((.expr63 * .expr150)/.expr152)))) + (((((
		.expr73 * .expr72) + (.expr16 * (.expr13 * 2))) * .expr23) - 
		.expr162) - (.expr162 + (.expr16 * ((.expr85 * .expr80) + (
		.expr33 * ((((.expr18 * 2)/.expr7) - (.expr167/.expr64)) + ((((
		.expr21 * .expr146) - .expr167)/.expr64) - ((.expr78 * .expr150
		)/.expr152)))))))))/.expr25) + ((.expr98 * .expr98)/.expr55)
	.hess[, "c1", "tau"] <- .expr214
	.hess[, "mu", "c1"] <- .expr136
	.hess[, "tau", "c1"] <- .expr214
	.hess[, "c1", "c1"] <- (((.expr125 * .expr124) + (((((.expr106 * 
		.expr105) + (.expr16 * (4 * .expr4))) * .expr23) - .expr223) - (
		.expr223 + (.expr16 * (.expr119 * .expr114)))))/.expr25) + ((
		.expr133 * .expr133)/.expr55)
	attr(.value, "gradient") <- .grad
	attr(.value, "hessian") <- .hess
	.value
}

