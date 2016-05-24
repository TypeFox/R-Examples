MultiStart <- function(obj, repeats, draw, spread, debug) {

	mFn <- c(1, 1, P3, 1, P5c, 1, P7c) # not sure if this is implemented !


	if (missing(repeats)) 
		repeats = 400
	if (missing(draw)) 
		draw = FALSE
	if (missing(spread)) 
		spread = 0.15
	if (missing(debug)) 
		debug = F


	if (debug) 
		print("+++ missing values assigned OK")

	resid = NULL
	fit = NULL
	val = NULL
	Pn = NULL
	AIC = NULL
	if (is.list(obj)) { # checks that obj is a list but doesnt report failure todo
		Res <- obj
		Res$call = NULL
		x <- obj$time
		y <- obj$thrs
		p <- obj$opt
		val <- obj$val
		if (is.null(obj$opt)) 
			p = obj$init
		Pn <- obj$Pn
		if (is.null(obj$Pn)) 
			Pn = 7
		AIC <- obj$AIC[1:7]
	}

	.GlobalEnv$x <- x
	.GlobalEnv$y <- y

	if (debug) 
		print("+++ object processed OK")
	Fn <- mFn[[Pn]]

	if (debug) 
		print(paste("+++ Using the ", Pn, " parameter model ", sep = ""))

	OptJK <- function(a) {
		tmp <- numeric(9)
		# two iterations of the optim fn could allow 1000 cycles, will compare later
		X = optim(a, Fn)
		X = optim(X$par, Fn)
		X = optim(X$par, Fn)
		# X = optim(X$par, Fn)
		tmp[1:Pn] = X$par
		tmp[8] = X$val
		tmp[9] = X$con
		tmp
	}

	Par <- matrix(p * rnorm(7 * repeats, 1, spread), 7, repeats)

	O <- t(apply(Par, 2, OptJK))


	input <- numeric(9)
	input[1:7] = p
	input[8] = val
	input[9] = 0
	O <- rbind(input, O)
	if (debug) 
		print(head(O))

	Test <- sum(O[, 9] == 0)
	if (debug) 
		print(paste("+++ the test has boolean value ", Test, sep = ""))
	if (Test) {
		idx <- O[, 9] == 0
		O <- O[idx, ]
	}
	if (debug) 
		print(head(idx))

	if (length(O) != 9) {
		idx <- order(O[, 8])
		O <- O[idx, ]
		val <- O[1, 8]
		p <- O[1, 1:Pn]
	} else {
		val <- O[8]
		p <- O[1:Pn]
	}
	if (length(O) == 9) 
		Res$warning <- "+++ None of the jittered values converged"

	if (debug) 
		print("+++ Ordered index made")


	val = val[[1]]

	fit <- Fn(p, x)
	resid <- (y - fit)

	if (draw) {
		plot(x, y)
		lines(x, fit)
	}
	#### create output object 
	Res$call <- match.call()
	Res$opt = p
	Res$time = x
	Res$thrs = y
	Res$resid = resid
	Res$fit = fit
	Res$val = val
	Res$data = obj$data
	Res$Mod = obj$Mod
	Res$Pn = obj$Pn
	Res$AIC = AIC
	Res$R2 <- 1 - (var(resid)/var(y))

	if (debug) 
		Res$O <- O

	#### clear up the global variables we created    
	on.exit(rm(list = c("x", "y"), envir = .GlobalEnv))
	class(Res) = "dark"
	return(Res)
}
