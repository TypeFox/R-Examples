Approx.cfpt.density <-
function (sfptl, variableStep = TRUE, from.t0 = FALSE, to.T = FALSE, skip = TRUE, n = 250, p = 0.2, alpha = 1, tol = 0.001, it.max = 50000L)
{
	if (!is.summary.fptl(sfptl)) 
      	stop(paste(sQuote("sfptl"), " object is not of class ", shQuote("summary.fptl"), ".", sep = ""))

	if (length(sfptl) > 1L)
		stop(paste(sQuote("sfptl"), "does not match a degenerate initial distribution."))

	if ((!missing(p)) & (!missing(alpha))) 
      	IS <- Integration.Steps(sfptl, variableStep, from.t0, to.T, n, p, alpha)
	else {
      	if (missing(p) & (!missing(alpha))) 
            	IS <- Integration.Steps(sfptl, variableStep, from.t0, to.T, n, , alpha)
        	else {
            	if ((!missing(p)) & missing(alpha)) 
                		IS <- Integration.Steps(sfptl, variableStep, from.t0, to.T, n, p)
            	else IS <- Integration.Steps(sfptl, variableStep, from.t0, to.T, n)
        	}
	}

	H <- IS$H
	Skip <- IS$skip

      Args <- as.list(attr(sfptl, "FPTLCall")[[1]])
      t0 <- Args$t0
      x0 <- Args$x0
	env <- Args$env
	e <- !is.null(env)
      if (e){
		if (!is.list(env)) env <- eval(env, attr(sfptl, "vars"))
		logic <- unlist(lapply(env, is.character))
      	if (any(logic))
			env[logic] <- lapply(env[logic], function(x) parse(text=x)[[1]])
	}

	it <- (H[, 2] - H[, 1])/H[, 3]     

      if ((sum(it) > it.max) & interactive()) {
        	label1 <- character(1)
		label2 <- character(1)
		if (!to.T){
			logic <- !unlist(Skip)
            	if (any(logic)){
				itStop <- cumsum(it)[logic]			
				cat("\nHowever,", ifelse(skip, " if no subinterval were avoided,", ""), " the algorithm can stop at ", sep = "")
                		if (length(itStop) == 1) label1 <- paste("iteration", itStop) else label1 <- paste("iterations ", paste(itStop[-length(itStop)], collapse = ", "), " and ", itStop[length(itStop)], sep = "")
				label1 <- paste(", although the algorithm can stop at", label1)
				itStop <- cumsum(it[logic])
				if (length(itStop) == 1) label2 <- paste("iteration", itStop) else label2 <- paste("iterations ", paste(itStop[-length(itStop)], collapse = ", "), " and ", itStop[length(itStop)], sep = "")
				label2 <- paste(", although the algorithm can stop at", label2)
            	}            	
        	}
        	itbyStep <- format(data.frame(matrix(H[, 1:2], ncol = 2), it), digits = 8)
        	logic <- unlist(Skip)
        	if ((skip) & (any(logic)))
            	itbyStep[, 3][logic] <- paste("0 or", itbyStep[, 3][logic])
        	row.names(itbyStep) <- paste("Subinterval", 1:nrow(itbyStep))
        	names(itbyStep) <- c("Lower end", "Upper end", "Iterations")
        	cat("\nThe table below shows the number of iterations of the approximation process step by step: \n\n")
        	print(itbyStep)
        	iterations <- sum(it)
		if ((skip) & (any(logic))){
			cat("\nIf no subinterval is avoided, the total number of iterations is ", iterations, label1, ".", sep = "")
			cat("\nIf the application of the numerical algorithm is avoided in all suitable subintervals, the total number of iterations is ", sum(it[!logic]), label2, ".", sep = "") 
		}
		else cat("\nThe total number of iterations is", iterations)
		cat("\n")       	
        	repeat {
            	answer <- readline("Do you want to continue? (y/n) ")
            	if ((answer == "n") | (answer == "y")) break
        	}
        	if (answer == "n")
            	stop("Approximation process stopped by the usuary.\n\n")

    	}
    
	exprS <- parse(text = Args$S)
      if (inherits(try(exprS, silent = TRUE), "try-error"))
		stop(paste("the mathematical expression of the boundary shows syntax errors.", sep = ""))
      if (inherits(try(D(exprS, "t"), silent = TRUE), "try-error"))
		stop("R can not compute the symbolic derivative with respect to 't' of the mathematical expression of the boundary")  
	
	dp <- attr(sfptl, "dp")

	if (e){
		A1 <- function(x, t) NULL
		expr <- as.expression(eval(substitute(substitute(expr, env), list(expr = parse(text = dp$mean)[[1]]))))
		if (any(is.element(c("x", "t"), all.names(expr)))){
      		vars <- setdiff(intersect(all.names(expr), names(env)), c("x", "t"))
      		if (length(vars) > 0) expr <- c(sapply(vars, function(x) as.expression(eval(substitute(substitute(expr, list(valor=env[[x]])), list(expr = parse(text = paste(x, "<- valor"))[[1]]))))), expr)
      		body(A1) <- as.call(c(as.name("{"), expr))
		}
		else body(A1) <- substitute(rep(a, max(length(x),length(t))), list(a = expr[[1]]))

      	DA2 <- function(x, t) NULL
		expr <- parse(text = dp$var)
		if (any(is.element(c("x", "t"), all.names(expr)))){
			exprD <- as.expression(eval(substitute(substitute(expr, env), list(expr = deriv(expr, "x")[[1]]))))
      		vars <- setdiff(intersect(all.names(exprD), names(env)), c("x", "t"))
      		if (length(vars) > 0) exprD <- c(sapply(vars, function(x) as.expression(eval(substitute(substitute(expr, list(valor=env[[x]])), list(expr = parse(text = paste(x, "<- valor"))[[1]]))))), exprD)
      		body(DA2) <- as.call(c(as.name("{"), exprD))
		}
		else body(DA2) <- as.call(c(as.name("{"), eval(substitute(expression({
    						.value <- rep(a, length(x))
    						.grad <- array(0, c(length(.value), 1L), list(NULL, c("x")))
    						.grad[, "x"] <- 0
    						attr(.value, "gradient") <- .grad
    						.value
						}), list(a = eval(substitute(substitute(expr, env), list(expr = expr[[1]]))))))))

      	Df <- function(x, t, y, s) NULL
		expr <- parse(text = dp$tpdf)
      	exprD <- as.expression(eval(substitute(substitute(expr, env), list(expr = deriv(expr, "x")[[1]]))))
      	vars <- setdiff(intersect(all.names(exprD), names(env)), c("x", "t", "y", "s"))
      	if (length(vars) > 0) exprD <- c(sapply(vars, function(x) as.expression(eval(substitute(substitute(expr, list(valor=env[[x]])), list(expr = parse(text = paste(x, "<- valor"))[[1]]))))), exprD)
      	body(Df) <- as.call(c(as.name("{"), exprD))

		B <- function(t) NULL
      	env <- c(t0 = t0, x0 = x0, env)
		exprS <- as.expression(eval(substitute(substitute(expr, env), list(expr = exprS[[1]]))))
		if (is.element("t", all.names(exprS))){	
			vars <- setdiff(intersect(all.names(exprS), names(env)), "t")
      		if (length(vars) > 0) expr <- c(sapply(vars, function(x) as.expression(eval(substitute(substitute(expr, list(valor=env[[x]])), list(expr = parse(text = paste(x, "<- valor"))[[1]]))))), exprS)
      		body(B) <- as.call(c(as.name("{"), exprS))
		}
		else body(B) <- substitute(rep(a, length(t)), list(a = eval(substitute(substitute(expr, env), list(expr = exprS[[1]])))))

		DB <- function(t) NULL
      	if (is.element("t", all.names(exprS))){
			exprD <- deriv(exprS, "t")
      		vars <- setdiff(intersect(all.names(exprD), names(env)), "t")
      		if (length(vars) > 0) exprD <- c(sapply(vars, function(x) as.expression(eval(substitute(substitute(expr, list(valor=env[[x]])), list(expr = parse(text = paste(x, "<- valor"))[[1]]))))), exprD)
      		body(DB) <- as.call(c(as.name("{"), exprD))
		}
		else body(DB) <- as.call(c(as.name("{"), eval(substitute(expression({
    						.value <- rep(a, length(t))
    						.grad <- array(0, c(length(.value), 1L), list(NULL, c("t")))
    						.grad[, "t"] <- 0
    						attr(.value, "gradient") <- .grad
    						.value
						}), list(a = eval(substitute(substitute(expr, env), list(expr = exprS[[1]]))))))))
	}
	else{
		A1 <- function(x, t) NULL
		expr <- parse(text = dp$mean)
		if (any(is.element(c("x", "t"), all.names(expr)))) body(A1) <- expr
		else body(A1) <- substitute(rep(a, max(length(x),length(t))), list(a = expr[[1]]))

		DA2 <- function(x, t) NULL
		expr <- parse(text = dp$var)
      	if (any(is.element(c("x", "t"), all.names(expr)))) body(DA2) <- deriv(expr, "x")
		else body(DA2) <- as.call(c(as.name("{"), eval(substitute(expression({
    						.value <- rep(a, length(x))
    						.grad <- array(0, c(length(.value), 1L), list(NULL, c("x")))
    						.grad[, "x"] <- 0
    						attr(.value, "gradient") <- .grad
    						.value
						}), list(a = expr[[1]])))))

		Df <- function(x, t, y, s) NULL
      	expr <- deriv(parse(text = dp$tpdf), "x")
      	body(Df) <- expr

		B <- function(t) NULL
      	env <- list(t0 = t0, x0 = x0)
      	exprS <- as.expression(eval(substitute(substitute(expr, env), list(expr = exprS[[1]]))))
		if (is.element("t", all.names(expr))) body(B) <- exprS
		else body(B) <- substitute(rep(a, length(t)), list(a = exprS[[1]]))

		DB <- function(t) NULL
      	if (is.element("t", all.names(exprS))) body(DB) <- deriv(exprS, "t")
		else body(DB) <- as.call(c(as.name("{"), eval(substitute(expression({
    						.value <- rep(a, length(t))
    						.grad <- array(0, c(length(.value), 1L), list(NULL, c("t")))
    						.grad[, "t"] <- 0
    						attr(.value, "gradient") <- .grad
    						.value
						}), list(a = exprS[[1]])))))

	}
	
	g <- alist(x = , y = , y.x0 = NULL)

	Call <- match.call()
	Args <- list(variableStep = variableStep, from.t0 = from.t0, to.T = to.T, skip = skip, n = n, p = p, alpha = alpha, tol = tol, it.max = it.max)
	label <- intersect(names(Args), names(Call))
	if (length(label) > 0L) Call[label] <- Args[label]	

	attr(g, "Call") <- Call
	attr(g, "Steps") <- H    					

      rho <- sign(B(t0) - x0)
      if (rho == 0) stop("S(t0) = x0")

	jumps <- replicate(nrow(H), integer(0))
      Integral <- numeric(nrow(H))
	now <- matrix(, nrow = nrow(H) + 1, ncol = 2)

	cat("\nComputing...\t")
      now[1, ] <- proc.time()[1:2]

	h <- numeric(0)
	F <- 0
      g0 <- H[1, 1]
      g$x <- H[1, 1] + H[1, 3]
      b <- DB(g$x)
      a2 <- DA2(b, g$x)
      ff <- Df(b, g$x, x0, t0)
      g$y <- max(0, - rho * (ff * (attr(b, "gradient") - A1(b, g$x) +  0.75 * attr(a2, "gradient")) + attr(ff, "gradient") * a2))
      Integral[1] <- H[1, 3] * g$y/2

      H[1, 1] <- g$x

      for (i in 1:nrow(H)) {
      	if (skip & Skip[[i]] & (g$y[length(g$y)] < 10^-10)) {
            	correction <- h[length(h)] * g$y[length(g$y)]/2
            	Integral[i - 1] <- Integral[i - 1] - correction
            	F <- F - correction
            	g$y[length(g$y)] <- 0
            	h <- c(h, H[i, 2] - H[i, 1])
            	g$x <- c(g$x, H[i, 2])
            	b <- c(b, DB(H[i, 2]))
            	g$y <- c(g$y, 0)
            	now[i + 1, ] <- proc.time()[1:2]
            	Integral[i] <- 0
            	jumps[[i]] <- 1L
        	}
        	else {
            	u2 <- seq(H[i, 1] + H[i, 3], H[i, 2], by = H[i, 3])
            	h2 <- rep(H[i, 3], length(u2))
            	b2 <- DB(u2)
            	a2 <- DA2(b2, u2)
            	a <- as.vector(attr(b2, "gradient")) - A1(b2, u2) + 0.75 * as.vector(attr(a2, "gradient"))
            	f0 <- Df(b2, u2, x0, t0)
            	g$y <- c(g$y, -rho * (f0 * a + as.vector(attr(f0, "gradient")) * a2))
            	n <- length(g$x)
            	h <- c(h, h2)
            	g$x <- c(g$x, u2)
            	if (length(b2) == 1) b2 <- rep(b2, length(u2))
            	b <- c(b, b2)
            	for (k in (1:length(u2))) {
                		index <- 1:(n + k - 1)
                		f1 <- Df(b2[k], u2[k], b[index], g$x[index])
                		g$y[n + k] <- max(0, g$y[n + k] + rho * sum(h[index] * g$y[index] * (f1 * a[k] + as.vector(attr(f1, "gradient")) * a2[k])))
            	}
            	Integral[i] <- Integral[i] + h2[1] * (sum(g$y[n:(n + length(u2))]) - (g$y[n] + g$y[n + length(u2)])/2)

            	F <- F + Integral[i]
            	now[i + 1, ] <- proc.time()[1:2]            	
           	}
		if ((!to.T) & (!Skip[[i]]) & (F >= (1 - tol))) break		
    	}

	cat("Done.\n")
      if (F < 1 - tol) {
      	cat(paste("\nThe value of the cumulative integral of the approximation is", F, "< 1 - tol.\n"))
      	if (to.T)
            	cat("It may be appropriate to extend the considered time interval.\n")
        	else {
            	cat("If the value of the cumulative integral is not high and the final stopping instant is less than T, it may be appropriate:")
            	cat("\n   - Check if the value of the final stopping instant increases using k argument to summary the fptl class object, or")
            	cat("\n   - Approximate the density again with to.T = TRUE.\n")
        	}
      }

      g$x <- c(g0, g$x)
      g$y <- c(0, g$y)
      if (i < nrow(H)) {
      	Integral <- Integral[1:i]
        	jumps <- jumps[1:i]
        	now <- now[1:(i + 1), ]
    	}
    	CPUTime <- apply(now, 2, diff)
    	if (!is.matrix(CPUTime)) CPUTime <- matrix(CPUTime, ncol = 2)
	attr(g, "cumIntegral") <- cumsum(Integral)
    	attr(g, "skips") <- jumps
    	attr(g, "CPUTime") <- CPUTime
	attr(g, "summary.fptl") <- sfptl
    	class(g) <- c("fpt.density", "list")
    	return(g)
}
