Approx.fpt.density <-
function (dp, t0, T, id, S, env = NULL, variableStep = TRUE, from.t0 = FALSE, to.T = FALSE, r = 4000, zeroSlope = 0.01, p0.tol = 8, k = 3, 
m = 100, n = 250, p = 0.2, alpha = 1, skip = TRUE, tol = 0.001, it.max)
{
	if (!is.diffproc(dp)) stop(paste(sQuote("dp"), "argument is not of class ", shQuote("diffproc"), ".", sep = ""))
      if (!is.numeric(t0)) stop(paste(sQuote("t0"), "argument is not of numeric type."))
    	if (length(t0) > 1) stop(paste(sQuote("t0"), "argument is not of length 1."))
    	if (!is.numeric(T)) stop(paste(sQuote("T"), "argument is not of numeric type."))
    	if (length(T) > 1) stop(paste(sQuote("T"), "argument is not of length 1."))
    	if (t0 >= T) stop("the final time instant is not greater than the initial time instant.")
	if (m < 25) {
      	warning("m can be small, the value 25 has been used", call. = FALSE)
      	m <- 25
    	}
	if (is.numeric(id)){
		if (length(id) > 1L) stop("id argument is not of length 1.")
	}
	else{
		if (is.list(id)){
			if (length(id) != 4L) stop("id argument is not of length 4.")
			if (!all(sapply(id, is.character))) stop("the components of id argument are not of character type.")
			if (!all(sapply(id, length) == 1L)) stop("the components of id argument are not of length 1.")
		}
		else stop(paste(sQuote("id"), "argument is not of numeric or list type."))
	}	
	if (is.numeric(id)){
		X0 <- id		
		df.X0 <- 1L
		cp <- TRUE
	}
	else{
		if (inherits(try(parse(text = id[[1]]), silent = TRUE), "try-error")) 
      		stop("the mathematical expression of the density function of the initial distribution process shows syntax errors.")
		expr <- parse(text = id[[1]])[[1]]
		v <- all.vars(expr)
		if(!identical(v, "x")) 
			stop("the mathematical expression of the density function of the initial distribution process is not an expression in 'x'.")
		distr <- substring(ls("package:stats", pattern="^d"),2)
		distr <- paste("d", intersect(distr, substring(ls("package:stats", pattern="^q"),2)), sep = "")
		if (!is.element(as.character(expr[[1]]), distr)) 
			stop("the initial distribution of the diffusion process is not available in the stats package.")
		Q <- eval(parse(text = sub("d", "q", id[[1]])), list(x = c(10^-4, 1 - 10^-4)))
		X0 <- seq(Q[1], Q[2], length.out = m)
		df.X0 <- eval(expr, list(x = X0))
		cp <- FALSE		
	}
	if (!is.numeric(S) & !is.character(S)) stop(paste(sQuote("S"), "argument is not of numeric or character type."))
    	if (length(S) > 1) stop(paste(sQuote("S"), "argument is not of length 1."))
    	exprS <- parse(text = S)
      if (inherits(try(exprS, silent = TRUE), "try-error")) stop(paste("the mathematical expression of the boundary shows syntax errors.", sep = ""))
      if (inherits(try(D(exprS, "t"), silent = TRUE), "try-error")) stop("R can not compute the symbolic derivative with respect to 't' of the mathematical expression of the boundary")     
	if ((!is.list(env)) & (!is.null(env))) stop(paste(sQuote("env"), "argument is not NULL or a list object."))
    	if (!all(is.element(sapply(env, mode), c("numeric", "character")))) stop(paste(sQuote("env"), "argument is not a list of numeric or character objects."))
	if (!is.numeric(n)) stop(paste(sQuote("n"), "argument is not of numeric type."))
    	if (round(n) != n) stop(paste(sQuote("n"), "argument is not an integer number."))
    	if (length(n) > 1) stop(paste(sQuote("n"), "argument is not of length 1.")) 
	if (missing(it.max)){
		if (cp) it.max <- 50000L else it.max <- 1000000L
	}      		
    	
	Call <- match.call()
	Args <- list(t0 = t0, T = T, S = S, variableStep = variableStep, from.t0 = from.t0, to.T = to.T, r = r, zeroSlope = zeroSlope, p0.tol = p0.tol, k = k, 
			m = m, n = n, p = p, alpha = alpha, skip = skip, tol = tol, it.max = it.max)
	label <- intersect(names(Args), names(Call))
	Call[label] <- Args[label]
	if (is.name(Call$id)) Call$id <- id	

	e <- is.null(env)
	if (e){
		B <- function(t, x0) NULL
      	env3 <- list(t0 = t0)
      	exprS <- as.expression(eval(substitute(substitute(expr, env3), list(expr = exprS[[1]]))))
		if (is.element("t", all.names(expr))) body(B) <- as.call(c(as.name("{"), exprS))
		else body(B) <- substitute(rep(a, length(t)), list(a = exprS[[1]]))
	}
	else{
		env2 <- env	
      	logic <- unlist(lapply(env2, is.character))
      	if (any(logic))
			env2[logic] <- lapply(env2[logic], function(x) parse(text=x)[[1]]) 
		env3 <- c(t0 = t0, env2)         

		B <- function(t, x0) NULL
      	exprS <- as.expression(eval(substitute(substitute(expr, env3), list(expr = exprS[[1]]))))			
		if (is.element("t", all.names(exprS))){
			vars <- setdiff(intersect(all.names(exprS), names(env3)), c("t", "t0", "x0"))
      		if (length(vars) > 0) expr <- c(sapply(vars, function(x) as.expression(eval(substitute(substitute(expr, list(valor=env3[[x]])), list(expr = parse(text = paste(x, "<- valor"))[[1]]))))), exprS)
			else expr <- exprS
      		body(B) <- as.call(c(as.name("{"), expr))
		}
		else body(B) <- substitute(rep(a, length(t)), list(a = exprS[[1]]))		
	}

      rho <- sign(sapply(X0, B, t = t0) - X0)
	
	logic <- rho != 0L

	if (any(logic)){
		if (e){
			A1 <- function(x, t) NULL
			expr <- parse(text = dp$mean)
			if (any(is.element(c("x", "t"), all.names(expr)))) body(A1) <- expr
			else body(A1) <- substitute(rep(a, max(length(x),length(t))), list(a = expr[[1]]))

			DA2 <- function(x, t) NULL
			expr <- parse(text = dp$var)
      		if (any(is.element(c("x", "t"), all.names(expr)))) body(DA2) <- deriv(expr, "x")
			else body(DA2) <- eval(substitute(expression({
    							.value <- rep(a, length(x))
    							.grad <- array(0, c(length(.value), 1L), list(NULL, c("x")))
    							.grad[, "x"] <- 0
    							attr(.value, "gradient") <- .grad
    							.value
							}), list(a = expr[[1]])))			

			Df <- function(x, t, y, s) NULL
      		expr <- deriv(parse(text = dp$tpdf), "x")
      		body(Df) <- as.call(c(as.name("{"), expr))

			DB <- function(t, x0) NULL
      		if (is.element("t", all.names(exprS))) body(DB) <- deriv(exprS, "t")
			else body(DB) <- eval(substitute(expression({
    						.value <- rep(a, length(t))
    						.grad <- array(0, c(length(.value), 1L), list(NULL, c("t")))
    						.grad[, "t"] <- 0
    						attr(.value, "gradient") <- .grad
    						.value
						}), list(a = exprS[[1]])))			
		}
		else{
			A1 <- function(x, t) NULL
			expr <- as.expression(eval(substitute(substitute(expr, env2), list(expr = parse(text = dp$mean)[[1]]))))
			if (any(is.element(c("x", "t"), all.names(expr)))){
      			vars <- setdiff(intersect(all.names(expr), names(env2)), c("x", "t"))
      			if (length(vars) > 0) expr <- c(sapply(vars, function(x) as.expression(eval(substitute(substitute(expr, list(valor=env2[[x]])), list(expr = parse(text = paste(x, "<- valor"))[[1]]))))), expr)
      			body(A1) <- as.call(c(as.name("{"), expr))
			}
			else body(A1) <- substitute(rep(a, max(length(x),length(t))), list(a = expr[[1]]))

			DA2 <- function(x, t) NULL
			expr <- parse(text = dp$var)
			if (any(is.element(c("x", "t"), all.names(expr)))){
				exprD <- as.expression(eval(substitute(substitute(expr, env2), list(expr = deriv(expr, "x")[[1]]))))
      			vars <- setdiff(intersect(all.names(exprD), names(env2)), c("x", "t"))
      			if (length(vars) > 0) exprD <- c(sapply(vars, function(x) as.expression(eval(substitute(substitute(expr, list(valor=env2[[x]])), list(expr = parse(text = paste(x, "<- valor"))[[1]]))))), exprD)
      			body(DA2) <- as.call(c(as.name("{"), exprD))
			}
			else body(DA2) <- as.call(c(as.name("{"), eval(substitute(expression({
    							.value <- rep(a, length(x))
    							.grad <- array(0, c(length(.value), 1L), list(NULL, c("x")))
    							.grad[, "x"] <- 0
    							attr(.value, "gradient") <- .grad
    							.value
							}), list(a = eval(substitute(substitute(expr, env2), list(expr = expr[[1]]))))))))

			Df <- function(x, t, y, s) NULL
			expr <- parse(text = dp$tpdf)
      		exprD <- as.expression(eval(substitute(substitute(expr, env2), list(expr = deriv(expr, "x")[[1]]))))
      		vars <- setdiff(intersect(all.names(exprD), names(env2)), c("x", "t", "y", "s"))
      		if (length(vars) > 0) exprD <- c(sapply(vars, function(x) as.expression(eval(substitute(substitute(expr, list(valor=env2[[x]])), list(expr = parse(text = paste(x, "<- valor"))[[1]]))))), exprD)
      		body(Df) <- as.call(c(as.name("{"), exprD))      		

			DB <- function(t, x0) NULL
      		if (is.element("t", all.names(exprS))){
				exprD <- deriv(exprS, "t")
      			vars <- setdiff(intersect(all.names(exprD), names(env3)), c("t", "x0"))
      			if (length(vars) > 0) exprD <- c(sapply(vars, function(x) as.expression(eval(substitute(substitute(expr, list(valor=env3[[x]])), list(expr = parse(text = paste(x, "<- valor"))[[1]]))))), exprD)
      			body(DB) <- as.call(c(as.name("{"), exprD))
			}
			else body(DB) <- as.call(c(as.name("{"), eval(substitute(expression({
    							.value <- rep(a, length(t))
    							.grad <- array(0, c(length(.value), 1L), list(NULL, c("t")))
    							.grad[, "t"] <- 0
    							attr(.value, "gradient") <- .grad
    							.value
							}), list(a = exprS[[1]])))))			
		}	      				
		
		X0 <- X0[logic]
		rho <- rho[logic]
		m <- length(X0)
		d <- diff(X0)

    		SFPTL <- lapply(X0, function(x, DP, T0, T, B, Env, N, zS, tol, K) unclass(summary(FPTL(DP, T0, T, x, B, Env, N), zS, tol, K)), 
						DP = dp, T0 = t0, T = T, B = S, Env = env, N = r, zS = zeroSlope, tol = p0.tol, K = k)
		SFPTLCall <- sapply(SFPTL, attr, "Call")
		FPTLCall <- sapply(SFPTL, attr, "FPTLCall")
		dp <- attr(SFPTL[[1]], "dp")
		vars <- attr(SFPTL[[1]], "vars")
		SFPTL <- unlist(SFPTL, recursive=FALSE)
		attr(SFPTL, "Call") <- SFPTLCall
		attr(SFPTL, "FPTLCall") <- FPTLCall
		attr(SFPTL, "dp") <- dp
		if (length(vars) > 0L) attr(SFPTL, "vars") <- vars
		if (is.list(id)) attr(SFPTL, "id") <- id
		class(SFPTL) <- "summary.fptl"
    		IS <- Integration.Steps(SFPTL, variableStep, from.t0, to.T, n, p, alpha)
		H <- IS$H
		Skip <- IS$skip

		it <- (H[, 2] - H[, 1])/H[, 3]

      	if ((m*sum(it) > it.max) & interactive()) {
			M <- m - sapply(Skip, sum)
			label1 <- character(1)
			if (!to.T){
				logic <- sapply(Skip, function(s) all(!s))            		
				if (any(logic)){
					itStop <- m*cumsum(it)[logic]
					if (length(itStop) == 1) label1 <- paste("iteration", itStop) else label1 <- paste("iterations ", paste(itStop[-length(itStop)], collapse = ", "), " and ", itStop[length(itStop)], sep = "")
					label1 <- paste(", although the algorithm can stop at", label1)
					itStop <- cumsum(M*it)[logic]
					if (length(itStop) == 1) label2 <- paste("iteration", itStop) else label2 <- paste("iterations ", paste(itStop[-length(itStop)], collapse = ", "), " and ", itStop[length(itStop)], sep = "")
					label2 <- paste(", although the algorithm can stop at", label2)
            		}            		
        		}
        		itbyStep <- format(data.frame(matrix(H[, 1:2], ncol = 2), m*it), digits = 8)        		
			logic <- M < m
        		if ((skip) & (any(logic)))
            		itbyStep[, 3][logic] <- paste(M[logic]*it[logic], ifelse(cp, "or", "to"), itbyStep[, 3][logic])
        		row.names(itbyStep) <- paste("Subinterval", 1:nrow(itbyStep))
        		names(itbyStep) <- c("Lower end", "Upper end", "Iterations")
        		cat("\nThe table below shows the number of iterations of the approximation process step by step: \n\n")
        		print(itbyStep)
        		iterations <- m*sum(it)
			if ((skip) & (any(logic))){
				cat("\nIf no subinterval is avoided, the total number of iterations is ", iterations, label1, ".", sep = "")
				cat("\nIf the application of the numerical algorithm is avoided in all suitable subintervals, the total number of iterations is ", sum(M*it), label2, ".", sep = "") 
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

		jumps <- vector("list", nrow(H))
    		Integral <- numeric(nrow(H))
		now <- matrix(, nrow = nrow(H) + 1, ncol = 2)

		cat("\nComputing...\t")
    		now[1, ] <- proc.time()[1:2]

    		h <- numeric(0)
    		F <- 0
    		g0 <- H[1, 1]
		gx <- H[1, 1] + H[1, 3]
   			
		cols <- 1:m

		if (skip) logic <- !Skip[[1]] else logic <- rep(TRUE, m)
		if (any(logic)){
    			b <- lapply(X0, DB, t = gx)
			B2 <- b[logic]
			dB2 <- sapply(B2, attr, "gradient")
			b <- unlist(b) 

			a1 <- sapply(B2, A1, t = gx)
    			a2 <- lapply(B2, DA2, t = gx)
			da2 <- sapply(a2, attr, "gradient")
			a2 <- unlist(a2)

    			f0 <- mapply(Df, B2, X0[logic], MoreArgs = list(t = gx, s = t0), SIMPLIFY = FALSE)
			df0 <- sapply(f0, attr, "gradient")
			f0 <- unlist(f0)

			gy0 <- rep(0L, m)
    			gy0[logic] <- pmax(0, -rho[logic] * (f0 * (dB2 - a1 +  0.75 * da2) + df0 * a2))
			gy0[gy0 < 10^-10] <- 0L
			if (cp) gy <- gy0
			else{
				gw <- df.X0 * gy0
				gy <- sum(d * (gw[-m] + gw[-1]))/2
			}
			Integral[1] <- H[1, 3] * gy/2
			jumps[[1]] <- cols[!logic]						
		}
		else{				
            	b <- sapply(X0, DB, t = gx)
            	gy0 <- rep(0L, m)
			gy <- 0L
            	Integral[1] <- 0
            	jumps[[1]] <- cols	
		}
						
		if (gx < H[1, 2]) rows <- 1:nrow(H)
		else{
			rows <- 2:nrow(H)
			now[2, ] <- proc.time()[1:2]
		}

		g1n <- gy0[1]		
		gy0 <- matrix(gy0, nrow = 1)
		H[1, 1] <- gx			
			
		logic <- rep(TRUE, m)
		
		if (cp) expr <- expression(g1n <- tail(gy0[,1], length(u2)))
		else expr <- expression({gw <- df.X0 * t(tail(gy0, length(u2))); g1n <- (d %*% (gw[-m, ] + gw[-1, ]))/2})

    		for (i in rows) {
			if (skip) logic <- !(Skip[[i]] & (gy0[nrow(gy0), ] < 10^-10))
			if (any(logic)){					
            		u2 <- seq(H[i, 1] + H[i, 3], H[i, 2], by = H[i, 3])
            		h2 <- rep(H[i, 3], length(u2))
				b2 <- lapply(X0, DB, t = u2)
				B2 <- b2[logic]
				dB2 <- as.matrix(as.data.frame(lapply(B2, attr, "gradient")))
				a1 <- as.matrix(as.data.frame(lapply(B2, A1, t = u2)))
				a2 <- lapply(B2, DA2, t = u2)

				da2 <- as.matrix(as.data.frame(lapply(a2, attr, "gradient")))
				f0 <- mapply(Df, B2, X0[logic], MoreArgs = list(t = u2, s=t0), SIMPLIFY = FALSE)
				df0 <- as.matrix(as.data.frame(lapply(f0, attr, "gradient")))
				f0 <- as.matrix(as.data.frame(f0))
				b2 <- as.matrix(as.data.frame(b2))
				a2 <- as.matrix(as.data.frame(a2))
					
				a <- dB2 - a1 + 0.75 * da2
				G <- matrix(0L, nrow = length(u2), ncol = m)
				G[, logic] <- - f0 * a - df0 * a2
				gy0 <- rbind(gy0, G)

				n <- length(gx)
            		h <- c(h, h2)
            		gx <- c(gx, u2)
				b <- rbind(b, b2)

				for (k in (1:length(u2))) {
                			index <- 1:(n + k - 1)
					f1 <- mapply(Df, b2[k, logic], as.list(as.data.frame(b[index, logic, drop=FALSE])), MoreArgs = list(t = u2[k], s = gx[index]), SIMPLIFY = FALSE)
                			df1 <- as.matrix(as.data.frame(lapply(f1, attr, "gradient")))
					f1 <- as.matrix(as.data.frame(f1))
					gy0[n + k, logic] <- pmax(0, rho[logic] * (gy0[n + k, logic] + h[index] %*% (gy0[index, logic, drop=FALSE] * t(a[k, ] * t(f1) + a2[k, ] * t(df1)))))
				}
			
				eval(expr)
				Integral[i] <- Integral[i] + h2[1] * ((gy[length(gy)] - g1n[length(g1n)])/2 + sum(g1n))
				F <- F + Integral[i]
				jumps[[i]] <- cols[!logic]
				gy <- c(gy, g1n)							
				now[i+1, ] <- proc.time()[1:2]									
			}
			else{
				correction <- h[length(h)] * gy[length(gy)]/2
            		Integral[i - 1] <- Integral[i - 1] - correction
            		F <- F - correction
            		gy0[nrow(gy0), ] <- 0L
            		h <- c(h, H[i, 2] - H[i, 1])
            		gx <- c(gx, H[i, 2])
            		b <- rbind(b, sapply(X0, DB, t = H[i, 2]))
            		gy0 <- rbind(gy0, rep(0L, m))
				gy <- c(gy, 0L)
            		now[i + 1, ] <- proc.time()[1:2]
            		Integral[i] <- 0
            		jumps[[i]] <- cols					
        		}
			if ((!to.T) & all(!Skip[[i]]) & (F >= (1 - tol))) break       				 			
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

		gx <- c(g0, gx)
		gy <- c(0L, gy)

		if (cp) out <- list(x = gx, y = gy, y.x0 = NULL)
		else{
			gy0 <- rbind(rep(0L, m), gy0)
			dimnames(gy0) <- NULL
			out <- list(x = gx, y = gy, y.x0 = gy0)
		}
		if (i < nrow(H)) {
        		Integral <- Integral[1:i]
        		jumps <- jumps[1:i]
        		now <- now[1:(i + 1), ]
    		}

    		CPUTime <- apply(now, 2, diff)

    		if (!is.matrix(CPUTime)) CPUTime <- matrix(CPUTime, ncol = 2)

		attr(out, "Call") <- Call
		attr(out, "Steps") <- IS$H						
    		attr(out, "cumIntegral") <- cumsum(Integral)
    		attr(out, "skips") <- jumps
    		attr(out, "CPUTime") <- CPUTime
		attr(out, "summary.fptl") <- SFPTL
    		class(out) <- c("fpt.density", "list")		
	}
	else stop("S(t0) = x0 for all x0 in X0")

	return(out)      
}
