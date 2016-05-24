report.fpt.density <-
function (obj, report.sfptl = FALSE, tex = FALSE, digits = 8, 
    ...) 
{
	if (!is.fpt.density(obj)) 
        stop(paste(sQuote("obj"), "is not of class", shQuote("fpt.density"))) 

	cp <- is.null(obj$y.x0)
    	A <- attr(obj, "Call")
	f <- A[[1]]
	Args <- formals(eval(f))
	logic <- (as.character(f) == "Approx.fpt.density")
    	Args[names(A)[-1]] <- A[-1]

	ArgsFPTL <- formals(FPTL)
	A <- attr(attr(obj, "summary.fptl"), "FPTLCall")[[1]][-1]
	ArgsFPTL[names(A)] <- A
	if (!cp) ArgsFPTL$x0 <- sapply(attr(attr(obj, "summary.fptl"), "FPTLCall"), "[[", "x0")
		
	dp <- attr(attr(obj, "summary.fptl"), "dp")
	S <- ArgsFPTL$S
	env <- ArgsFPTL$env
	if (is.name(env)) env <- attr(attr(obj, "summary.fptl"), "vars")[[as.character(env)]] else env <- as.list(env)[-1]

	if (!is.null(env)){
		if (cp) env <- c(t0 = ArgsFPTL$t0, x0 = ArgsFPTL$x0, env) else env <- c(t0 = ArgsFPTL$t0, env)
		A1 <- parse(text = dp$mean)[[1]]
		A2 <- parse(text = dp$var)[[1]]
		if (is.character(S)) S <- parse(text = S)[[1]]								
		l1 <- (sapply(env, length) == 1L)
		if (any(l1)){
			env1 <- env[l1]
			logic <- sapply(env1, is.character)
			if (any(logic)){
				env2 <- lapply(env1[logic], function(x) eval(parse(text = paste("substitute(", x, ")", sep = ""))))
				A1 <- eval(substitute(substitute(expr, env2), list(expr = A1)))
	      		A2 <- eval(substitute(substitute(expr, env2), list(expr = A2)))
				if (!is.numeric(S)) S <- eval(substitute(substitute(expr, env2), list(expr = S)))
			}
			logic <- sapply(env1, is.numeric)
			if (any(logic)){
				env2 <- env1[logic]
				A1 <- eval(substitute(substitute(expr, env2), list(expr = A1)))
	      		A2 <- eval(substitute(substitute(expr, env2), list(expr = A2)))
				if (!is.numeric(S)) S <- eval(substitute(substitute(expr, env2), list(expr = S)))
			}
		}
		l2 <- !l1
		if (!is.null(attr(attr(obj, "summary.fptl"), "vars")) | any(l2)){
			env2 <- c(attr(attr(obj, "summary.fptl"), "vars"), env[l2])
			v <- intersect(union(all.names(A1), all.names(A2)), names(env2))
			if (length(v) > 0L){
				w <- lapply(env2[v], function(x) paste("(", paste(x, collapse = ", "), "),", sep = "")) 
				w[[length(w)-1L]] <- substring(w[[length(w)-1L]], 1L, nchar(w[[length(w)-1L]])-1L)			
			}			
		}
		else v <- character(0)
		
		A1 <- as.character(parse(text = deparse(A1)))
		A2 <- as.character(parse(text = deparse(A2)))
		if (!is.numeric(S)) S <- as.character(parse(text = deparse(S)))
	}
	else{
		A1 <- as.character(parse(text = dp$mean))
		A2 <- as.character(parse(text = dp$var))		
		if (!is.numeric(S)){
			if (cp) env2 <- list(t0 = ArgsFPTL$t0, x0 = ArgsFPTL$x0) else env2 <- list(t0 = ArgsFPTL$t0)	
			S <- as.character(parse(text = deparse(eval(substitute(substitute(expr, env2), list(expr = parse(text = S)[[1]]))))))
		}
		v <- character(0)
	}
	
	if (tex) { 
		noindent <- "\\noindent "
        	vskip <- "\n\\vskip 10pt "
        	dollar <- "$"
		sim <- " \\sim "
        	x0.label <- "x_{0}"
        	ti.label <- "t_{i}^{\\phantom{*}}"
        	labels <- c("$I_{i}^{\\phantom{+}}$", "$t_{i}^{*}$", "$t_{max,i}^{-}$", "$t_{max,i}^{\\phantom{+}}$", "$t_{max,i}^{+}$")
		alpha <- "$\\alpha$"
		Q <- "$q_{\\alpha}$"
		IQ <- "[$q_{\\mbox{\\tiny $10^{-4}$}}$, $q_{\\mbox{\\tiny $1-10^{-4}$}}$]"      	
    	}
    	else {
        	noindent <- ""
        	vskip <- "\n"
        	dollar <- ""
		sim <- " ~ "
        	x0.label <- "x0"
        	ti.label <- "t[i]"
        	labels <- c("I[i]", "t[i]*", "tmax[i]^-", "tmax[i]", "tmax[i]^+")
		alpha <- "alpha"
		Q <- "q[alpha]"
		IQ <- "[q[10^(-4)], q[1-10^(-4)]]"		  	
	}
	
	aux <- function(paragraph, m){
			N <- nchar(paragraph)
			if (m < N){
				index <- unlist(gregexpr(" ", paragraph))
				J <- numeric(0)
				i <- m
				while(i < N){
					i <- index[tail(which(index <= i),1)]
					J <- c(J, i)
					i <- m + i
				}
				return(substring(paragraph, c(1, J+1L), c(J, N)))
			}
			else return(paragraph)
		}

	pardp <- "stores an approximation of the first-passage-time density through the boundary"
	if (is.name(match.call()$obj)) pardp <- paste("The fpt.density class object", ifelse(tex, paste("{\\tt ", match.call()$obj, "}", sep=""), 
									shQuote(match.call()$obj)), pardp)
	else pardp <- pardp <- paste("This fpt.density class object", pardp)		
	cat(aux(pardp, 125L), sep="\n")
	if (tex){
		S <- paste("\\verb#", aux(S, 60L), "#", sep = "")
		if (length(S) > 1L){
			cat("$$\\arraycolsep 2pt \\begin{array}{ll}")
			cat(paste(c("\n\\mbox{\\tt S(t)} =", character(length(S)-1)), S, sep = " & "), sep = " \\\\ \n")
			cat("\\end{array}$$")
		}
		else cat("$$\\mbox{\\tt S(t)} = ", S, "$$", sep = "")
		cat(aux(paste("\nof the diffusion process $\\{X(t) \\thinspace ; \\thinspace ", ArgsFPTL$t0, " \\leq t \\leq ", ArgsFPTL$T, 
				" \\}$ with infinitesimal moments \\medskip", sep = ""), 125L), sep = "\n")
		A1 <- paste("\\verb#", aux(A1, 60L), "#", sep = "")
		if (length(A1) > 1L){
				cat("\n\\qquad $\\arraycolsep 2pt \\begin{array}{ll}")
				cat(paste(c("\n\\mbox{\\tt A}_1\\mbox{\\tt (x,t)} =", character(length(A1)-1)), A1, sep = " & "), sep = " \\\\ \n")
				cat("\\end{array}$ \\medskip")
		}
		else cat("\n\\qquad $\\mbox{\\tt A}_1\\mbox{\\tt (x,t)} = ", A1, "$ \\medskip", sep = "")
		cat("\n\n\\noindent and \\medskip \n")
		A2 <- paste("\\verb#", aux(A2, 60L), "#", sep = "")
		if (length(A2) > 1L){
				cat("\n\\qquad $\\arraycolsep 2pt \\begin{array}{ll}")
				cat(paste(c("\n\\mbox{\\tt A}_2\\mbox{\\tt (x,t)} = ", character(length(A2)-1)), A2, sep = " & "), sep = " \\\\ \n")
				cat("\\end{array}$ \\medskip")
		}
		else cat("\n\\qquad $\\mbox{\\tt A}_2\\mbox{\\tt (x,t)} = ", A2, "$ \\medskip", sep = "")
		cat("\n")
		if (length(v) > 0L){
			cat("\n\\noindent where \\medskip \n")
			v <- paste(v, "=")
			w <- lapply(w, aux, m = 60L)
			logic <- (sapply(w, length) > 1L)
			w[logic] <- mapply(function(x, y) c("\n\\qquad {\\tt \\tabcolsep 2pt \\begin{tabular}{ll}", paste(c(x, character(length(y)-1L)), 
						" & ", y, "\\\\", sep = ""), "\\end{tabular}} \\medskip"), v[logic], w[logic], SIMPLIFY = FALSE)
			w[!logic] <- mapply(function(x, y) paste("\n\\qquad {\\tt ", x, y, "} \\medskip", sep = ""), v[!logic], w[!logic], SIMPLIFY = FALSE)
			lapply(w[-length(w)], cat, sep = "\n")
			cat("\n\\noindent and\n")
			cat(w[[length(w)]], sep = "\n")	
		}
		cat("\n")
	}
	else{
		n <- nchar("S(t) = ")
		S <- aux(S, 125L - 8L - n) 
		S <- paste("\t", c("S(t) = ", rep(paste(rep(" ", n), collapse = ""), length(S)-1L)), S, sep = "")
		cat(S, sep = "\n")
		cat(aux(paste("of the diffusion process {X(t); ", ArgsFPTL$t0, " <= t <= ", ArgsFPTL$T, "} with infinitesimal moments", sep = ""), 125L), sep = "\n")
		n <- nchar("A1(x,t) = ")
		A1 <- aux(A1, 125L - 8L - n) 
		A1 <- paste("\t", c("A1(x,t) = ", rep(paste(rep(" ", n), collapse = ""), length(A1)-1L)), A1, sep = "")
		cat(A1, sep = "\n")
		cat("and\n")
		n <- nchar("A2(x,t) = ")
		A2 <- aux(A2, 125L - 8L - n) 
		A2 <- paste("\t", c("A2(x,t) = ", rep(paste(rep(" ", n), collapse = ""), length(A2)-1L)), A2, sep = "")
		A2[length(A2)] <- paste(A2[length(A2)], ",", sep = "")
		cat(A2, sep = "\n")
		if (length(v) > 0L){
			cat("where\n")
			v <- paste(v, "= ")
			w <- lapply(w, aux, m = 125L - 8L - max(nchar(v)))
			w <- mapply(function(x,y) paste("\t", c(x, rep(paste(rep(" ", nchar(x)+1L), collapse = ""), length(y)-1L)), y, sep = ""), v, w)
			lapply(w[-length(w)], cat, sep = "\n")
			cat("and\n")
			cat(w[[length(w)]], sep = "\n")	
		}		
	}
		
	if (cp) cat(aux(paste(noindent, "conditioned to ", dollar, "X(", ArgsFPTL$t0, ") = ", ArgsFPTL$x0, dollar, ".", sep = ""), 125L), sep = "\n")
	else{
		cat(aux(paste(noindent, "and initial distribution ", dollar, "X(", ArgsFPTL$t0, ")", sim, ifelse(tex, attr(attr(obj, "summary.fptl"), "id")[[3]], 
			attr(attr(obj, "summary.fptl"), "id")[[4]]), dollar, ".", sep = ""), 125L), sep = "\n")
		
		cat(aux(paste(vskip, noindent, "The first-passage-time density function has been obtained by numerical integration, in the range of variation of ", 
			dollar, "X(", ArgsFPTL$t0, ")", dollar, ", of the corresponding first-passage-time densities conditioned to values of ", dollar, 
			"X(", ArgsFPTL$t0, ")", dollar, ", weighted by the initial density function. Concretely, ", Args$m, " equally spaced values of ", 
			dollar, "X(", ArgsFPTL$t0, ")", dollar, " in ", IQ, " has been considered, where ", Q, " is the ", alpha, "-quantile for the distribution of ", 
			dollar, "X(", ArgsFPTL$t0, ")", dollar, ".", sep = ""), 125L), sep = "\n")
	}

	cat(aux(paste(vskip, noindent, "The approximation process makes use of the First-Passage-Time Location (FPTL) function to locate the first-passage-time variable.", 
		sep = ""), 125L), sep = "\n")

    	if (report.sfptl) report(attr(obj, "summary.fptl"), tex, digits, title="", heading=FALSE)

	to.T <- Args$to.T

	m <- length(ArgsFPTL$x0)
	if (Args$skip) jumps <- which(sapply(attr(obj, "skips"), identical, 1:m)) else jump <- integer(0)
	
	nI <- length(attr(obj, "cumIntegral"))
    	index <- 1:nI
    	y <- data.frame(matrix(, nrow = nI, ncol = 5))
    	names(y) <- c("Subinterval", "Integration step", "Cumulative integral", "Iterations", "User time")
    
    	lower <- attr(obj, "Steps")[index, 1]
    	upper <- attr(obj, "Steps")[index, 2]

	y[, 2] <- attr(obj, "Steps")[index, 3]
    	y[, 3] <- attr(obj, "cumIntegral")
    	y[, 4] <- (upper - lower)/y[, 2]
    	y[, 5] <- attr(obj, "CPUTime")[,1]

	it <- sum(y[setdiff(index, jumps), 4])
	ut <- sum(y[, 5])
	y <- format(y, digits = digits, ...)
    	y[,1] <- paste("(", format(lower, digits = digits, ...), ", ", format(upper, digits = digits, ...), "]", sep="")

	paragraph <- paste(vskip, noindent, "The f.p.t. density has been approximated from ", ArgsFPTL$t0, " to ", ArgsFPTL$T, " using a ", ifelse(Args$variableStep, "variable", "fixed"), 
		" integration step", ifelse(Args$skip, " and avoiding the application of the numerical algorithm in those subintervals in which this is possible", 
		""), ".", sep = "")

	j <- length(jumps)
	if (j > 0L){		
		paragraph <- paste(paragraph, " In particular, in this case it has avoided the application of the numerical algorithm on the subinterval", rep("s:", j > 1L), " ", paste(y[jumps, 1], c(rep(", ", max(0L, j - 2)), 
					rep(" and ", j > 1L), character(1)), sep = "", collapse = ""), ".", sep = "")
		y[jumps, 4] <- "0"
	}
	else paragraph <- paste(paragraph, " In particular, in this case no interval has been avoided in the application of the numerical algorithm.", sep = "")

	cat(aux(paragraph, 125L), sep = "\n")

    	cumI <- attr(obj, "cumIntegral")[nI]
	tstop <- paste(dollar, "t", dollar, " = ", format(upper[nI], digits = digits, ...), sep = "")
    		 
	if (to.T) paragraph <- paste(vskip, noindent, "The value of the cumulative integral of the approximation is ", format(cumI, digits = digits), ".", sep = "")
	else {        
        	if ((cumI > Args$tol) & (nI < nrow(attr(obj, "Steps")))) 
            	paragraph <- paste(vskip, noindent, "The algorithm was stopped at ", tstop, ", since the value of the cumulative integral of the approximation is ", format(cumI, digits = digits), 
				ifelse(tex, " $\\geq$ ", " >= "), "1 - tol.", sep = "")
        	else paragraph <- paste(vskip, noindent, "The algorithm was stopped at ", tstop, " and the value of the cumulative integral of the approximation is ", format(cumI, digits = digits), ".", sep = "")
    	}
    	paragraph <- paste(paragraph, " The total number of iterations is ", format(it,...), " and the user time employed was ", format(ut,...), " (in seconds).", sep = "")

	cat(aux(paragraph, 125L), sep = "\n")
    	
	cat(vskip, noindent, "The table below shows the approximation process step by step:", sep = "")
	labels <- c("Subintervals", "Integration steps", "Cumulative integral", "Iterations", "User time")

	if (tex){
		y[,1] <- paste("{", y[,1], "}", sep = "")
		labels <- paste("\\fbox{", labels, "}", sep = "")
		cat("\n\\fboxrule 0pt \\fboxsep 2pt")
		cat("\n\\begin{longtable}{|", rep("r|", ncol(y)), "}", sep = "")
		cat("\n\\hline ")
		endhead <- paste("\\multicolumn{1}{", c("|", character(ncol(y)-1L)), ">{\\columncolor[gray]{0.7}}c|}{", labels, "}", sep = "")
		cat(endhead, sep = " & \n")
		cat("\\\\ \\hline \\endfirsthead \\multicolumn{", ncol(y), "}{c}{\\tiny Continued from previous page} \\\\ \\hline ", "\n", sep = "")
		cat(endhead,  sep = " & \n")
		cat("\\endhead \\multicolumn{", ncol(y), "}{c}{\\tiny Continued on next page} \\endfoot \\endlastfoot", sep = "")
		cat(paste("\n", apply(y, 1, paste, collapse = " & "), " \\\\", sep = ""))
		cat(" \\hline")
		cat("\n\\end{longtable}")
		cat("\n")		
	}
	else{			
		y <- rbind(labels, y)
		y <- apply(y, 2, format, digits = digits, justify = "right")
		cat("\n\n")
		cat(apply(y, 1, paste, collapse = "  "), sep = "\n")						
	}	
}
