report.summary.fptl <-
function (obj, tex = FALSE, digits = 8, heading = TRUE, ...) 
{
	if (!is.summary.fptl(obj)) 
      	stop(paste(sQuote("obj"), "is not of class", shQuote("summary.fptl")))

  	cp <- (length(obj) == 1L)
	Args <- formals(summary.fptl)
	ArgsFPTL <- formals(FPTL)
	A <- attr(obj, "Call")[[1]][-1]
	Args[names(A)] <- A
	A <- attr(obj, "FPTLCall")[[1]][-1]
	ArgsFPTL[names(A)] <- A
	if (!cp) ArgsFPTL$x0 <- sapply(attr(obj, "FPTLCall"), "[[", "x0")
	
	dp <- attr(obj, "dp")
	S <- ArgsFPTL$S		
	env <- ArgsFPTL$env
	if (is.name(env)) env <- attr(obj, "vars")[[as.character(env)]] else env <- as.list(env)[-1]

	if (!is.null(env)){
		A1 <- parse(text = dp$mean)[[1]]
		A2 <- parse(text = dp$var)[[1]]		
		if (is.character(S)) S <- parse(text = S)[[1]]
		if (cp) env <- c(t0 = ArgsFPTL$t0, x0 = ArgsFPTL$x0, env) else env <- c(t0 = ArgsFPTL$t0, env)		
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
		if (!is.null(attr(obj, "vars")) | any(l2)){
			env2 <- c(attr(obj, "vars"), env[l2])
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
    	}
    	else {
        	noindent <- ""
        	vskip <- "\n"
        	dollar <- ""
		sim <- " ~ "
        	x0.label <- "x0"
        	ti.label <- "t[i]"
        	labels <- c("I[i]", "t[i]*", "tmax[i]^-", "tmax[i]", "tmax[i]^+")	  	
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

	if (heading){
		pardp <- "stores the information provided by the First-Passage-Time Location (FPTL) function about the location of the variation range of the first-passage-time variable through the boundary"
		if (is.name(match.call()$obj)) pardp <- paste("The summary.fptl class object", ifelse(tex, paste("{\\tt ", match.call()$obj, "}", sep=""), shQuote(match.call()$obj)), pardp)
	  	else pardp <- pardp <- paste("This summary.fptl class object", pardp)		
		cat(aux(pardp, 125L), sep="\n")
		if (tex){
			S <- paste("\\verb#", aux(S, 60L), "#", sep = "")
			if (length(S) > 1L){
				cat("$$\\arraycolsep 2pt \\begin{array}{ll}")
				cat(paste(c("\n\\mbox{\\tt S(t)} =", character(length(S)-1)), S, sep = " & "), sep = " \\\\ \n")
				cat("\\end{array}$$")
			}
			else cat("$$\\mbox{\\tt S(t)} = ", S, "$$", sep = "")
			cat(aux(paste("\nof the diffusion process $\\{X(t) \\thinspace ; \\thinspace ", ArgsFPTL$t0, " \\leq t \\leq ", ArgsFPTL$T, " \\}$ with infinitesimal moments \\medskip", sep = ""), 125L), sep = "\n")
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
		else cat(aux(paste(noindent, "conditioned to ", dollar, "X(", ArgsFPTL$t0, ") = ", x0.label, dollar, " for ", length(obj), 
			   " equally spaced values in the range of the initial distribution ", dollar, "X(", ArgsFPTL$t0, ")", sim, attr(obj, "id")[[3]], dollar, ".", sep = ""), 125L), sep = "\n")		
	}
	
	cat(aux(paste(vskip, noindent, "The FPTL function was evaluated at ", ArgsFPTL$n, " points in [", format(ArgsFPTL$t0, digits=digits), ", ", 
			format(ArgsFPTL$T, digits=digits), "] and it was considered constant in those growth subintervals where the slope of the function between the endpoints is less than ", 
			dollar, Args$zeroSlope, dollar, " degrees.", sep = ""), 125L), sep = "\n")

	cat(aux(paste(vskip, noindent, "In order to determine the time instants ", labels[2], " it has been considered that the value of the FPTL function is significantly bigger than the value at ", dollar, ti.label, dollar, 
		" if its difference is over ", dollar, "10^{", - Args$p0.tol, "}", dollar, " times the increment of the function between ", 
		dollar, ti.label, dollar, " and ", labels[4], ".", sep = ""), 125L), sep = "\n")
	
	x0 <- format(ArgsFPTL$x0, digits = digits)
	X <- matrix(numeric(0), ncol = 5)
	upper <- numeric(0)
	for (i in 1:length(obj)){
		upper <- c(upper, obj[[i]]$instants[,1][-1], ArgsFPTL$T)
		X <- rbind(X, obj[[i]]$instants)			
	}
	X <- apply(X, 2, format, digits = digits, ...)
	if (!is.matrix(X)) X <- matrix(X, ncol = 5)
	X[,1] <- paste("[", X[,1], ", ", format(upper, digits = digits), "]", sep = "")
	if (!cp){
		X <- cbind(paste(rep(" ", nchar(x0[1])), collapse = ""), X)
		m <- cumsum(sapply(obj, function(x) nrow(x$instants)))[-length(obj)]
		X[c(1, m + 1L), 1] <- x0		
	}

	if (tex){
		if (cp){
			item <- " "
			X[,1] <- paste("{", X[,1], "}", sep = "")
		}
		else{
			x0.label <- paste("$", x0.label, "^{\\phantom{*}}$", sep = "")
			item <- paste(", for each ", x0.label, ", ", sep = "") 
			labels <- c(x0.label, labels)
			X[,2] <- paste("{", X[,2], "}", sep = "")
		}
		cat(aux(paste(vskip, noindent, "The table below shows", item, "the interesting time instants in the subintervals determined by the growth points of the FPTL function and T = ", 
			format(ArgsFPTL$T, digits = digits), ":", sep = ""), 125L), sep = "\n")		 
		labels <- paste("\\fbox{", labels, "}", sep = "")
		cat("\n\\fboxrule 0pt \\fboxsep 2pt")
		cat("\n\\begin{longtable}{|", rep("r|", ncol(X)), "}", sep = "")
		cat("\n\\hline ")
		endhead <- paste("\\multicolumn{1}{", c("|", character(ncol(X)-1L)), ">{\\columncolor[gray]{0.7}}c|}{", labels, "}", sep = "")
		cat(endhead, sep = " & \n")
		cat("\\\\ \\hline \\endfirsthead \\multicolumn{", ncol(X), "}{c}{\\tiny Continued from previous page} \\\\ \\hline ", "\n", sep = "")
		cat(endhead,  sep = " & \n")
		cat("\\endhead \\multicolumn{", ncol(X), "}{c}{\\tiny Continued on next page} \\endfoot \\endlastfoot", sep = "")
		endline <- rep(" \\\\", nrow(X))
		if (!cp) endline[m] <- paste(endline[m], "\\hline")	
		cat(paste("\n", apply(X, 1, paste, collapse = " & "), endline, sep = ""))
		cat(" \\hline")
		cat("\n\\end{longtable}")
		cat("\n")		
	}
	else{	
		cat("\n")		
		if (cp) item <- " "
		else{
			item <- ", for each x0, " 
			labels <- c(x0.label, labels)
		}
		cat(aux(paste("The table below shows", item, "the interesting time instants in the subintervals determined by the growth points of the FPTL function and T = ", 
			format(ArgsFPTL$T, digits = digits), ":", sep = ""), 125L), sep = "\n")		 
		cat("\n")
		X <- rbind(labels, X)
		X <- apply(X, 2, format, justify = "right")
		cat(apply(X, 1, paste, collapse = "  "), sep = "\n")						
	}	
}
