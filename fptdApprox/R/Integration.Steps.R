Integration.Steps <-
function (sfptl, variableStep = TRUE, from.t0 = FALSE, to.T = FALSE, n = 250, p = 0.2, alpha = 1)
{
	Call <- match.call()
      if (!is.summary.fptl(sfptl)) stop(paste(sQuote("sfptl"), " object is not of class ", shQuote("summary.fptl"), ".", sep = ""))

      if (n <= 0) stop("n <= 0.")
     
      if (variableStep) {
            if ((p <= 0) | (p > 1)) stop("p is not in (0,1].")
      	if ((alpha <= 0) | (alpha > 1)) stop("alpha is not in (0,1].")
        	cotes <- c(n = 50, p = 0.1, alpha = 0.75)
        	logic <- (c(n, p, alpha) < cotes)
        	if (any(logic)) {
            	text <- paste(names(cotes)[logic], cotes[logic], sep = " >= ")
            	if (length(text) > 1)
                		text <- paste(paste(text[-length(text)], collapse = ", "), text[length(text)], sep = " and ")
            	text <- paste(text, "is recommended. If not, some integration steps can be too large.")
            	cat(text, "\n")
            	repeat {
                		answer <- readline("Do you want to continue? (y/n) ")
                		if ((answer == "n") | (answer == "y")) break
            	}
            	if (answer == "n") stop("\n\n")
            	else cat("\n")
        	}
      }
      else{
        	if (missing(p)) {
            	if (!missing(alpha)) cat("alpha argument is ignored.\n")
        	}
        	else {
            	if (missing(alpha)) cat("p argument is ignored.\n")
            	else cat("p and alpha arguments are ignored.\n")
        	}
      }

	if (length(sfptl) == 1L){
		t0 <- attr(sfptl, "FPTLCall")[[1]]$t0
      	T <- attr(sfptl, "FPTLCall")[[1]]$T
		SFPTL <- sfptl[[1]]$instants
		e <- c(t(SFPTL[, c(2,5)]))
		m <- 2 * nrow(SFPTL) - 1 
	    	skip <- rep(c(FALSE, TRUE), len = m)
		A <- diff(e)
		a <- rep(SFPTL[, 3] - SFPTL[, 2], each = 2, len = m)
		q <- A/a		
		if (any(skip)) q[skip] <- p*ifelse(q[skip] > 1, q[skip]^alpha, q[skip])
		H <- matrix( , nrow = m, ncol = 3)
	    	H[,1] <- e[-length(e)]
	    	H[,2] <- e[-1]	    	
	    	H[,3] <- A/ceiling(n*q)
		if (any(skip)){
			index <- which(skip)
			logic <- H[index, 3] < H[index - 1L, 3]
			if (any(logic)){
				J1 <- index[logic]
				J2 <- J1 - 1L
				H[J2, 2] <- H[J1, 2]
				w <- H[J2, 2] - H[J2, 1]
				H[J2, 3] <- w/ceiling(w/H[J2, 3])
				H <- H[-J1, ]
				e <- e[-J1]
				skip <- skip[-J1]
			}
		} 
		if (from.t0 & (t0 < e[1])){
			A <- e[1] - t0
			q <- A/a[1]
			h <- A/ceiling(n*p*ifelse(q > 1, q^alpha, q))
			if (h >= H[1,3]){
				H <- rbind(c(t0, e[1], h), H)
				skip <- c(FALSE, skip)
				e <- c(t0, e)
			}
			else{
				e[1] <- t0
				H[1, 1] <- t0
				w <- e[2] - t0
				H[1, 3] <- w/ceiling(n*w/a[1])
			} 
		}
		if (to.T & (e[length(e)] < T)){
			A <- T - e[length(e)]
			q <- A/a[length(a)]
			h <-  A/ceiling(n*p*ifelse(q > 1, q^alpha, q))
			if (h >= H[nrow(H), 3]){
				H <- rbind(H, c(e[length(e)], T, h))
				skip <- c(skip, FALSE)
				e <- c(e, T)
			}
			else{
				e[length(e)] <- T
				H[nrow(H), 2] <- T
				w <- T - e[length(e)-1L]
				H[1, 3] <- w/ceiling(n*w/a[length(a)])
			}
		}	
		skip <- as.list(skip)
      }
      else{		
		t0 <- attr(sfptl, "FPTLCall")[[1]]$t0
		T <- attr(sfptl, "FPTLCall")[[1]]$T
		SFPTL <- lapply(sfptl, "[[", "instants")
		m <- sapply(SFPTL, nrow)
		S <- lapply(2*m - 1, rep_len, x = c(FALSE, TRUE))
		E <- lapply(SFPTL, function(x) c(t(x[,c(2,5)])))
		a <- lapply(SFPTL, function(x) rep(x[,3] - x[,2], each = 2, len = 2*nrow(x)-1))
		e <- sort(unique(unlist(E)))
		from.t0 <- from.t0 & (t0 < e[1])
		to.T <- to.T & (e[length(e)] < T)
		if (from.t0) e <- c(t0, e)
		if (to.T) e <- c(e, T)
		l1 <- (sapply(E, head, n = 1L) > t0)
		l2 <- (sapply(E, tail, n = 1L) < T) 
		if (any(l1)){
			S[l1] <- lapply(S[l1], function(x) c(TRUE, x))
			E[l1] <- lapply(E[l1], function(x, y) c(y, x), y = t0)
			a[l1] <- lapply(a[l1], function(x) c(x[1], x))
		}
		if (any(l2)){
			S[l2] <- lapply(S[l2], function(x) c(x, TRUE))
			E[l2] <- lapply(E[l2], function(x, y) c(x, y), y = T)
			a[l2] <- lapply(a[l2], function(x) c(x, x[length(x)]))
		}
		A <- lapply(E, diff)
		Q <- mapply("/", A, a, SIMPLIFY = FALSE)
		logic <- sapply(S, any)
		if (any(logic)) 
			Q[logic] <- mapply(function(q, s, P, Alpha){q[s] <- P*ifelse(q[s]>1, q[s]^Alpha, q[s]); q}, Q[logic], S[logic], MoreArgs = list(P=p, Alpha=alpha), SIMPLIFY = FALSE)
		Q <- lapply(lapply(Q, "*", n), ceiling)
		P <- mapply("/", A, Q, SIMPLIFY = FALSE)
		J <- mapply(function(h, s) if (any(s)){i <- setdiff(which(s), 1L); i[h[i] < h[i-1]]} else integer(0), P, S, SIMPLIFY = FALSE)
		logic <- (sapply(J, length) > 0L)
		if (any(logic)){
			J <- J[logic]
			P[logic] <- mapply(function(e, h, b, j, N){k <- j-1L; w <- e[j+1L] - e[k]; h[k] <- w/ceiling(N*w/b[j]); h <- h[-j]; h}, E[logic], P[logic], a[logic], J, MoreArgs = list(N = n), SIMPLIFY = FALSE)
			E[logic] <- mapply(function(e, j) e[-j], E[logic], J, SIMPLIFY = FALSE)			
			S[logic] <- mapply(function(s, j) s[-j], S[logic], J, SIMPLIFY = FALSE)
			e <- intersect(e, sort(unique(unlist(E))))
		}
		if (any(l1)){
			logic <- sapply(P[l1], function(h) h[1] < h[2])
			if (any(logic)){
				J <- which(l1)[logic]
				P[J] <- mapply(function(e, h, b, N){w <- e[3] - e[1]; h[1] <- w/ceiling(N*w/b[1]); h <- h[-2]; h}, E[J], P[J], a[J], MoreArgs = list(N = n), SIMPLIFY = FALSE)
				E[J] <- lapply(E[J], function(e) e[-2])			
				S[J] <- lapply(S[J], function(s) s[-2])	
				e <- intersect(e, sort(unique(unlist(E))))
			}
		}		
		H <- matrix(, nrow = length(e) - 1, ncol = 3)
		H[,1] <- e[-length(e)]
		H[,2] <- e[-1]
		index <- lapply(E, findInterval, x=H[,1])
		skip <- as.list(as.data.frame(t(mapply("[", S, index))))
		P <- mapply("[", P, index)
		H[,3] <- apply(P, 1, min, na.rm = TRUE)		
		if (from.t0) skip[[1]] <- rep(FALSE, length(SFPTL))
		if (to.T) skip[[length(skip)]] <- rep(FALSE, length(SFPTL))
		d <- rev(which(c(0, diff(H[,3])) == 0)[-1])
		for (i in d){
			if (all(skip[[i-1]] == skip[[i]])){
				H[i-1,2] <- H[i,2]
				H <- H[-i,]
				skip <- skip[-i]
			}
		}
	} 

	w <- H[,2]-H[,1]
	H[,3] <- w/ceiling(w/H[,3])

	if (!variableStep){		
		l1 <- sapply(skip, all)
		l2 <- !l1
		h <- min(H[l2, 3])
		w <- H[nrow(H),2] - H[1,1]
		h <- w/ceiling(w/h)
		upper <- e[-1]
		upper[l2] <- e[1] + h * ceiling((upper[l2] - e[1])/h)
		upper[l1] <- e[1] + h * floor((upper[l1] - e[1])/h)
		lower <- c(e[1], upper[-length(upper)])		
		H[,1] <- lower
	      H[,2] <- upper
		H[ ,3] <- h
		l1 <- (upper == lower)
		if (any(l1)){
			H <- H[!l1, ]
			skip <- skip[!l1]
		}

	}
	colnames(H) <- c("lower end", "upper end", "integration step")
	rownames(H) <- paste("Subinterval", 1:nrow(H))
	if (length(skip) > 1L) names(skip) <- rownames(H)

	out <- list(H = H, skip = skip)	
	return(out)
}
