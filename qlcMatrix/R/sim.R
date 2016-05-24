# ==================================================
# some shortcuts for computing similarities directly
# ==================================================

# similarities between nominal attributes, i.e nominal variables
# this code could use some clean-up and harmonization :-)

sim.att <- function(D, method = "chuprov", sparse = TRUE) {
	
	X <- splitTable(D)
	
	# Chuprov's T, almost the same as Cramér's V, but easier to implement
	if (!is.na(pmatch(method,"chuprov"))) {

		r <- assocCol(X$OV, X$AV, method = res, sparse = sparse)
		if (!sparse) {
			r@x[is.na(r@x)] <- 0 # residuals can be NA when E==zero
		}

		X2 <- (X$AV*1) %*% r^2 %*% t(X$AV*1)
		N <- crossprod(tcrossprod(X$OV*1,X$AV*1))
		D <- Diagonal( x = sqrt(rowSums(X$AV) - 1) )
		R <- D %*% (as(N,"nMatrix")*1) %*% D

		if (sparse) {
			X2 <- as( X2 , "symmetricMatrix" )
			R <- as( R , "symmetricMatrix" )		
		} else {
			X2 <- as(as( X2, "dgCMatrix"), "symmetricMatrix" )
			R <- as(as( R, "dgCMatrix"), "symmetricMatrix" )
		}
		
		result <- N # just to get the right sparsity structure
		result@x <- sqrt( X2@x/(N@x * R@x) )
	}

	# The following options are highly similar, using the same base functions
	get_wpmi_assoc <- function(X) {
		r <- assocCol(X$OV, X$AV, method = wpmi)
		g <- (X$AV*1) %*% r %*% t(X$AV*1)
		g <- as(as( g, "dgCMatrix"), "symmetricMatrix" )
		return(g)
	}
	get_N <- function(X,g) {
		N <- crossprod(tcrossprod(X$OV*1,X$AV*1))
		if ( length(g@x) != length(N@x) ) {
			N <- N * (as(g,"nMatrix")*1)
		}
		return(N)
	}

	# G-test from Sokal and Rohlf (1981), also known as 'Dunning's G'
	# related to Mutual Information by a factor N
	if (!is.na(pmatch(method,"g-test"))) {
		g <- get_wpmi_assoc(X)		
		g@x <- 2*g@x
		result <- g
	}
	
	# Mutual Information
	if (!is.na(pmatch(method,"mutual information"))) {
		g <- get_wpmi_assoc(X)
		N <- get_N(X,g)
		g@x <- g@x/N@x
		result <- g
	}
	
	# Variation of Information = Mutual information metric
	if (!is.na(pmatch(method,"variation of information"))) {
		g <- get_wpmi_assoc(X)
		N <- get_N(X,g)
		g@x <- g@x/N@x

		O <- crossprod(X$OV*1)
		H1 <- (X$AV*1) %*% O %*% t(X$AV*1)
		H1 <- as(H1, "symmetricMatrix")
		H1@x <- H1@x * log(N@x) / N@x
		
		O@x <- O@x * log(O@x)
		H2 <- (X$AV*1) %*% O %*% t(X$AV*1)
		H2 <- as(H2, "symmetricMatrix")
		H2@x <- H2@x / N@x
		
		H <- g # just to get the right sparsity structure
		H@x <- (H1@x - H2@x - g@x)
		result <- H
	}
	
	rownames(result) <- X$attributes
	return(result)

}

# similarities between observations from nominal data
# this is a very simple wrapper around cosRow and assocRow

sim.obs <- function(D, method = "hamming", sparse = TRUE) {
	
	X <- splitTable(D)
	
	# Relative Hamming similarity (Goebl's "Relativer Identitätswert"), i.e. the number of similarities divided by the number of comparisons made
	if (!is.na(pmatch(method,"hamming"))) {
		result <- cosRow(t(X$OV), t(X$AV), norm = norm1)

	# weighted similarity very similar to Goebl's "Gewichteter Identitätswert". Note that his definition is slightly different, but that one is tricky to replicate
	} else	if (!is.na(pmatch(method,"weighted"))) {
		result <- cosRow(t(X$OV), t(X$AV), norm = norm2, weight = isqrt)

	#assoc methods
	} else {
		result <- assocRow(t(X$OV), t(X$AV), method = method)
	}

	rownames(result) <- X$observations
	return(result)
	
}

# similarity for words in parallel text. If weight is specified, method is ignored: cosSparse is used with norm2 and specified weight
# best uses rowMax/colMax, which is not very quick

sim.words <- function(text1, text2 = NULL, method = res, weight = NULL, 
						lowercase = TRUE, best = FALSE, tol = 0) {

	if (is.null(text2)) {
		T1 <- splitText( text1, simplify = TRUE, lowercase = lowercase )
		# compute co-occurrence statistics
		if (!is.null(weight)) {
			R <- cosSparse( t(T1), weight =  weight )
		} else {
			R <- assocSparse( t(T1), method =  method )
		}
	} else {
		globalID <- union(names(text1), names(text2))

		T1 <- splitText( text1, globalID, simplify = TRUE, lowercase = lowercase )
		T2 <- splitText( text2, globalID, simplify = TRUE, lowercase = lowercase )

		# collapse verses in which one of the translation is empty (i.e. combined translation of multiple verses into one verse)
		m1 <- which(text1 == "")
		m2 <- which(text2 == "")
		m <- union(m1, m2)
		
		M <- Diagonal(n = length(globalID))
		while ( sum(M[,m]) > 0 ) {
			tmp <- M[,m]
			M[,m] <- 0
			M[,(m-1)] <- M[,(m-1)] + tmp	
		}
		M <- M[,colSums(M)>0]
		
		# remap to collapse verses
		T1 <- T1 %*% M
		T2 <- T2 %*% M
	
		# compute co-occurrence statistics
		if (!is.null(weight)) {
			R <- cosSparse( t(T1), t(T2), weight =  weight )
		} else {
			R <- assocSparse( t(T1), t(T2), method =  method )
		}
	}
	
	if (tol > 0) {
		R <- drop0(R, tol = tol)
	}
	
	if (best) {
		choice <-	colMax(R, which = TRUE, ignore.zero = FALSE)$which +
					rowMax(R, which = TRUE, ignore.zero = FALSE)$which
		choice <- as(choice, "nMatrix")
		return(list(sim = R, best = choice))
	} else {
		return(R)
	}	
}

# quick string comparison based on cosine similarity between bigrams

sim.strings <- function(strings1, strings2 = NULL, sep = "", boundary = TRUE) {

	S1 <- splitStrings(strings1, sep = sep, boundary = boundary, simplify = TRUE)

	if (is.null(strings2)) {
		sim <- cosSparse(S1)
	} else {
		S2 <- splitStrings(strings2, sep = sep, boundary = boundary, simplify = TRUE)
		M <- jMatrix( rownames(S1), rownames(S2) )
		sim <- cosSparse( (M$M1*1) %*% S1, (M$M2*1) %*% S2 )
	}
	return(drop(sim))
}

# various similarities for wordlists
# sim.graph: similarity between graphemes, based on cooccurrences in context

sim.graph <- function(
		wordlist,
		doculects = "DOCULECT", concepts = "CONCEPT", counterparts = "TOKENS",
		method = "cooccurrence", assoc.method = poi, weight = NULL, sep = " "
		) {	
			
	W <- splitWordlist(
		wordlist, doculects =  doculects, concepts = concepts, counterparts = counterparts, sep = sep
		)
	CG <- (W$CW*1) %*% t(W$SW*1) %*% t(W$GS*1)		
	if (!is.null(weight)) {
			sim <- cosSparse( CG, weight =  weight )
		} else {
			sim <- assocSparse( CG, method = assoc.method )
		}
	rownames(sim) <- W$graphemes

	# additional matrix to identify the graphemes per language
	# without needing to parse the rownames...
	GD <- W$GS %*% W$SW %*% t(W$DW)
	colnames(GD) <- W$doculects
	rownames(GD) <- W$graphemes

	return(list(GG = sim, GD = GD))
}

# sim.lang: similarity between languages

sim.lang <- function(
		wordlist, 
		doculects = "DOCULECT", concepts = "CONCEPT", counterparts = "COUNTERPART",
		method = "parallel", assoc.method =  res, weight = NULL, sep = ""
		) {	
			
	W <- splitWordlist(
		wordlist, doculects = doculects, concepts = concepts, counterparts = counterparts, sep = sep
		)	

	if (!is.na(pmatch(method,"global"))) {
		BD <- (W$BS*1) %*% (W$SW*1) %*% t(W$DW*1)
		if (!is.null(weight)) {
			sim <- cosSparse( BD, weight =  weight )
		} else {
			sim <- assocSparse( BD, method = assoc.method )
		}
	}	
	if (!is.na(pmatch(method,"parallel"))) {
		BW <- (W$BS*1) %*% (W$SW*1)
		CBxW <- KhatriRao(BW, (W$CW*1))
		CBxD <- CBxW %*% t(W$DW*1)
		if (!is.null(weight)) {
			sim <- cosSparse( CBxD, weight =  weight )
		} else {
			sim <- assocSparse( CBxD, method = assoc.method )
		}
	} 	
	colnames(sim) <- rownames(sim) <- W$doculects
	return(sim)
}

# sim.con: Similarity between concepts

sim.con <- function(
		wordlist,
		doculects = "DOCULECT", concepts = "CONCEPT", counterparts = "COUNTERPART",
		method = "bigrams", assoc.method = res, weight = NULL, sep = ""
		) {
	if (!is.na(pmatch(method,"colexification"))) {
		W <- splitWordlist(
			wordlist, doculects = doculects, concepts = concepts, counterparts = counterparts, 
			splitstrings = FALSE, simplify = FALSE
			)
		sim <- tcrossprod(W$CW*1)
	}
	if (!is.na(pmatch(method,"global"))) {
		W <- splitWordlist(
			wordlist, doculects = doculects, concepts = concepts, counterparts = counterparts, 
			sep = sep
			)
		BC <- (W$BS*1) %*% (W$SW*1) %*% t(W$CW*1)
		if (!is.null(weight)) {
			sim <- cosSparse( BC, weight =  weight )
		} else {
			sim <- assocSparse( BC, method = assoc.method )
		}		
	}	
	if (!is.na(pmatch(method,"bigrams"))) {
		W <- splitWordlist(
			wordlist, doculects = doculects, concepts = concepts, counterparts = counterparts, 
			sep = sep
			)
		TC <- (W$TS*1) %*% (W$SW*1) %*% t(W$CW*1)
		if (!is.null(weight)) {
			sim <- cosSparse( TC, weight =  weight )
		} else {
			sim <- assocSparse( TC, method = assoc.method )
		}		
	}		
	colnames(sim) <- rownames(sim) <- W$concepts
	return(sim)		
}