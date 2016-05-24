unidTab <- function(lev, fac, algo = 'cd2') {
	##########################################
	# parameter initialize 
	# algo <- cd2 // sd2

	##############################################
	## centered L2-discrepancy 
	uniForm <- function(U, algo = 'cd2') {
		size <- dim(U)
		lev <- size[1]
		fac <- size[2]
		
		if (algo == 'cd2'){
			U <- (U - 1/2) / lev
			csum.one <- 0
			csum.two <- 0
			for (i in seq(lev)){
				cprod <- prod((2 + abs(U[i, ] - 1/2) - (U[i, ] - 1/2)^2))
				csum.one <- csum.one + cprod	 
			}
			
			for (i in seq(lev)){
				for(j in seq(lev)){
					cprod <- prod(1 + 1/2 * ((abs(U[i, ] - 1/2)) + abs(U[j, ] - 1/2) -
													abs(U[i, ] - U[j, ])))
					csum.two <- csum.two + cprod
				}
			}
			csum <- sqrt((13/12)^fac - 2^(1 - fac) / lev * csum.one + 1/lev^2 * csum.two)
			return(csum)
			
		}else if (algo == 'sd2') {
			U <- (U - 1/2) / lev	
			csum.one <- 0
			csum.two <- 0
			
			for (i in seq(lev)){
				cprod <- prod(1 + 2 * U[i, ] - 2 * U[i, ]^2)
				csum.one <- csum.one + cprod
			}
			
			for (i in seq(lev)){
				for (j in seq(lev)){
					cprod <- prod(1 - abs(U[i, ] - U[j, ]))
					csum.two <- csum.two + cprod
				}
			}
			csum <- sqrt((4/3)^fac - 2 / lev * csum.one + 2^fac / lev^2 * csum.two)
			return(csum)
		}
	}

	##############################################
	splitMod <- function(fac, MARK){
		# generate the split mod matrix
		# mod(x^n, z) <- mod(mod(x, z)^n, z)
		# mod(x^2n, z) <- mod(mod(x^n, z) * mod(x^n, z), z)
		sec <- round(seq(0, (fac - 1), length.out = MARK))
		pst <- as.integer(sec + 1)
		secLen <- length(sec)
		mat <- matrix(0, MARK, fac)
		temp <- rep(0, fac)
		temp.one <- temp
		temp.all <- temp
		
		for (i in seq(secLen - 1)){
		  temp.one[pst[i] : (pst[i + 1] - 1)] <- seq(sec[i], (sec[i + 1] - 1))
		  temp.one[pst[i + 1] : fac] <- sec[i + 1]
		  temp.two <- temp.one - temp
		  mat[i, ] <- temp.two
		  temp <- temp.one
		  temp.all <- temp.all + temp.two
		  
		  if (i == secLen - 1){
			mat[secLen, ] <- seq(0, (fac - 1)) - temp.all
			break
		  }
		}
		
		return(mat)
	}

	##############################################
	## mod
	gcd <- function(a, b) ifelse (b == 0, a, gcd(b, a %% b))

	##############################################
	makeTab <- function(lev, fac, adj = FALSE){
		# generate the uni-table
		
		levSeq <- seq(lev)
		# first row in uni-table
		rowInit <- which(sapply(levSeq, function(x) gcd(x, lev)) == 1)
		rowLen <- length(rowInit)
		tab <- crossprod(t(levSeq), t(rowInit)) %% lev
		
		if (adj == TRUE) {
			tab <- tab[-lev, ]
		}else {
			tab[lev, ] <- lev
		}
		
		return(tab)
	}

	##############################################
	getTab <- function(lev, fac, UPLIMIT, algo, adj = FALSE){
		# choose optimal uni-table
		# MARK <- 5
		levSeq <- seq(lev)
		count <- choose(lev, fac)
		
		if (count < UPLIMIT){
			tab <- makeTab(lev, fac, adj)
			rowLen <- ncol(tab)
			option <- as.matrix(combn(seq(rowLen), fac))
			# choose options start with 1 (first level)
			option <- as.matrix(option[ , (option[1, ] == 1)])
			optionNum <- dim(option)[2]
			dev <- rep(0, optionNum)

			for (i in seq(optionNum)){
				tab.x <- tab[, option[, i]]
				# two different metrics
				dev[i] <- uniForm(tab.x, algo)
			}
			
			idx <- which(dev == min(dev))
			minDev <- dev[idx]
			bestOption = as.matrix(option[, idx])
			idxNum <- length(idx)
			
			if (adj == TRUE) {
				uniTab <- array(0, dim = c((lev - 1), fac, idxNum))
			}else {
				uniTab <- array(0, dim = c(lev, fac, idxNum))
			}
			
			for (k in seq(idxNum)){
				uniTab[, , k] <- tab[, bestOption[, k]]
			}
			
		}else {
			# MARK: stepwise power, use factor as an indicator
			
			if (fac <= 30) {
				MARK <- 8
			}else if (fac <= 60) {
				MARK <- 15
			}else if (fac <= 100) {
				MARK <- 25
			} else {
				MARK <- 98
			}
			
			# store the deviance
			dev <- rep(0, lev)
			
			for (i in seq(lev)){
				# mod(x^n, z) <- mod(mod(x, z)^n, z)
				# mod(x^2n, z) <- mod(mod(x^n, z) * mod(x^n, z), z)
				## this is a unique point ****
				
				if (fac < 7) {
					tab <- crossprod(t(levSeq) %% lev, t(i^(seq(0, (fac - 1)))) %% lev) %% lev
				}else {
					# mod(x^2n, z) <- mod(mod(x^n, z) * mod(x^n, z), z)
					# split the mod of power into 4 parts
					mat <- splitMod(fac, MARK)
					powermod <- 1
					for (j in seq(MARK)){
						powermod <- powermod * (i^mat[j, ])
						#print(powermod)
						powermod <- powermod %% lev
						#print(powermod)
					}
					tab <- crossprod(t(levSeq) %% lev, t(powermod) %% lev) %% lev
				}
				
				if (adj == TRUE) {
					tab <- tab[-lev, ]
				}else {
					tab[lev, ] <- lev
				}
				dev[i] <- uniForm(tab, algo)
			}
			
			idx <- which(dev == min(dev))
			idxNum <- length(idx)
			minDev <- dev[idx]
			
			if (adj == TRUE) {
				uniTab <- array(0, dim = c((lev - 1), fac, idxNum))
			}else {
				uniTab <- array(0, dim = c(lev, fac, idxNum))
			}
			
			for (i in seq(idxNum)){
				# extreme case of level 100 and fac 100
				
				if (fac < 7) {
					tab <- crossprod(t(levSeq) %% lev, t(idx[i]^(seq(0, (fac - 1)))) %% lev) %% lev
				}else {
					# mod(x^2n, z) <- mod(mod(x^n, z) * mod(x^n, z), z)
					# split the mod of power into 4 parts
					mat <- splitMod(fac, MARK)
					powermod <- 1
					for (j in seq(MARK)){
						powermod <- powermod * idx[i]^mat[j, ]
						powermod <- powermod %% lev
					}
					tab <- crossprod(t(levSeq) %% lev, t(powermod) %% lev) %% lev
				}
				
				if (adj == TRUE) {
					tab <- tab[-lev, ]
				}else {
					tab[lev, ] <- lev
				}
				uniTab[, , i] <- tab
			}
		}
		
		list(T = uniTab, D = minDev)
	}

##########################################
	# output the standard uni-table
	levo <- lev + 1
	levSeq.orn <- seq(lev)
	levSeq.adj <- seq(levo)
	# first row in uni-table
	rowInit.orn <- which(sapply(levSeq.orn, function(x) gcd(x, lev)) == 1)
	rowInit.adj <- which(sapply(levSeq.adj, function(x) gcd(x, levo)) == 1)
	rowLen.orn <- length(rowInit.orn)
	rowLen.adj <- length(rowInit.adj)

	flag = ""
	UPLIMIT = 10000
	
	if (rowLen.orn == rowLen.adj){
		# cases like level 3 and 15
		
		if (rowLen.orn > fac){
			# check the factor # and actual level #
			flag <- 'b'
			tab.orn <- getTab(lev, fac, UPLIMIT, algo)
			tab.adj <- getTab(levo, fac, UPLIMIT, algo, adj = TRUE)
			
		}else if (rowLen.orn == fac){
			flag <- 'b'
			tab.orn.T <- makeTab(lev, fac)
			tab.orn.D <- uniForm(tab.orn.T, algo)
			tab.adj.T <- makeTab(levo, fac, adj = TRUE)
			#tab.adj.T <- tab.adj.T[-levo, ]
			tab.adj.D <- uniForm(tab.adj.T, algo)
			tab.orn <- list(T = tab.orn.T, D = tab.orn.D)
			tab.adj <- list(T = tab.adj.T, D = tab.adj.D)
			
		}else {
			flag <- 'n'
			cat('The factor: \"', fac, '\"','can not match level: \"', lev, '\"')
			stop("problem in the number of factor")
		}
		
	}else if (rowLen.orn > rowLen.adj){
		# cases like level 5 and level 7
		if (rowLen.adj >= fac){
			# check the factor # and actual level #
			flag <- 'b'
			tab.orn <- getTab(lev, fac, UPLIMIT, algo)
			
			if (rowLen.adj == fac) {
				tab.adj.T <- makeTab(levo, fac, adj = TRUE)
				# tab.adj.T <- tab.adj.T[-levo, ]
				tab.adj.D <- uniForm(tab.adj.T, algo)
				tab.adj <- list(T = tab.adj.T, D = tab.adj.D)
			}else {
				tab.adj <- getTab(levo, fac, UPLIMIT, algo, adj = TRUE)
			}
			
		}else if (rowLen.orn >= fac){
			flag <- 'o'	
			
			if (rowLen.orn == fac){
				tab.orn.T <- makeTab(lev, fac)
				tab.orn.D <- uniForm(tab.orn.T, algo)
				tab.orn <- list(T = tab.orn.T, D = tab.orn.D)		
			}else {
				tab.orn <- getTab(lev, fac, UPLIMIT, algo)
			}
			
		}else {
			flag <- 'n'
			cat('The factor: \"', fac, '\"','can not match level: \"', lev, '\"')
			stop("problem in the number of factor")	
		}
		
	}else if (rowLen.orn < rowLen.adj){
		# cases like level 4 and level 6
		if (rowLen.orn >= fac){
			# check the factor # and actual level #
			flag <- 'b'
			
			if (rowLen.orn == fac) {
				tab.orn.T <- makeTab(lev, fac)
				tab.orn.D <- uniForm(tab.orn.T, algo)
				tab.orn <- list(T = tab.orn.T, D = tab.orn.D)		
			}else {
				tab.orn <- getTab(lev, fac, UPLIMIT, algo)
			}
			
			tab.adj <- getTab(levo, fac, UPLIMIT, algo, adj = TRUE)
			
		}else if (rowLen.adj >= fac){
			flag <- 'a'
			
			if (rowLen.adj == fac) {
				tab.adj.T <- makeTab(levo, fac, adj = TRUE)
				tab.adj.D <- uniForm(tab.adj.T, algo)
				tab.adj <- list(T = tab.adj.T, D = tab.adj.D)
			}else {
				tab.adj <- getTab(levo, fac, UPLIMIT, algo, adj = TRUE)
			}
			
		}else {
			flag <- 'n'
			cat('The factor: \"', fac, '\"','can not match level: \"', lev, '\"')
			stop("problem in the number of factor")	
		}
	}
	
	if (exists("tab.orn") && length(dim(tab.orn$T)) == 3){
	
		if (dim(tab.orn$T)[3] == 1) {
			tab.orn$T <- tab.orn$T[, , 1]
		}else {
			tab.orn$D <- tab.orn$D[1]
		}
	}

	if (exists("tab.adj") && length(dim(tab.adj$T)) == 3){
	
		if (dim(tab.adj$T)[3] == 1) {
			tab.adj$T <- tab.adj$T[, , 1]
		}else {
			tab.adj$D <- tab.adj$D[1]
		}
	}
	
	# select the final standard uni-table with smaller deviance 
	if (flag == "b") {
	
		if (tab.orn$D <= tab.adj$D)
			return(tab.orn)
		else
			return(tab.adj)
			
	}else if (flag == "o"){
		return(tab.orn)
		
	}else if (flag == "a") {
		return(tab.adj)
	}
}
