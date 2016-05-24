
.oneSample <- function(index, subsample, seed, Sx, Mx, Ax, Bx, cx, pattern, Ineq) {
  set.seed(seed[index])
  .determinant(Sx, Mx, Ax, Bx, cx, pattern, Ineq, subsample[index])
}

.oneSample_c <- function(index, subsample, seed, Sx, Mx, Ax, Bx, cx, pattern, Ineq, mConst) {
  set.seed(seed[index])
  .determinant_c(Sx, Mx, Ax, Bx, cx, pattern, Ineq, subsample[index], mConst)
}

bmpt.fia <- function(s, parameters, category, N, ineq0 = NULL, Sample = 2e+05, multicore = FALSE, split = NULL, mConst = NULL) {

	t0 <- Sys.time()
	print(paste("Computing FIA: Iteration begins at ", t0, sep = ""))
	flush.console()

	s <- strsplit(s, "")

	s <- s[[1]]

	if (any(!(s %in% c("p", "C")))) stop("Characters in the string must be either p or C.")

	type <- ifelse(s == "p", 1, 0)

	if (sum(type) != length(parameters)) stop("number of parameters not consistent in input arguments.")

	fixid <- (parameters <= 0)
	if (any(fixid)) parameter0 <- -parameters[fixid]
	else parameter0 <- 0
	parameters <- parameters[!fixid]

	if(any(floor(parameters) != parameters)) stop("use positive integers for free parameter assignment.")
	if(any(parameter0 < 0 | parameter0 > 1)) stop ("use negative values between 0 and 1 to fix a parameter.")
	if(sum(!type)!=length(category)) stop("number of categories not consistent in input arguments.")
	if (!is.null(ineq0)) if (dim(ineq0)[2] != 2) stop("ineq0 should have two columns.")

	######################

	L <- length (s)
	code <- matrix(0, L, L)

	if (type[1] == 0 & L == 1) return()
	if (type[1] == 0 & L != 1) stop("Not a BMPT model.")

	p <- 1
	u <- 1
	for (i in 2:L) {
		code[i,] <- code[p,]
		code[i,p] <- u
		if (type[i] == 1) {
			u <- 1
			p <- i
		} else {
			u <- -1
			ind <- i-1
			while (ind > 0) {
				if (ind <= 0 & i < L) stop( "Not a BMPT model.")
				
				if (type[ind] == 1) {
					p <- ind
					break
				} else {
					if (type[ind] == 0) {
						if (type[ind-1] !=1) stop("Not a BMPT model.")
						type[c(ind-1,ind, ind+1)] <- -1
						ind <- ind-2
					} else {
						if (type[ind] == -1) {
							type[ind+1] <- -1
							while (type[ind] == -1) ind <- ind-1
							if (type[ind] != 1) stop ("Not a BMPT model.")
							type[ind] <- -1
							ind <- ind-1
						}
					}
				}
			}
		}
	}

	if (ind > 0) stop ("Not a BMPT model.")

	code <- code[s == "C", s == "p"]

	##################

	code1 <- code[,!fixid, drop = FALSE]
	code0 <- code[,fixid, drop = FALSE]

	###################

	P1 <- length(parameters)
	assignment <- matrix(0, P1, P1)
	id = 1:P1
	pos <- matrix(0, 1, P1)

	i <- 1
	while (length(parameters) > 0) {
		a <- min(parameters)
		ind <- (parameters==a)
		assignment[id[ind],i] <- 1
		parameters <- parameters[-which(ind==TRUE)]
		id <- id[-which(ind==TRUE)]
		pos[i] <- a
		i <- i+1
	}
	assignment <- assignment[,1:(i-1)]
	pos <- pos[1:(i-1)]

	A <- (code1 == 1) %*% assignment 
	B <- (code1 == -1) %*% assignment 

	##################
	# fixed parameters
	# (this part was slightly changed for dealing with cases without fixed parameters, v 0.6.3)
	if (any(fixid)) {
		As <- matrix(1, dim(code0)[1], 1) %*% parameter0
		c <- apply((As^(code0==1)),1,prod)*apply((1-As)^(code0==-1),1,prod)
	} else {
		c <- rep(1, dim(code)[1])
	}

	#################
	if (!is.null(ineq0)) {
		nineq <- dim(ineq0)[1]
		Ineq <- matrix(0, nineq, dim(assignment)[2])
		for (i in 1:nineq) {
			Ineq[i, pos==ineq0[i,1]] <- -1  ####### erased transpose operator from originial code
			Ineq[i, pos==ineq0[i,2]] <- 1  ####### erased transpose operator from original code
		}
	} else Ineq <- matrix(0, 1, 1)

	###

	t1 <- dim(A)
	M <- t1[1]
	S <- t1[2]
	C <- max(category)
	pattern<-matrix(0, C, M)
	for (i in 1:C) {
		if (!any(category==i)) stop("argument category should involve consecutive numbers from 1 to the number of categories")
		pattern[i, category==i] <- 1
	}

  storage.mode(A) <- "integer"
  storage.mode(B) <- "integer"
	storage.mode(Sample) <- "integer"
  if (!is.null(split)) {
    subsamples <- rep(ceiling(Sample/split), split)
    seeds <- runif(split) * 10e7
  } else {
    if (multicore && eval(call("require", package = "snowfall", character.only = TRUE))) {
      subsamples <- vapply(snowfall::sfClusterSplit(1:Sample), length, 0)
      seeds <- runif(length(subsamples)) * 10e7
    } else {
      subsamples <- Sample
      seeds <- runif(1) * 10e7
    }
  }
  #browser()
  if (multicore) {
    snowfall::sfLibrary("MPTinR", character.only = TRUE)
    if (is.null(mConst)) tmp.dat <- snowfall::sfClusterApplyLB(1:length(seeds),.oneSample, Sx = S, Mx = M, Ax = A, Bx = B, cx = c, pattern = pattern, Ineq = Ineq, seed = seeds, subsample = subsamples)  
    else {
      storage.mode(mConst) <- "integer"
      tmp.dat <- snowfall::sfClusterApplyLB(1:length(seeds),.oneSample_c, Sx = S, Mx = M, Ax = A, Bx = B, cx = c, pattern = pattern, Ineq = Ineq, seed = seeds, subsample = subsamples, mConst = mConst)
    }
  } else 
    if (is.null(mConst)) tmp.dat <- lapply(1:length(seeds), .oneSample, Sx = S, Mx = M, Ax = A, Bx = B, cx = c, pattern = pattern, Ineq = Ineq, seed = seeds, subsample = subsamples) 
  else {
    storage.mode(mConst) <- "integer"
    tmp.dat <- lapply(1:length(seeds), .oneSample_c, Sx = S, Mx = M, Ax = A, Bx = B, cx = c, pattern = pattern, Ineq = Ineq, seed = seeds, subsample = subsamples, mConst = mConst) 
  }
  #browser()
  #str(tmp.dat)
	count <- sum(vapply(tmp.dat, "[", 0, i = 1))
	if (is.null(mConst)) {
    integral <- sum(vapply(tmp.dat, "[", 0, i = 3))
    vr <- sum(vapply(tmp.dat, "[", 0, i = 2))
	} else {
	  integral <- as.brob(sum(vapply(tmp.dat, "[", 0, i = 3))) / sqrt(as.brob(mConst) ^ dim(A)[2])
	  vr <- as.brob(sum(vapply(tmp.dat, "[", 0, i = 2))) / (as.brob(mConst) ^ dim(A)[2])
	}
	#sum(sapply(tmp.dat[2,], function(x) isTRUE(all.equal(x, 0))))
	# if (detI != 0) sample <- sample + 1	
  integral <- integral/Sample
	vr <- vr/Sample
	vr1 <- abs(vr-integral^2)/Sample
	lnInt <- log(integral)
	d1 <- sqrt(vr1)/integral*qnorm(.975)
	CI1 <- c(as.numeric(lnInt-d1), as.numeric(lnInt+d1))
	const <- Sample/count
	lnconst <- log(const)
	d2 <- sqrt((1-const)/Sample)* qnorm(.975)
	CI2 <- c(lnconst-d2, lnconst+d2)
	CFIA <- lnInt + lnconst+S/2*log(N/2/pi)
	d <- as.numeric(sqrt(d2^2 + d1^2))
	CI = c(CFIA-d, CFIA+d)
	out <- c(CFIA, CI, lnInt, CI1, lnconst, CI2)
	names(out) <- c("CFIA", "CI.l", "CI.u", "lnInt", "CI.lnint.l", "CI.lnint.u", "lnconst", "CI.lnconst.l", "CI.lnconst.u")
	t1 <- Sys.time()
	print(paste("Computing FIA: Iteration stopped at ", t1, sep = ""))
	print(t1-t0)
	
	out

}

