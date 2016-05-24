"kappa2" <-
function(ratings, weight = c("unweighted", "equal", "squared"), sort.levels = FALSE) {
	ratings <- as.matrix(na.omit(ratings))
	if (is.character(weight))
        weight = match.arg(weight)

	ns <- nrow(ratings)
	nr <- ncol(ratings)

  if (nr>2) {
    stop("Number of raters exeeds 2. Try kappam.fleiss or kappam.light.")
  }

	r1 <- ratings[,1]; r2 <- ratings[,2]

	if ((is.numeric(r1)) | (is.numeric(r2))) sort.levels <- TRUE
	
	if (!is.factor(r1)) r1 <- factor(r1)
	if (!is.factor(r2)) r2 <- factor(r2)

	#Find factor levels
	if (length(levels(r1)) >= length(levels(r2))) {
    lev <- c(levels(r1), levels(r2))
  } else { 
    lev <- c(levels(r2), levels(r1))
  }

	if (sort.levels) lev <- sort(lev)
	lev <- lev[!duplicated(lev)]
	r1 <- factor(ratings[,1],levels=lev)
	r2 <- factor(ratings[,2],levels=lev)

  #Compute table
	ttab <- table(r1, r2)

	#Compute weights
	nc <- ncol(ttab)

	if (is.numeric(weight))
		w <- 1-(weight-min(weight))/(max(weight)-min(weight))
  else if (weight == "equal")
    w <- (nc-1):0/(nc-1)
  else if (weight == "squared")
  	w <- 1 - (0:(nc-1))^2/(nc - 1)^2
  else #unweighted
    w <- c(1, rep(0,nc-1))

	wvec <- c(sort(w, decreasing=FALSE), w[2:length(w)])
	nw <- length(w)
	weighttab <- matrix(0, nrow=nw, ncol=nw)
	for (i in 1:nw) {
		weighttab[i,] <- wvec[(nw-(i-1)):(2*nw-i)]
	}

	agreeP <- sum(ttab*weighttab)/ns

	tm1 <- apply(ttab, 1, sum)
	tm2 <- apply(ttab, 2, sum)

	eij <- outer(tm1, tm2)/ns
	chanceP <- sum(eij*weighttab)/ns

	#Kappa for 2 raters
	value <- (agreeP - chanceP)/(1 - chanceP)

	#Compute statistics
	w.i <- apply(rep(tm2/ns,nc)*weighttab,2,sum)
	w.j <- apply(rep(tm1/ns,each=nc)*weighttab,1,sum)

	var.matrix <- (eij/ns)*(weighttab-outer(w.i,w.j,'+'))^2

  varkappa <- (sum(var.matrix)-chanceP^2)/(ns*(1-chanceP)^2)
  
	SEkappa <- sqrt(varkappa)
	u <- value/SEkappa
	p.value <- 2 * (1 - pnorm(abs(u)))

  rval <- structure(list(method = paste("Cohen's Kappa for 2 Raters (Weights: ",paste(weight,collapse=","),")",sep=""),
                         subjects = ns, raters = nr,
                         irr.name = "Kappa", value = value,
                         stat.name = "z", statistic = u, p.value = p.value),
                    class="irrlist")
  return(rval)
}

