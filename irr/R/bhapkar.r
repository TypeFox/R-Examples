"bhapkar" <-
function(ratings) {
	ratings <- as.matrix(na.omit(ratings))

	ns <- nrow(ratings)
	nr <- ncol(ratings)

  if (nr>2) {
    stop("Number of raters exeeds 2. Coefficient cannot be calculated!")
  }

	r1 <- ratings[,1]; r2 <- ratings[,2]

	if (!is.factor(r1)) r1 <- factor(r1)
	if (!is.factor(r2)) r2 <- factor(r2)

	#Find factor levels
	if (length(levels(r1)) >= length(levels(r2))) {
    lev <- c(levels(r1), levels(r2))
  } else { 
    lev <- c(levels(r2), levels(r1))
  }

	lev <- lev[!duplicated(lev)]
	r1 <- factor(ratings[,1],levels=lev)
	r2 <- factor(ratings[,2],levels=lev)

  #Compute table
	ttab <- table(r1, r2)
  
  # get the marginals
  rowsums<-apply(ttab,1,sum)[1:(nrow(ttab)-1)]
  colsums<-apply(ttab,2,sum)[1:(ncol(ttab)-1)]
  # compute d matrix
  dmat <- matrix(rowsums-colsums,nrow=length(rowsums),ncol=length(rowsums),byrow=TRUE)
  # setup delta matrix
  delta <- matrix(0,nrow=nrow(ttab)-1,ncol=ncol(ttab)-1)
  diag(delta) <- rowsums+colsums
  # Dump last category from smx table
  smx <- ttab[1:(nrow(ttab)-1),1:(ncol(ttab)-1)]
  
  # Compute w matrix
  w  <- delta-smx-t(smx)-(dmat*t(dmat))/ns
  w1 <- solve(w)
  # Compute Chisq-value
  chimat <- dmat*t(dmat)*w1

	#test statistics
	Xvalue  <- sum(chimat)
	df1     <- nrow(ttab)-1
	p.value <- pchisq(Xvalue, df1, lower.tail = FALSE)

  rval <- list(method = paste("Bhapkar marginal homogeneity"),
               subjects = ns, raters = nr,
               irr.name = "Chisq", value = Xvalue,
               stat.name = paste("Chisq(",df1,")",sep=""), statistic = Xvalue, p.value = p.value)
 	class(rval) <- "irrlist"
  return(rval)
}

