"maxwell" <-
function(ratings) {
  ratings <- as.matrix(na.omit(ratings))
  ns <- nrow(ratings)
  nr <- ncol(ratings)

  if (nr>2) {
    stop("Number of raters exeeds 2")
  }
  
	r1 <- ratings[,1]; r2 <- ratings[,2]

	if (!is.factor(r1)) r1 <- factor(r1)
	if (!is.factor(r2)) r2 <- factor(r2)

	#Find factor levels
	if (length(levels(r1)) >= length(levels(r2))) lev <- c(levels(r1), levels(r2))
	else lev <- c(levels(r2), levels(r1))

	lev <- lev[!duplicated(lev)]
  if (length(lev)>2) {
    stop("Ratings are not binary")
  }
  
	r1 <- factor(ratings[,1],levels=lev)
	r2 <- factor(ratings[,2],levels=lev)

  # Compute table
	ttab <- table(r1, r2)

  # Compute coefficient
  coeff <- 2*sum(diag(ttab))/ns-1

  rval <- structure(list(method = "Maxwell's RE",
                         subjects = ns, raters = nr,
                         irr.name = "RE", value = coeff),
                    class="irrlist")
  return(rval)
}

