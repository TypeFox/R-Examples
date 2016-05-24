"kendall" <-
function(ratings, correct=FALSE) {
	ratings <- as.matrix(na.omit(ratings))

	ns <- nrow(ratings)
	nr <- ncol(ratings)

	#Without correction for ties
	if (!correct) {
		#Test for ties
		TIES = FALSE
		testties <- apply(ratings, 2, unique)
		if (!is.matrix(testties)) TIES=TRUE
		else { if (length(testties) < length(ratings)) TIES=TRUE }

		ratings.rank <- apply(ratings,2,rank)

		coeff.name <- "W"
		coeff <- (12*var(apply(ratings.rank,1,sum))*(ns-1))/(nr^2*(ns^3-ns))
	}
	else { #With correction for ties
		ratings <- as.matrix(na.omit(ratings))

		ns <- nrow(ratings)
		nr <- ncol(ratings)

		ratings.rank <- apply(ratings,2,rank)

		Tj <- 0
		for (i in 1:nr) {
			rater <- table(ratings.rank[,i])
			ties  <- rater[rater>1]
			l 	  <- as.numeric(ties)
			Tj	  <- Tj + sum(l^3-l)
		}

		coeff.name <- "Wt"
		coeff <- (12*var(apply(ratings.rank,1,sum))*(ns-1))/(nr^2*(ns^3-ns)-nr*Tj)
	}

	#test statistics
	Xvalue  <- nr*(ns-1)*coeff
	df1     <- ns-1
	p.value <- pchisq(Xvalue, df1, lower.tail = FALSE)

  rval <- list(method = paste("Kendall's coefficient of concordance",coeff.name),
               subjects = ns, raters = nr,
               irr.name = coeff.name, value = coeff,
               stat.name = paste("Chisq(",df1,")",sep=""), statistic = Xvalue, p.value = p.value)
  if (!correct && TIES) rval <- c(rval, error="Coefficient may be incorrect due to ties")
 	class(rval) <- "irrlist"
  return(rval)
}

