"robinson" <-
function(ratings) {
	ratings <- as.matrix(na.omit(ratings))
	ns <- nrow(ratings)
	nr <- ncol(ratings)

	SStotal <- var(as.numeric(ratings))*(ns*nr-1)
	SSb <- var(apply(ratings,1,mean))*nr*(ns-1)
	SSw <- var(apply(ratings,2,mean))*ns*(nr-1)
	SSr <- SStotal-SSb-SSw

	coeff  <- SSb/(SSb+SSr)

  rval <- structure(list(method = "Robinson's A",
                         subjects = ns, raters = nr,
                         irr.name = "A", value = coeff),
                    class="irrlist")
  return(rval)
}

