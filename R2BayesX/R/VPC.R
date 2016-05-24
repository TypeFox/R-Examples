VPC <- function(x, term = NULL)
{
  stopifnot(inherits(x, "bayesx"))
  
  vpc <- NULL

  if(!is.null(sigma2 <- x$variance[, 1])) {
    j <- if(is.null(term)) 1 else {
      if(is.character(term)) pmatch(term, names(x$effects)) else term 
    }

    vpc <- if(!is.null(sigma2term <- attr(x$effects[[j]], "variance"))) {
      sigma2term <- sigma2term[, 1]
      sigma2 / (sigma2 + sigma2term)
    }
  }

  vpc
}

#data("ZambiaBnd")
#data("ZambiaNutrition")

#f <- stunting ~ memployment + urban + gender + meducation + sx(mbmi) +
#  sx(agechild) + sx(district, bs = "mrf", map = ZambiaBnd) +
#  sx(district, bs = "re")

#m <- bayesx(f, family = "gaussian", method = "MCMC", iterations = 12000,
#  burnin = 2000, step = 10, seed = 123, data = ZambiaNutrition)

