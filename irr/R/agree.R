"agree" <-
function(ratings, tolerance=0) {
  ratings <- as.matrix(na.omit(ratings))
  ns <- nrow(ratings)
  nr <- ncol(ratings)

  if (is.numeric(ratings)) {
    rangetab <- apply(ratings,1,max)-apply(ratings,1,min)
    coeff <- 100*sum(rangetab<=tolerance)/ns
  }
  else {
    rangetab <- as.numeric(sapply(apply(ratings,1,table),length))
    coeff <- 100*(sum(rangetab==1)/ns)
    tolerance <- 0
  }
	
  rval <- structure(list(method = paste("Percentage agreement (Tolerance=",tolerance,")",sep=""),
                         subjects = ns, raters = nr,
                         irr.name = "%-agree", value = coeff),
                    class="irrlist")
  return(rval)
}

