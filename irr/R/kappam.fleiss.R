"kappam.fleiss" <-
function(ratings, exact = FALSE, detail = FALSE) {
	ratings <- as.matrix(na.omit(ratings))

	ns <- nrow(ratings)
	nr <- ncol(ratings)

	#Build table
	lev <- levels(as.factor(ratings))

	for (i in 1:ns) {
		frow <- factor(ratings[i,],levels=lev)

		if (i==1)
			ttab <- as.numeric(table(frow))
		else
			ttab <- rbind(ttab, as.numeric(table(frow)))
	}

	ttab <- matrix(ttab, nrow=ns)

	agreeP <- sum((apply(ttab^2,1,sum)-nr)/(nr*(nr-1))/ns)

  if (!exact) {
    method  <- "Fleiss' Kappa for m Raters"
  	chanceP <- sum(apply(ttab,2,sum)^2)/(ns*nr)^2
	} else {
    method  <- "Fleiss' Kappa for m Raters (exact value)"
  	for (i in 1:nr) {
  		rcol <- factor(ratings[,i],levels=lev)

  		if (i==1)
  			rtab <- as.numeric(table(rcol))
  		else
  			rtab <- rbind(rtab, as.numeric(table(rcol)))
  	}

  	rtab <- rtab/ns

  	chanceP <- sum(apply(ttab,2,sum)^2)/(ns*nr)^2 - sum(apply(rtab,2,var)*(nr-1)/nr)/(nr-1)
  }

 	#Kappa for m raters
 	value <- (agreeP - chanceP)/(1 - chanceP)

  if (!exact) {
    pj <- apply(ttab,2,sum)/(ns*nr)
    qj <- 1-pj
    
    varkappa <- (2/(sum(pj*qj)^2*(ns*nr*(nr-1))))*(sum(pj*qj)^2-sum(pj*qj*(qj-pj)))
  	SEkappa <- sqrt(varkappa)

  	u <- value/SEkappa
  	p.value <- 2 * (1 - pnorm(abs(u)))

    if (detail) {
    	pj  <- apply(ttab,2,sum)/(ns*nr)
    	pjk <- (apply(ttab^2,2,sum)-ns*nr*pj)/(ns*nr*(nr-1)*pj)

    	kappaK    <- (pjk-pj)/(1-pj)

    	varkappaK <- 2/(ns*nr*(nr-1))
    	SEkappaK  <- sqrt(varkappaK)

    	uK <- kappaK/SEkappaK
    	p.valueK <- 2 * (1 - pnorm(abs(uK)))

    	tableK <- as.table(round(cbind(kappaK, uK, p.valueK), digits=3))

    	rownames(tableK) <- lev
    	colnames(tableK) <- c("Kappa", "z", "p.value")
    }
  }

  if (!exact) {
    if (!detail) {
      rval <- list(method = method,
                   subjects = ns, raters = nr,
                   irr.name = "Kappa", value = value)
    }
    else {
      rval <- list(method = method,
                   subjects = ns, raters = nr,
                   irr.name = "Kappa", value = value,
                   detail = tableK)
    }
    rval <- c(rval, stat.name = "z", statistic = u, p.value = p.value)
  }
  else {
    rval <- list(method = method,
                 subjects = ns, raters = nr,
                 irr.name = "Kappa", value = value)
  }
  class(rval) <- "irrlist"

  return(rval)
}

