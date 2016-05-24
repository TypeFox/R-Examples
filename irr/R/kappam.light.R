"kappam.light" <-
function(ratings) {
	ratings <- as.matrix(na.omit(ratings))

	ns <- nrow(ratings)
	nr <- ncol(ratings)

  for (i in 1:(nr-1))
    for (j in (i+1):nr) {
      if ((i==1) & (j==(i+1))) kappas <- kappa2(ratings[,c(i,j)], weight="u")$value
      else kappas <- c(kappas, kappa2(ratings[,c(i,j)], weight="u")$value)
    }

  value <- mean(kappas)

  #Variance & Computation of p-value
  lev    <- levels(as.factor(ratings))
  levlen <- length(levels(as.factor(ratings)))

  for (nri in 1:(nr-1))
    for (nrj in (nri+1):nr) {
      for (i in 1:levlen)
        for (j in 1:levlen) {
          if (i!=j) {
            r1i <- sum(ratings[,nri]==lev[i])
            r2j <- sum(ratings[,nrj]==lev[j])
            if (!exists("dis")) dis <- r1i*r2j
            else dis <- c(dis,r1i*r2j)
          }
        }
        if (!exists("disrater")) disrater <- sum(dis)
        else disrater <- c(disrater,sum(dis))
        rm(dis)
      }

  B <- length(disrater) * prod(disrater)

  chanceP  <- 1-B/ns^(choose(nr,2)*2)
  varkappa <- chanceP/(ns*(1-chanceP))

	SEkappa <- sqrt(varkappa)
	u <- value/SEkappa
	p.value <- 2 * (1 - pnorm(abs(u)))

  rval <- structure(list(method = "Light's Kappa for m Raters",
                         subjects = ns, raters = nr,
                         irr.name = "Kappa", value = value,
                         stat.name = "z", statistic = u, p.value = p.value),
                    class="irrlist")
  return(rval)
}

