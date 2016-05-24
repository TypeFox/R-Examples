# 2.18.2009
# Compute trait metrics from a bcnt.otu file

trait.stat <- function(bcnt.otu, coefs) {

  # identify unique traits
  if (is.factor(coefs$TRAITVAL)) {
    temp <- levels(coefs$TRAITVAL)
  }
  else {
    if (is.character(coefs$TRAITVAL)) {
      temp <- unique(coefs$TRAITVAL)
    }
  }
  
  traitval.u <- character(0)
  repeat{
    w <- regexpr("[A-Za-z]+", temp)
    if (sum(w!=-1) == 0) break
    traitval.u <- c(traitval.u,substring(temp, w, w+attributes(w)$match.length-
                                          1))
    temp <- substring(temp, w+attributes(w)$match.length+1, nchar(temp))
  }
  traitval.u[traitval.u == ""] <- NA
  traitval.u <- sort(unique(traitval.u))

  if (!is.factor(bcnt.otu[,1])) {
    bcnt.otu[,1] <- factor(bcnt.otu[,1])
  }

  bcnt.t <- merge(bcnt.otu, coefs, by.x = "OTU", by.y="TAXON", all.x = T)

  matout <- matrix(NA, ncol = 3*length(traitval.u),
                   nrow = length(levels(bcnt.otu[,1])))

  if (is.factor(coefs$TRAITVAL)) {
    all.traits <- levels(bcnt.t$TRAITVAL)[bcnt.t$TRAITVAL]
  }
  else {
    if (is.character(coefs$TRAITVAL)) {
      all.traits <- bcnt.t$TRAITVAL
    }
  }

  fsite <- names(bcnt.otu)[1]
  fabund <- names(bcnt.otu)[3]
  ftaxon <- names(bcnt.otu)[2]

  for (i in 1:length(traitval.u)) {
    incvec <- regexpr(traitval.u[i], all.traits) != -1
    tr.rich <- tapply(bcnt.t$OTU[incvec], bcnt.t[incvec,fsite],
                      function(x) length(unique(x)))
    tr.rich[is.na(tr.rich)] <- 0
    matout[,(i-1)*3+1] <- tr.rich
    trich <- tapply(bcnt.otu[, ftaxon], bcnt.otu[,fsite],
                    function(x) length(unique(x)))
    matout[,(i-1)*3+2] <- tr.rich/trich
    totabund <- tapply(bcnt.otu[,fabund], bcnt.otu[,fsite], sum)
    tr.abund <- tapply(bcnt.t[incvec,fabund], bcnt.t[incvec,fsite], sum)
    tr.abund[is.na(tr.abund)] <- 0
    matout[, (i-1)*3+3] <- tr.abund/totabund
  }

  names0 <- paste(rep(traitval.u, times = rep(3, times = length(traitval.u))),
                  rep(c(".RICH", ".PTAX", ".PABN"), times = length(traitval.u)),
                  sep = "")
  dfout <- data.frame(id = levels(bcnt.otu[,1]), matout)
  names0 <- c(fsite, names0)
  names(dfout) <- names0
  return(dfout)
}

