qdc <- function(dataset, nfactors, zsc, sed) {
  if (nfactors==1) {
    qdc.res <- "Warning: Only one factor selected. No distinguishing and consensus statements will be calculated."
    warning(qdc.res)
  } else {
    # Distinguishing and consensus statements
    # create data frame
    comparisons <- combn(nfactors, 2, simplify=F)
    comp <- vector()
    for (i in 1:length(comparisons)) {
      comp <- append(comp, paste("f", comparisons[[i]], collapse="_", sep=""), after = length(comp))
    }
    qdc1 <- data.frame(matrix(data=as.numeric(NA), ncol=length(comp), nrow=nrow(dataset), dimnames=list(row.names(dataset), comp)))
    # differences in zsc between factors
    for (n in 1:length(comp)) {
      first <-  names(zsc)[grep(paste0("f", comparisons[[n]][1]), 
                                names(zsc))]
      second <- names(zsc)[grep(paste0("f", comparisons[[n]][2]), 
                                names(zsc))]
      qdc1[n] <- zsc[first] - zsc[second]
    }
    qdc2 <- as.data.frame(qdc1)
    # significant differences
    for (n in 1:length(comp)) {
      # find the threshold for the pair of factors
      sed <- data.frame(sed)
      first <-  names(sed)[grep(paste0("f", comparisons[[n]][1]), 
                                names(sed))]
      second <- names(sed)[grep(paste0("f", comparisons[[n]][2]), 
                                names(sed))]
      sedth.01 <- sed[first, second]*2.58
      sedth.05 <- sed[first, second]*1.96 # differences are significant when > 2.58*SED for p < .01, or the same value rounded upwards (Brown, 1980, pp.245)
      qdc2[which(abs(qdc1[[n]]) <= sedth.05), n] <- ""
      qdc2[which(abs(qdc1[[n]]) >  sedth.05), n] <- "*"
      qdc2[which(abs(qdc1[[n]]) >  sedth.01), n] <- "**"
    }
    names(qdc2) <- paste0("sig_",names(qdc2))
    qdc2$dist.and.cons <- as.character(apply(qdc2, 1, function(x) sum(x!="")==0))
    qdc2[which(qdc2$dist.and.cons == T), "dist.and.cons"] <- "Consensus"
    if (nfactors == 2) {
      qdc2[which(qdc2$dist.and.cons != "Consensus"), "dist.and.cons"] <- "Distinguishing"
    }
    if (nfactors > 2) {
      qdc2[which(qdc2$dist.and.cons != "Consensus"), "dist.and.cons"] <- ""
      for (i in 1:nfactors) {
        varsin  <- names(qdc2)[grep(i, names(qdc2))]
        varsout <- names(qdc2)[-grep(i, names(qdc2))]
        varsout <- varsout[-which(varsout=="dist.and.cons")]
        for (s in 1:nrow(qdc2)) {
          if (sum(qdc2[s, varsin] != "") == length(varsin) & sum(qdc2[s, varsout] != "") == 0) qdc2[s, "dist.and.cons"] <- paste0("Distinguishes f",i, " only") else if (sum(qdc2[s, c(varsin, varsout)] != "") == length(qdc1)) qdc2[s, "dist.and.cons"] <- "Distinguishes all" else if (sum(qdc2[s, varsin] != "") == length(varsin) & sum(qdc2[s, varsout] != "") != 0 & sum(qdc2[s, c(varsin, varsout)] != "") != length(qdc1)) qdc2[s, "dist.and.cons"] <- paste0(qdc2[s, "dist.and.cons"], "Distinguishes f",i, " ", collapse="")
        }
        #The above loop assigns these values in the column dist.and.cons, according to the following rules:
        # -- "Distinguishes f* only" when the differences of f* with all other factors are significant, AND all other differences are not.
        # -- "Distinguishes all" when all differences are significant.
        # -- "Distinguishes f*" when the differences of f* and all other factors are significant, AND some (but not all) of the other differences are significant.
        # -- "" leaves empty those which do not fullfil any of the above conditions, i.e. are not consensus neither are clearly distinguishing any factor.
      }
    }
    qdc.res <- cbind(qdc1, qdc2)
    ord <- rep(1:length(qdc1), each=2)
    ord[which(1:(length(qdc1)*2) %% 2 == 0)] <- ord[which(1:(length(qdc1)*2) %% 2 == 0)] + length(qdc1)
    qdc.res <- qdc.res[c(length(qdc.res), ord)]
  }
  return(qdc.res)
}