strans <- function(M){
  # Check stochastic transitivities in a paired-comparison matrix (abs. freq.)
  # last mod: 08/Oct/2003 works now for unbalanced design
  #           21/Aug/2007 log the triples that violate transitivity
  #           14/Feb/2009 re-work strans: chkdf, viol.tab, re-ordered pcm
  #           18/Feb/2009 approx. LRT of WST
  #           30/Sep/2009 replace viol.tab by violdf
  #                       obj.names for I > 26
  #           24/May/2013 pre-allocate vectors in chkdf
  # author: Florian Wickelmaier (wickelmaier@web.de)

  I <- sqrt(length(as.matrix(M)))    # number of stimuli
  R <- as.matrix( M / (M + t(M)) )   # pcm rel. freq.
  n <- (M + t(M))[lower.tri(R)]      # number of obs per pair
  R[which(is.na(R), arr.ind=TRUE)] <- 0  # replace diag (and missing cells)
  obj.names <- if(!is.null(colnames(M))) colnames(M)
               else make.unique(rep(letters, length.out=I), sep="")

  triple <- rep(seq_len(choose(I, 3)), each=6)
  perm <- character(length(triple))
  p12 <- p23 <- p13 <- numeric(length(triple))

  ## For each triple, go through the 6 permutations
  m <- 1
  for(ii in 1:(I-2)){ for(jj in (ii+1):(I-1)){ for(kk in (jj+1):I){ # tpl loop
    for(i in c(ii, jj, kk)){  # permutation loops
      for(j in c(ii, jj, kk)){
        for(k in c(ii, jj, kk)){
          if(i != j && j != k && i != k){
            perm[m] <- paste(obj.names[i], obj.names[j], obj.names[k],
                             sep=".")
            p12[m] <- R[i, j]
            p23[m] <- R[j, k]
            p13[m] <- R[i, k]
            m <- m + 1
          }
        }
      }
    }
  } } }
  chkdf <- data.frame(triple, perm, p12, p23, p13)

  ## Data frame of permutations that satisfy the premise
  predf <- chkdf[chkdf$p12 >=.5 & chkdf$p23 >= .5,]
  predf$wst <- predf$p13 >= .5
  predf$mst <- predf$p13 >= pmin(predf$p12, predf$p23)
  predf$sst <- predf$p13 >= pmax(predf$p12, predf$p23)
  
  ## Order to minimize the magnitude of violations
  predf <- predf[order(predf$p13, decreasing=TRUE),]
  
  ## Reduce predf to a single entry per triple
  violdf <- predf[gsub("^.+\\.(.+)$", "\\1",
    names(sapply(   # get the rownames of chkdf
      split(rowSums(predf[,c("wst","mst","sst")]), predf$triple),  # split
    which.max))),]  # which permutation satisfies the most transitivities?

  wv <- .5 - violdf[!violdf$wst, "p13"]                   # deviation from
  mv <- with(violdf[!violdf$mst,], pmin(p12, p23) - p13)  #   minimum
  sv <- with(violdf[!violdf$sst,], pmax(p12, p23) - p13)  #   permis prob

  # Re-order R such that (as much as possible given violations) entries
  # increase from left to right and from bottom to top
  rownames(R) <- colnames(R) <- obj.names
  # idx <- names(sort(rowSums(R),       decreasing=TRUE))
    idx <- names(sort(rowSums(R >= .5), decreasing=TRUE))  # binary matrix
  R <- R[idx, idx]

  ## Likelihood ratio test of WST
  p.obs <- t(R)[lower.tri(R)]  # observed proportions
  p.wst <- p.obs               # wst corrected proportions
  p.wst[p.obs < .5] <- .5

  wst.mat <- matrix(0, I, I)
  wst.mat[lower.tri(wst.mat)] <- round(p.wst*n)
  wst.mat <- t(wst.mat)
  wst.mat[lower.tri(wst.mat)] <- n - round(p.wst*n)
  dimnames(wst.mat) <- dimnames(R)

  wst.fit <- c(
    G2   <- 2*( sum(dbinom(round(p.obs*n), n, p.obs, log=TRUE)) - 
                sum(dbinom(round(p.obs*n), n, p.wst, log=TRUE)) ),
    df   <- sum(p.obs < .5),
    pval <- 1 - pchisq(G2, df)
  )

  z <- list(weak=sum(!violdf$wst), moderate=sum(!violdf$mst),
    strong=sum(!violdf$sst), n.tests=max(chkdf$triple),
    wst.violations=wv, mst.violations=mv, sst.violations=sv,
    pcm=R, ranking=idx, chkdf=chkdf, violdf=violdf, wst.fit=wst.fit,
    wst.mat=wst.mat)
  class(z) <- "strans"
  z
}


print.strans <- function(x, digits = max(3, getOption("digits") - 4), ...){
  # Last mod: 21/Aug/2007: replace printCoefmat by print
  #           14/Feb/2009: checks if ?st.violations is empty

  cat("\nStochastic Transitivity\n\n")
  tran <- c(x$weak, x$moderate, x$strong)
  ntst <- x$n.tests
  ttab <- cbind(tran/ntst,
    c(if(length(x$wst.violations)) mean(x$wst.violations) else 0,
      if(length(x$mst.violations)) mean(x$mst.violations) else 0,
      if(length(x$sst.violations)) mean(x$sst.violations) else 0),
    c(if(length(x$wst.violations))  max(x$wst.violations) else 0,
      if(length(x$mst.violations))  max(x$mst.violations) else 0,
      if(length(x$sst.violations))  max(x$sst.violations) else 0),
    c(x$wst.fit[1], NA, NA),
    c(x$wst.fit[2], NA, NA),
    c(x$wst.fit[3], NA, NA))
 
  ttran <- cbind(tran, ttab)
  rownames(ttran) <- c("Weak", "Moderate", "Strong")
  colnames(ttran) <- c("Violations", "ErrorRatio", "MeanDev", "MaxDev",
    "Deviance", "Df", "Pr(>Chi)")

  print(ttran, digits=digits)
  cat("---\nNumber of Tests:", ntst, "\n")
  cat("\n")
  invisible(x)
}

