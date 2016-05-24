export.qm <- function(qmobject, file, style = c("R", "PQMethod")) {
  style <- match.arg(style)
  if (style=="R") capture.output(print(qmobject, length=max(qmobject$brief$nstat, qmobject$brief$nqsorts)), file=file)
  else if (style=="PQMethod") {
    pqmethod.output <- function(qmobject) {
      cat(qmobject$brief$info, sep="\n")
      cat("\n")
      pqout <- as.list(rep(NA, 15))
      # Correlation Matrix Between Sorts
      cat("\n\nCorrelation Matrix Between Sorts", sep="\n")
      print(round(cor(qmobject[[2]]), digits=4))
      # Unrotated Factor Matrix
      cat("\n\nUnrotated Factor Matrix", sep="\n")
      print(principal(qmobject[[2]], nfactors=8, rotate="none")$loadings, cutoff=0)
      # Cumulative Communalities Matrix
      cat("\n\nCumulative Communalities Matrix", sep="\n")
      comm <- data.frame(matrix(as.numeric(NA), ncol=8, 
                                nrow=length(qmobject[[2]])))
      for (i in 1:8) {
        comm[[i]] <- principal(qmobject[[2]], nfactors=i, 
                               rotate="none")$communality
      }
      print(round(comm, digits=4))
      # Factor Matrix with an X Indicating a Defining Sort
      cat("\n\nFactor Matrix and Defining Sorts", sep="\n")
      fmx <- list(round(qmobject$loa, digits=4),       qmobject$flagged)
      names(fmx) <- c("Q-sort factor loadings", "Flagged Q-sorts")
      print(fmx)
      # Free Distribution Data Results
      cat("\n\nFree Distribution Data Results -- not calculated", sep="\n")
      # Factor Scores with Corresponding Ranks
      cat("\n\nFactor Scores (z-scores)", sep="\n")
      print(round(qmobject$zsc, digits=2))
      # Correlations Between Factor Scores
      cat("\n\nCorrelations Between Factor Scores", sep="\n")
      print(round(qmobject$f_char[[2]], digits=4))
      # Factor Scores -- For Factor *
      cat("\n\nFactor Scores ", sep="\n")
      nfactors <- length(qmobject$zsc_n)
      fsco <- as.list(1:nfactors)
      for (i in 1:nfactors) {
        names(fsco)[i] <- paste0("-- For Factor ", i)
        fsco[[i]] <- round(qmobject$zsc[order(qmobject$zsc[i], decreasing=T), c(i, i)][1], digits=3)
      }
      print(fsco)
      # Descending Array of Differences Between Factors x and y
      cat("\n\nDescending Array of Differences Between Factors ", sep="\n")
      comparisons <- combn(nfactors, 2, simplify=F)
      daf <- as.list(1:length(comparisons))
      for (i in 1:length(comparisons)) {
        names(daf)[i] <- paste(comparisons[[i]], collapse=" and ", sep="")
        zsc <- qmobject$zsc[comparisons[[i]]]
        dif <- qmobject[[8]][c(names(qmobject[[8]])[grep(paste(paste(comparisons[[i]][1], comparisons[[i]][2], sep=".*"), paste(comparisons[[i]][2], comparisons[[i]][1], sep=".*"), sep="|"),names(qmobject[[8]]))],"dist.and.cons")]
        dad <- cbind(zsc, dif)
        daf[[i]] <- format(dad[order(dad[3], decreasing = T), ], digits=2)
      }
      print(daf)
      # Factor Q-Sort Values for Each Statement
      cat("\n\nFactor Q-Sort Values for Each Statement", sep="\n")
      print(qmobject$zsc_n)
      # Factor Q-Sort Values for Statements sorted by Consensus vs. Disagreement (Variance across Factor Z-Scores)
      cat("\n\nFactor Q-Sort Values for Statements sorted by Consensus vs. Disagreement (Variance across Factor Z-Scores)", sep="\n")
      zsc <- qmobject$zsc
      zsc.ord <- order(apply(zsc, 1, var))
      print(qmobject$zsc_n[zsc.ord,])
      # Factor Characteristics
      cat("\n\nFactor Characteristics", sep="\n")
      fch <- t(qmobject$f_char[[1]])
      rownames(fch) <- c(
        "Average reliability coefficient", 
        "Number of loading Q-sorts", 
        "Eigenvalues", 
        "Percentage of explained variance", 
        "Composite reliability", 
        "Standard error of factor scores")
      print(round(fch, digits=2))
      # Standard Errors for Differences in Factor Z-Scores
      cat("\n\nStandard Errors for Differences in Factor Z-Scores", sep="\n")
      print(qmobject$f_char[[3]])
      # Distinguishing Statements for Factor *
      cat("\n\nDistinguishing Statements ", sep="\n")
      dc <- qmobject$qdc
      zsc <- qmobject$zsc
      dsf <- as.list(1:nfactors)
      for (i in 1:nfactors) {
        names(dsf)[i] <- paste0("for Factor ", i)
        d <- grep(paste0("f",i, "|all"), dc$dist.and.cons)
        dsf[[i]] <- cbind(round(zsc[d, ], digits=2),dc[d, c(1,1+(2*(1:length(comparisons))))])
      }
      print(dsf)
      # Consensus Statements  --  Those That Do Not Distinguish Between ANY Pair of Factors.
      cat("\n\nConsensus Statements  --  Those That Do Not Distinguish Between ANY Pair of Factors.", sep="\n")
      dcon <- which(dc$dist.and.cons == "Consensus")
      print(cbind(round(qmobject$zsc[dcon, ], digits=2),dc[dcon, c(1,1+(2*(1:length(comparisons))))]))
      #return(pqout)
    }
    capture.output(pqmethod.output(qmobject), file=file)
  }
}