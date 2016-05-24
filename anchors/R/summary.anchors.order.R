#######################################################################
##
## Function: summary.anchors.rank() 
## Author  : Jonathan Wand <wand(at)stanford.edu>
## Created : 2006-10-02
##
## DESCRIPTION : Print frequency of vignette ordering
## 
## MODIFIED:  
## 
#######################################################################
summary.anchors.order <- function(object, top, digits=getOption("digits"),verbose=FALSE,...) {
  if (missing(top)) { top <- length(object$freq) }

  cat("ANCHORS: SUMMARY OF VIGNETTE ORDERING\n\n")

  cat(object$main,"\n\n")

  N <- object$N

  if (object$ties=="set") {

    y0 <- xtabs( ~ object$labels )
    r0 <- rownames(y0)
    
    y1 <- xtabs( ~ object$labels + object$ngroup)
    m1 <- as.numeric(colnames(y1))[max.col(y1)]
    r1 <- rownames(y1)

    y2 <- xtabs( ~ object$labels + object$nviol)
    m2 <- as.numeric(colnames(y2))[max.col(y2)]
    r2 <- rownames(y2)

    if ( any(r0 != r1) || any(r0 != r2) ) stop("mismatch in summary\n")
    
    z <- as.data.frame(list(Frequency  = as.numeric(y0),
                            Proportion = as.numeric(y0) / N,
                            Ndistinct  = m1,
                            Nviolation = m2
                            )
                       )
    rownames(z) <- r0
    z <- z[ order( as.numeric(y0), decreasing=TRUE), ]    
  } else {
    z <- as.data.frame(list(Frequency  = object$freq,
                            Proportion = object$freq / N))
  }

    cat("Number of cases with at least two distinct vignette responses:",
        sum(object$ngroup > 1),"\n")
    cat("and with no violations of natural ordering:",
        sum(object$ngroup > 1 & object$nviol == 0),"\n")
    cat("and with no more than 1 violation of natural ordering:",
        sum(object$ngroup > 1 & object$nviol <= 1),"\n")
    cat("and with no more than 2 violation of natural ordering:",
        sum(object$ngroup > 1 & object$nviol <= 2),"\n")
    cat("\n")

    if (verbose) {
      cat("\n")
      Nviolation <- as.character(object$nviol)
      Nviolation[ object$nviol > 4 ] <- "5+"
      Ngroup <- object$ngroup
      print(xtabs(~ Ngroup + Nviolation))
      cat("\n")
    }

      cat("Proportion of cases a vignette (row) is less than another (column):\n")
      print(round( object$compare.matrix, digits))
      cat("\n")

      cat("Upper tri =     p_{ij} - p_{ji} (negative values suggest misorderings)\n")
      cat("Lower tri = 1 - p_{ij} - p_{ji} (big numbers means many ties)\n")
      tt <- idx <- upper.tri(object$compare.matrix)
      lidx <- lower.tri(object$compare.matrix)                   
      tt[idx] <- object$compare.matrix[idx] - t(object$compare.matrix)[idx]
      tt[lidx] <- 1 - object$compare.matrix[idx] - t(object$compare.matrix)[idx]
      diag(tt) <- NA
      colnames(tt) <- rownames(object$compare.matrix)
      rownames(tt) <- rownames(object$compare.matrix)
      print(round(tt,digits))
      cat("\n")

  
  cat("Top",top,"orderings (out of",NROW(z),"unique orderings):\n\n")
  print(format(z[1:top,],digits=digits))
  cat("\n")
  
  return(invisible(z))
}

