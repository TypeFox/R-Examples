### This file contains the function for the gaps plot.

plotgaps <- function(X, code.type = .code.type[1],
    main = "Gaps Plot", xlab = "Sites", ylab = "Proportion", ...){
  if(sum(code.type %in% .code.type) != 1){
    stop("The code.type is not found.")
  }

  if(code.type == .code.type[1]){
    GAP <- .nucleotide$nid[.nucleotide$code == "-"]
  } else if(code.type == .code.type[2]){
    GAP <- .snp$sid[.snp$code == "-"]
  } else{
    stop("code.type is not implemented.")
  }
  my.col <- "gray"

  N <- nrow(X)
  L <- ncol(X)
  xlim <- c(1, L)
  ylim <- c(0, 1)

  prop <- colSums(X == GAP) / N

  plot(1:L, prop, type = "l",
       xlim = xlim, ylim = ylim, col = my.col,
       main = main, xlab = xlab, ylab = ylab)
} # End of plotdots().

