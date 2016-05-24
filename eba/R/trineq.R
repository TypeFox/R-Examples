# Trinary inequality
# Last mod: May/26/2013, FW


## Check inclusion rule (Tversky & Sattath, 1979, p. 546)
inclusion.rule <- function(A){
  inclusion <- logical(length(A) * choose(length(A) - 1, 2))
  i <- 1
  idx <- seq_along(A)
  for (x in idx){                             # for each stimulus combined
    for (y in idx[-c(x, length(idx))]){       #   with all pairs 
      for (z in idx[idx != x & idx > y]){     #   of remaining stimuli
        inclusion[i] <- 
          all(intersect(A[[x]], A[[y]]) %in% intersect(A[[x]], A[[z]])) |
          all(intersect(A[[x]], A[[z]]) %in% intersect(A[[x]], A[[y]]))
        i <- i + 1
      }
    }
  }
  all(inclusion)
}


## Check trinary inequality (Tversky & Sattath, 1979, p. 554)
trineq <- function(M, A = 1:I){
  I <- sqrt(length(as.matrix(M)))    # number of stimuli
  if (!inclusion.rule(A)) warning("family of aspect sets is not a tree")

  R <- as.matrix( M / (M + t(M)) )   # pcm rel. freq.
  obj.names <- if(!is.null(colnames(M))) colnames(M)
               else make.unique(rep(letters, length.out=I), sep="")

  Pxy <- Rxyz <- triple <- NULL
  idx <- seq_len(I)
  for (x in 1:(I - 1)){
    for (y in (x + 1):I){
      for (z in idx[!idx %in% c(x, y)]){
        if (length(unlist(A)) == I){  # BTL model
          Rxyz <- c(Rxyz,
                    R[x, y]/R[y, x] * R[y, z]/R[z, y] * R[z, x]/R[x, z])
          triple <- c(triple, paste0("(", obj.names[x], ",", obj.names[y],
                      ")", obj.names[z]))
          Pxy <- c(Pxy, R[x, y])
        } else {                      # preference tree
          x_and_y <- intersect(A[[x]], A[[y]])
          x_and_z <- intersect(A[[x]], A[[z]])
          y_and_z <- intersect(A[[y]], A[[z]])
          if (all(x_and_z == y_and_z) &
              all(x_and_z %in% x_and_y) &
              length(x_and_z) < length(x_and_y) &
              R[x, y] != 1/2){
            if (R[x, y] > 1/2){
              Rxyz <- c(Rxyz,
                        R[x, y]/R[y, x] * R[y, z]/R[z, y] * R[z, x]/R[x, z])
              triple <- c(triple, paste0("(", obj.names[x], ",", obj.names[y],
                          ")", obj.names[z]))
              Pxy <- c(Pxy, R[x, y])
            } else {
              Rxyz <- c(Rxyz,
                        R[y, x]/R[x, y] * R[x, z]/R[z, x] * R[z, y]/R[y, z])
              triple <- c(triple, paste0("(", obj.names[y], ",", obj.names[x],
                          ")", obj.names[z]))
              Pxy <- c(Pxy, R[y, x])
            }
          }
        }
      }
    }
  }
  chkdf <- data.frame(triple, Pxy, Rxyz, trineq=Rxyz > 1)
  chkdf <- chkdf[Rxyz != 0 & is.finite(chkdf$Rxyz),]  # exclude 0/1 probs
  rownames(chkdf) <- NULL
  out <- list(n=nrow(chkdf), prop=mean(chkdf$Rxyz > 1),
              quant=quantile(chkdf$Rxyz, prob=c(.25, .5, .75)),
              chkdf=chkdf)
  class(out) <- "trineq"
  out
}


print.trineq <- function(x, digits = max(3, getOption("digits") - 4), ...){
  cat("\nTrinary Inequality\n")
  cat("\nNumber of Tests:", x$n, "\n")
  cat("% triples confirming trinary inequality:",
      format(x$prop, digits=digits, ...), "\n")
  cat("Quantiles of R(xyz):\n")
  print(x$quant, digits=digits, ...)
  cat("\n")
  invisible(x)
}

