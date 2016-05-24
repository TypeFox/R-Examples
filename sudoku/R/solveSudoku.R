solveSudoku <- function(z, verbose=FALSE, map=c(1:9,letters), level=0,
                        print.it=TRUE) {
  if (length(z)==1) z <- readSudoku(z, map)
  N <- nrow(z);  Ns <- sqrt(N)                              # Traditionally N=9
  oldngot <- sum(z > 0)
  if (verbose) cat("Known:", oldngot, ", level:", level, "\n")
  isok <- TRUE
  a <- array(NA, c(N,N,N))                       # T=this num, F=not this, NA=?
  fill <- function(i, j, k, txt="") {
    if (length(i)!=1 || length(j)!=1 || length(k)!=1) {isok<<-FALSE; return()}
    if (verbose && txt != "") cat(i, j, "=", k, txt, "\n")
    z[i,j] <<- k
    ain <- a
    a[ i, j,  ] <<- seq(1:N)==k
    a[ i,-j, k] <<- FALSE                         # No other k's in this row
    a[-i, j, k] <<- FALSE                         # No other k's in this column
    for (ii in Ns*((i-1) %/% Ns) + 1:Ns) for (jj in Ns*((j-1) %/% Ns) + 1:Ns)
      if (!(ii==i && jj==j)) a[ii,jj, k] <- FALSE # No other k's in this box
    if (any(a != ain, na.rm=TRUE)) isok <<- FALSE   # You turned a T into an F!
  }
  for (i in 1:N) for (j in 1:N) if (k <- z[i,j]) fill(i, j, k)

  repeat {
    for (i in 1:N) for (j in 1:N)      # Check each cell for only 1 possibility
      if (sum(!a[i,j, ], na.rm=TRUE)==N-1 & sum(a[i,j, ], na.rm=TRUE)==0)
        fill(i, j, which(is.na(a[i,j, ])), "by elimination")
    for (k in 1:N) {                # Now explore each digit (a[ , ,k]) in turn
      for (i in which(rowSums(!a[ , ,k],TRUE)==N-1 & !rowSums(a[ , ,k],TRUE)))
        fill(i, which(is.na(a[i, ,k])), k, "each row has a k")
      for (j in which(colSums(!a[ , ,k],TRUE)==N-1 & !colSums(a[ , ,k],TRUE)))
        fill(which(is.na(a[ ,j,k])), j, k, "each col has a k")
      for (bi in seq(0, N-Ns, Ns)) for (bj in seq(0, N-Ns, Ns)) {
        idx <- cbind(data.matrix(expand.grid(bi + 1:Ns, bj + 1:Ns)), k)
        if (sum(!a[idx], na.rm=TRUE)==N-1 && sum(a[idx], na.rm=TRUE)==0) {
          m <- which(is.na(a[idx]))
          fill(idx[m,1], idx[m,2], idx[m,3], "each box has a k")
        }
      }
    }

    if (!isok) {if (verbose) cat("Inconsistent level", level, "\n"); return()}
    ngot <- sum(z > 0)
    if (verbose) cat("Known:", ngot, ", level:", level, "\n\n")
    if (ngot==N^2) {
      if (print.it) print(matrix(map[z],N), quote=FALSE, right=TRUE)
      return(invisible(z))
    }
    if (ngot==oldngot) {                               # Failed.  Take a guess!
      poss <- rowSums(is.na(a), ,2)                # Number of possible guesses
      if (!any(poss > 0)) {if (verbose) cat("No possibilities left\n"); return()}
      ij <- which(poss == min(setdiff(poss,0)), TRUE)[1, ]
      k <- which(is.na(a[ij[1], ij[2], ]))[1]              # 1st possible guess
      if (verbose) cat("Guessing:", ij[1], ij[2], "=", k, "\n")
      zg <- z
      zg[ij[1], ij[2]] <- k
      res <- Recall(zg, verbose, map, level+1, print.it)
      if (is.null(res)) a[ij[1], ij[2], k] <- FALSE else return(invisible(res))
    }
    oldngot <- ngot
  }
}
