adiag <- function (..., pad = as.integer(0), do.dimnames = TRUE) # function from package 'magic'
{
    args <- list(...)
    if (length(args) == 1) {
        return(args[[1]])
    }
    if (length(args) > 2) {
        jj <- do.call("Recall", c(args[-1], list(pad = pad)))
        return(do.call("Recall", c(list(args[[1]]), list(jj), 
            list(pad = pad))))
    }
    a <- args[[1]]
    b <- args[[2]]
    if (is.null(b)) {
        return(a)
    }
    if (is.null(dim(a)) & is.null(dim(b))) {
        dim(a) <- rep(1, 2)
        dim(b) <- rep(1, 2)
    }
    if (is.null(dim(a)) & length(a) == 1) {
        dim(a) <- rep(1, length(dim(b)))
    }
    if (is.null(dim(b)) & length(b) == 1) {
        dim(b) <- rep(1, length(dim(a)))
    }
    if (length(dim.a <- dim(a)) != length(dim.b <- dim(b))) {
        stop("a and b must have identical number of dimensions")
    }
    s <- array(pad, dim.a + dim.b)
    s <- do.call("[<-", c(list(s), lapply(dim.a, seq_len), list(a)))
    ind <- lapply(seq(dim.b), function(i) seq_len(dim.b[[i]]) + 
        dim.a[[i]])
    out <- do.call("[<-", c(list(s), ind, list(b)))
    n.a <- dimnames(a)
    n.b <- dimnames(b)
    if (do.dimnames & !is.null(n.a) & !is.null(n.b)) {
        dimnames(out) <- mapply(c, n.a, n.b, SIMPLIFY = FALSE)
        names(dimnames(out)) <- names(n.a)
    }
    return(out)
}

vec2mat2 <- function (x, sep = "-") 
{
    splits <- strsplit(x, sep)
    n.spl <- sapply(splits, length)
    if (any(n.spl != 2)) 
        stop("Names must contain exactly one '", sep, "' each;  instead got ", 
            paste(x, collapse = ", "))
    x2 <- t(as.matrix(as.data.frame(splits)))
    dimnames(x2) <- list(x, NULL)
    x2
}

multcompLetters <- function (x, compare = "<", threshold = 0.05,   # function from package 'multcompView'
  Letters = c(letters, LETTERS, "."), reversed = FALSE) 
{
  x.is <- deparse(substitute(x))
  if (class(x) == "dist") 
    x <- as.matrix(x)
  if (!is.logical(x)) 
    x <- do.call(compare, list(x, threshold))
  dimx <- dim(x)
{
    if ((length(dimx) == 2) && (dimx[1] == dimx[2])) {
      Lvls <- dimnames(x)[[1]]
      if (length(Lvls) != dimx[1]) 
        stop("Names requred for ", x.is)
      else {
        x2. <- t(outer(Lvls, Lvls, paste, sep = ""))
        x2.n <- outer(Lvls, Lvls, function(x1, x2) nchar(x2))
        x2.2 <- x2.[lower.tri(x2.)]
        x2.2n <- x2.n[lower.tri(x2.n)]
        x2a <- substring(x2.2, 1, x2.2n)
        x2b <- substring(x2.2, x2.2n + 1)
        x2 <- cbind(x2a, x2b)
        x <- x[lower.tri(x)]
      }
    }
    else {
      namx <- names(x)
      if (length(namx) != length(x)) 
        stop("Names required for ", x.is)
      x2 <- vec2mat2(namx)
      Lvls <- unique(as.vector(x2))
    }
  }
  n <- length(Lvls)
  LetMat <- array(TRUE, dim = c(n, 1), dimnames = list(Lvls, 
                                                       NULL))
  k2 <- sum(x)
  if (k2 == 0) {
    Ltrs <- rep(Letters[1], n)
    names(Ltrs) <- Lvls
    dimnames(LetMat)[[2]] <- Letters[1]
    return(Ltrs)
  }
  distinct.pairs <- x2[x, , drop = FALSE]
  absorb <- function(A.) {
    k. <- dim(A.)[2]
    if (k. > 1) {
      for (i. in 1:(k. - 1)) for (j. in (i. + 1):k.) {
        if (all(A.[A.[, j.], i.])) {
          A. <- A.[, -j., drop = FALSE]
          return(absorb(A.))
        }
        else {
          if (all(A.[A.[, i.], j.])) {
            A. <- A.[, -i., drop = FALSE]
            return(absorb(A.))
          }
        }
      }
    }
    A.
  }
  for (i in 1:k2) {
    dpi <- distinct.pairs[i, ]
    ijCols <- (LetMat[dpi[1], ] & LetMat[dpi[2], ])
    if (any(ijCols)) {
      A1 <- LetMat[, ijCols, drop = FALSE]
      A1[dpi[1], ] <- FALSE
      LetMat[dpi[2], ijCols] <- FALSE
      LetMat <- cbind(LetMat, A1)
      LetMat <- absorb(LetMat)
    }
  }
  sortCols <- function(B) {
    firstRow <- apply(B, 2, function(x) which(x)[1])
    B <- B[, order(firstRow)]
    firstRow <- apply(B, 2, function(x) which(x)[1])
    reps <- (diff(firstRow) == 0)
    if (any(reps)) {
      nrep <- table(which(reps))
      irep <- as.numeric(names(nrep))
      k <- dim(B)[1]
      for (i in irep) {
        i. <- i:(i + nrep[as.character(i)])
        j. <- (firstRow[i] + 1):k
        B[j., i.] <- sortCols(B[j., i., drop = FALSE])
      }
    }
    B
  }
  LetMat. <- sortCols(LetMat)
  if (reversed) 
    LetMat. <- LetMat.[, rev(1:ncol(LetMat.))]
  k.ltrs <- dim(LetMat.)[2]
  makeLtrs <- function(kl, ltrs = Letters) {
    kL <- length(ltrs)
    if (kl < kL) 
      return(ltrs[1:kl])
    ltrecurse <- c(paste(ltrs[kL], ltrs[-kL], sep = ""), 
                   ltrs[kL])
    c(ltrs[-kL], makeLtrs(kl - kL + 1, ltrecurse))
  }
  Ltrs <- makeLtrs(k.ltrs, Letters)
  dimnames(LetMat.)[[2]] <- Ltrs
  LetVec <- rep(NA, n)
  names(LetVec) <- Lvls
  for (i in 1:n) LetVec[i] <- paste(Ltrs[LetMat.[i, ]], collapse = "")
  nch.L <- nchar(Ltrs)
  blk.L <- rep(NA, k.ltrs)
  for (i in 1:k.ltrs) blk.L[i] <- paste(rep(" ", nch.L[i]), 
                                        collapse = "")
  monoVec <- rep(NA, n)
  names(monoVec) <- Lvls
  for (j in 1:n) {
    ch2 <- blk.L
    if (any(LetMat.[j, ])) 
      ch2[LetMat.[j, ]] <- Ltrs[LetMat.[j, ]]
    monoVec[j] <- paste(ch2, collapse = "")
  }
  return(monoVec)
}


