.listMat <- function(M, f) {
  l <- vector("list", ncol(M))
  names(l) <- dimnames(M)[[2]]
  for(i in 1:ncol(M))
      l[[i]] <- f(M[, i])
  l
}

.paste0 <- function(...) paste(..., sep = "")

.zf <- function(x) {
  x <- as.character(x)
  m <- max(n <- nchar(x))
  z <- paste(rep(0, m), collapse="")
  .paste0(substring(z, 0, m-n), x)
}

conf.design <- function(G, p, block.name = "Blocks",
                        treatment.names = NULL) {
  if(!is.matrix(G))
        G <- rbind(G)
  if(is.null(treatment.names))
      treatment.names <- if(is.null(nam <- dimnames(G)[[2]]))
          .paste0("T", .zf(1:ncol(G))) else nam

  stopifnot(is.character(treatment.names), length(treatment.names) == ncol(G),
            is.numeric(G), all(G >= 0), all(G %% 1 == 0),
            is.numeric(p), length(p) == 1, p > 0, p %% 1 ==0,
            p %in% primes(p),
            is.character(block.name), length(block.name) == 1)

  D <- as.matrix(expand.grid(rep(list(0:(p-1)), ncol(G))))
  B <- .listMat((D %*% t(G)) %% p, format)
  o <- do.call(order, B)
  B <- do.call(".paste0", B)
  D <- cbind(B, format(D))[o,  ]
  dimnames(D) <- list(B[o], c(block.name, treatment.names))
  data.frame(.listMat(D, factor))
}

.space <- function(G, p) {
###
### Generate all distinct linear combinations of the rows of G over GF(p)
###
  x <- 0:(p - 1)
  M <- as.matrix(x)
  k <- nrow(G)
  if(k > 1) {
    for(i in 2:k) {
      N <- NULL
      for(j in x)
          N <- rbind(N, cbind(M, j))
      M <- N
    }
  }
  M <- (M %*% G) %% p
###
### if the rows of G can be assumed linearly independent the rest can
### be omitted.
###
  m <- 0
  for(j in 1:ncol(M))
      m <- p * m + M[, j]
  M[!duplicated(m),  , drop = FALSE]
}

conf.set <- function(G, p) {
  S <- .space(G, p)
  dimnames(S)[[2]] <- dimnames(G)[[2]]
  S[apply(S, 1, function(x)
          any(t <- x > 0) && x[t][1] == 1),  ]
}

direct.sum <- function(D1, ..., tiebreak = letters) {
  l <- list(...)
  if(length(l) == 0) return(D1)
  D2 <- if(length(l) == 1) l[[1]] else {
    Recall(..., tiebreak = c(tiebreak[-1], tiebreak[1]))
  }
  stopifnot(is.data.frame(D1), is.data.frame(D2))
  n1 <- nrow(D1)
  n2 <- nrow(D2)
  D1 <- D1[rep(seq_len(n1), each  = n2), ]
  D2 <- D2[rep(seq_len(n2), times = n1), ]
  if(any(dups <- names(D2) %in% names(D1)))
      names(D2)[dups] <- .paste0(names(D2)[dups], tiebreak[1])
  D <- cbind(D1, D2)
  row.names(D) <- format(seq_len(nrow(D)))
  D
}

factorize <- function(x, ...)
  UseMethod("factorize")

factorize.default <- function(x, divisors = primes(max(x)), ...) {
  stopifnot(is.numeric(x))
  if (length(x) > 1) {
    l <- vector("list", length(x))
    names(l) <- as.character(x)
    for (i in seq_along(x))
      l[[i]] <- Recall(x[i], divisors = divisors, ...)
    return(l)
  }
  if (x %% 1 > 0 || x < 2)
    return(x)
  tab <- divisors
  fac <- numeric(0)
  while(length(tab <- tab[x %% tab == 0]) > 0) {
    x <- x/prod(tab)
    fac <- c(fac, tab)
  }
  sort(fac)
}

factorize.factor <- function(x, name = deparse(substitute(x)),
                             extension = letters, ...) {
  llev <- factorize.default(length(levels(x)))
  if(length(llev) == 1)
    return(structure(data.frame(x), names = name))

  D <- expand.grid(lapply(llev, seq_len))[x, ] - 1

  D <- lapply(D, factor)
  names(D) <- .paste0(name, rep(extension, length.out = length(D)))
  D <- data.frame(D)
  row.names(D) <- format(seq_len(nrow(D)))
  D
}

.makeList <- function(x) if(is.list(x)) x else list(x)

join <- function(...) {
  l <- list(...)
  l <- do.call(c, lapply(l, .makeList))
  l <- lapply(l, function(f)
              format(as.character(f)))
  factor(do.call(paste, c(l, list(sep = ":"))))
}

rjoin <- function(..., part.name = "Part") {
  l <- lapply(list(...), as.data.frame)	### for some safety...
  checkNames <- all(sapply(l, function(bit)
                           all(names(bit) == names(l[[1]]))))
  if(!checkNames) stop("The column names differ between components!")
  bf <- factor(.paste0(part.name, .zf(rep(1:length(l), sapply(l, nrow)))))
  D <- data.frame(Part = bf, do.call("rbind", l))
  names(D)[1] <- part.name
  D
}

primes <- function(n) {
### Find all primes less than n (or max(n) if length(n) > 1).
### Uses an obvious sieve method.  Nothing flash.
###
### 2013: This now uses a slightly improved coding of the version of
###       the algorithm used in the pracma package primes() function
###
  stopifnot(is.numeric(n), all(n %% 1 == 0))
  if ((M2 <- max(n)) <= 1)
    return(numeric(0))
  x <- seq(1, M2, 2)
  np <- length(x)
  x[1] <- 2
  if(M2 > 8) {
    top <- floor(sqrt(M2))
    p <- 1
    while((p <- p+2) <= top)
        if(x[(p + 1)/2] > 0)
            x[seq((p*p + 1)/2, np, p)] <- 0
  }
  x[x > 0]
}



