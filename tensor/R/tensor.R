"tensor" <-
function(A, B, alongA = integer(0), alongB = integer(0))
{
  A <- as.array(A)
  dimA <- dim(A)
  dnA <- dimnames(A)
  if (nnA <- is.null(dnA))
    dnA <- rep(list(NULL), length(dimA))

  B <- as.array(B)
  dimB <- dim(B)
  dnB <- dimnames(B)
  if (nnB <- is.null(dnB))
    dnB <- rep(list(NULL), length(dimB))

  if (length(alongA) != length(alongB))
    stop("\"along\" vectors must be same length")

  # special case of both length zero

  if (length(alongA) == 0) {
    R <- as.vector(A) %*% t(as.vector(B))
    dim(R) <- c(dimA, dimB)
    if (!(nnA && nnB))
    dimnames(R) <- c(dnA, dnB)
    return(R)
  }

  mtch <- dimA[alongA] == dimB[alongB]
  if (any(is.na(mtch)) || !all(mtch))
    stop("Mismatch in \"along\" dimensions")

  seqA <- seq(along=dimA)
  allA <- length(seqA) == length(alongA)
  permA <- c(seqA[-alongA], alongA)
  if (!all(seqA == permA))
    A <- aperm(A, permA)
  dim(A) <- c(
    if (allA) 1 else prod(dimA[-alongA]),
    prod(dimA[alongA])
  )

  seqB <- seq(along=dimB)
  allB <- length(seqB) == length(alongB)
  permB <- c(alongB, seqB[-alongB])
  if (!all(seqB == permB))
    B <- aperm(B, permB)
  dim(B) <- c(
    prod(dimB[alongB]),
    if (allB) 1 else prod(dimB[-alongB])
  )

  R <- A %*% B

  if (allA && allB)
    R <- drop(R)
  else {
    dim(R) <- c(
      if (allA) integer(0) else dimA[-alongA],
      if (allB) integer(0) else dimB[-alongB]
    )
    if (!(nnA && nnB))
      dimnames(R) <- c(dnA[-alongA], dnB[-alongB])
  }
  R
}

"%*t%" <- function(x, y) tensor(x, y, 2, 2)

"%t*%" <- function(x, y) tensor(x, y, 1, 1)

"%t*t%" <- function(x, y) tensor(x, y, 1, 2)
