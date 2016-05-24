## Calls Fortran implentation of:
## R.E. BURKARD and F. RENDL. A thermodynamically motivated
## simulation procedure for combinatorial optimization problems.
## European Journal of Operations Research, 17(2):169-174, 1984.

qapSA <- function(A, B, rep = 1L, miter = 2*nrow(A), fiter = 1.1, ft = .5,
   maxsteps = 50L, verbose = FALSE) {
  A <- unname(as.matrix(A))
  B <- unname(as.matrix(B))

  storage.mode(A) <- "double"
  storage.mode(B) <- "double"
  n <- nrow(A)
  if(any(dim(A) != n) || any(dim(B) != n)) stop("Matrix do not conform!")

  if(!isSymmetric(A) || !isSymmetric(B))
    stop("Heuristic only available for symmetric QAP.")
  if(any(A<0) || any(B<0))
    stop("All values in A and B need to be positive.")

  if(ft<0 || ft>=1) stop("ft needs to be in (0 ,1).")

  if(verbose) cat("Simulated Annealing Heuristic by Burkart and Rendl.\n")
  if(verbose) cat(sprintf("%5s %10s %10s\n", "rep", "best_obj", "current_obj"))

  best_perm <- NULL
  best_obj <- Inf

  ## we do repetitions in R (not in the FORTRAN code)
  ## we start with a random permutation in perm
  for(i in 1:rep) {
    res <- .Fortran("qaph4", n = n, a = A, b = B,
      miter = as.integer(miter), fiter = as.double(fiter),
      ft = as.double(ft), ope = integer(n), ol = double(1), perm = sample(n),
      maxsteps = as.integer(maxsteps), PACKAGE = "qap")


    if(res$ol < best_obj) {
      best_obj <- res$ol
      best_perm <- res$ope
    }

    if(verbose) cat(sprintf("%5i %10.0f %10.0f\n", i, best_obj, res$ol))

  }

  attr(best_perm, "obj") <- best_obj
  best_perm
}
