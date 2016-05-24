#' MPI Demonstrations
#' 
#' These functions are examples of simple statistics via MPI calls.
#' 
#' \code{x.gbd} and \code{y.gbd} are vectors with length \code{N.gbd}.
#' \code{X.gbd} is a matrix with dimension \code{N.gbd * p} and exists on all
#' processors. \code{N.gbd} may be vary across processors.
#' 
#' For demonstration purpose, these objects should not contains weird values
#' such \code{NA}.
#' 
#' @param x.gbd
#' gbd a GBD vector.
#' @param breaks 
#' a set to break data in groups.
#' @param prob 
#' a desired probability for quantile.
#' @param y.gbd 
#' a GBD vector.
#' @param X.gbd 
#' a GBD matrix.
#' 
#' @return 
#' \code{mpi.stat} returns sample mean and sample variance.
#' \code{mpi.bin} returns binning counts for the given breaks.
#' \code{mpi.quantile} returns a quantile.
#' \code{mpi.ols} returns ordinary least square estimates (beta_hat).
#' 
#' @examples
#' \dontrun{
#' ### Under command mode, run the demo with 4 processors by
#' ### (Use Rscript.exe for windows system)
#' mpiexec -np 4 Rscript -e "demo(sample_stat,'pbdDEMO',ask=F,echo=F)"
#' mpiexec -np 4 Rscript -e "demo(binning,'pbdDEMO',ask=F,echo=F)"
#' mpiexec -np 4 Rscript -e "demo(quantile,'pbdDEMO',ask=F,echo=F)"
#' mpiexec -np 4 Rscript -e "demo(ols,'pbdDEMO',ask=F,echo=F)"
#' mpiexec -np 4 Rscript -e "demo(gbd2dmat,'pbdDEMO',ask=F,echo=F)"
#' mpiexec -np 4 Rscript -e "demo(balance,'pbdDEMO',ask=F,echo=F)"
#' }
#' 
#' @keywords programming
#' @name mpi_example
#' @rdname mpi_example
NULL



#' @rdname mpi_example
#' @export
mpi.stat <- function(x.gbd){
  ### For mean(x).
  N <- allreduce(length(x.gbd), op = "sum")
  bar.x.gbd <- sum(x.gbd / N)
  bar.x <- allreduce(bar.x.gbd, op = "sum")

  ### For var(x).
  s.x.gbd <- sum(x.gbd^2 / (N - 1))
  s.x <- allreduce(s.x.gbd, op = "sum") - bar.x^2 * (N / (N - 1))

  list(mean = bar.x, s = s.x)
} # End of mpi.stat().


#' @rdname mpi_example
#' @export
mpi.bin <- function(x.gbd, breaks = pi / 3 * (-3:3)){
  bin.gbd <- table(cut(x.gbd, breaks = breaks))
  bin <- as.array(allreduce(bin.gbd, op = "sum"))
  dimnames(bin) <- dimnames(bin.gbd)
  class(bin) <- class(bin.gbd)
  bin
} # End of mpi.bin().


#' @rdname mpi_example
#' @export
mpi.quantile <- function(x.gbd, prob = 0.5){
  if(sum(prob < 0 | prob > 1) > 0){ 
    stop("prob should be in (0, 1)")
  }

  N <- allreduce(length(x.gbd), op = "sum")
  x.max <- allreduce(max(x.gbd), op = "max")
  x.min <- allreduce(min(x.gbd), op = "min")

  f.quantile <- function(x, prob = 0.5){
    allreduce(sum(x.gbd <= x), op = "sum") / N - prob
  }

  uniroot(f.quantile, c(x.min, x.max), prob = prob[1])$root
} # End of mpi.quantile().


#' @rdname mpi_example
#' @export
mpi.ols <- function(y.gbd, X.gbd){
  if(length(y.gbd) != nrow(X.gbd)){
    stop("length(y.gbd) != nrow(X.gbd)")
  }
  
  t.X.gbd <- t(X.gbd)
  A <- allreduce(t.X.gbd %*% X.gbd, op = "sum")
  B <- allreduce(t.X.gbd %*% y.gbd, op = "sum")

  solve(matrix(A, ncol = ncol(X.gbd))) %*% B
} # End of mpi.ols().

