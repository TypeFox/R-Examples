# # @rdname mvspectrum
# # 
# # @description
# # \code{.SDF2mvspectrum} converts the output of \code{\link[sapa]{SDF}} to a 3D array with
# # frequencies in the firs and the spectral density matrix at a given frequency in the latter
# # two dimensions.
# # 
# # @details
# # The \code{\link[sapa]{SDF}} function only returns the upper diagonal matrix
# # (including diagonal), since spectrum matrices are Hermitian. For fast
# # vectorized computations, however, the full matrices are required.
# # Thus \code{.SDF2mvspectrum} converts SDF output to a \eqn{3D} array with
# # number of frequencies in the first dimension and the spectral density matrix
# # in the latter two.
# # 
# # \code{.SDF2mvspectrum} is typically not called by the user, but by
# # \code{\link{mvspectrum}}.
# # 
# # @param sdf.output an object of class \code{"SDF"} from \code{\link[sapa]{SDF}}
# # @keywords ts manip

.SDF2mvspectrum <- function(sdf.output) {
  num.series <- attr(sdf.output, "n.series")
  num.freqs <- length(attr(sdf.output, "frequency"))
  f.lambda <- array(NA, dim = c(num.freqs, num.series, num.series))
  
  if (num.series > 1) {
    col.sels <- seq_len(num.series)
    # TODO: make this faster/vectorized
    for (ii in seq_len(num.series)) {
      f.lambda[, ii, ii:num.series] <- sdf.output[, col.sels]
      col.sels <- col.sels + (num.series - ii)
      col.sels <- col.sels[-1]
    }
    f.lambda <- apply(f.lambda, 1, fill_hermitian)
    f.lambda <- array(t(f.lambda), dim = c(num.freqs, num.series, num.series))
  } else {
    f.lambda[, 1, 1] <- c(sdf.output)
  }
  # remove frequency 0
  f.lambda <- f.lambda[-1, , ]
  attr(f.lambda, "frequency") <- attr(sdf.output, "frequency")[-1] * pi
  
  return(f.lambda)
} 