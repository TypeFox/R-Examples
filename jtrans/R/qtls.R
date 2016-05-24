#' @rdname fit_functions
qtls <- function(x, z) {
  qtls <- quantile(x, probs=pnorm(c(-3 * z, -z, z, 3 * z)))
  q <- list(xl = qtls[2] - qtls[1],
            xm = qtls[3] - qtls[2],
            xu = qtls[4] - qtls[3],
            QR = (qtls[4] - qtls[3]) * (qtls[2] - qtls[1]) / 
              (qtls[3] - qtls[2])^2, 
            z  = z,
            x2 = qtls[2], 
            x3 = qtls[3])
  return(q)
}

