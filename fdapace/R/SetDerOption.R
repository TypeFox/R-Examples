SetDerOptions <- function(fpcaObject = NULL, derOptns = list()) {
  if (is.null(derOptns)) {
    derOptns <- list()
  }
  # These are relevant for fitted.FPCA
  derOptns$method <- ifelse (is.null(derOptns$method), 'EIG',
                            derOptns$method)
  #derOptns$k <- ifelse (is.null(derOptns$k), length(fpcaObject$lambda), derOptns$k)
  # derOptns$GCV <- ifelse (is.null(derOptns$GCV), FALSE, TRUE)
  
  derOptns$p <- ifelse (is.null(derOptns$p), 0, derOptns$p)
  derOptns$kernelType <-  ifelse(is.null(derOptns$kernelType), 'gauss',
                                 derOptns$kernelType)
  derOptns$bw <- ifelse( is.null(derOptns$bw), 
                         ifelse( !is.null(fpcaObject$sigma2) && (fpcaObject$sigma2 / sum(fpcaObject$lambda)) >= 0.01, 
                                  derOptns$p * 0.10 * diff(range(fpcaObject$workGrid)), derOptns$p * 0.05 * diff(range(fpcaObject$workGrid))),
                         derOptns$bw) 
  return(derOptns)
}
