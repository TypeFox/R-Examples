kseq <- function(from = NULL, to = NULL, f = 0.05, ...)
{
  if (is.null(from)) {
    stop(sQuote("from"), " must not be NULL!", call. = FALSE)
  }

  if (!is.wholenumber(from)) {
    stop(sQuote("from"), " integer is requested!", call. = FALSE)
  }
  
  if (is.null(to)) {
    stop(sQuote("to"), " must not be NULL!", call. = FALSE)
  }

  if (!is.wholenumber(to)) {
    stop(sQuote("to"), " integer is requested!", call. = FALSE)
  }  
  
  if (to <= from) {
    stop(sQuote("to"), " must be greater than ", sQuote("from"), "!", call. = FALSE)
  }  
  
  if (from < 1) {
    stop(sQuote("from"), " must be greater than 0!", call. = FALSE)
  }    
  
  if (!is.numeric(f)) {
    stop(sQuote("f"), " numeric is requested!", call. = FALSE)
  }

  if ((f <= 0.0) || (f >= 1.0)) {
    stop(sQuote("f"), " must be greater than 0.0 and less than 1.0!", call. = FALSE)
  }  

  K <- array(0)

  i <- 1; K[i] <- from

  while (to > K[i]) {
    K[i + 1] <- ceiling(K[i] / (1.0 - f)); i <- i + 1
  }
  
  K[i] <- to

  rm(list = ls()[!(ls() %in% c("K"))])
  
  invisible(K)
} ## kseq
