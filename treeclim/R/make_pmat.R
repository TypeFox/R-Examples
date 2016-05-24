##' Reformat climate input data into separate lists
##' 
##' Make separate lists for each parameter in the climate data input, and make
##' them accessible by the parameter names. The single lists correspond to the
##' scheme years in rows, observations from previous december to current
##' december in columns.
##' @title Make list of parameter matrices from climate input data
##' @param x climate data as returned from as_tcclimate
##' @param pad should the previous year be padded with NAs? This is useful when
##'   we want to use the full range of available data when all predictors derive
##'   from the current year
##' @return a list of matrices
##' @keywords manip internal
make_pmat <- function(x, pad = FALSE) {
  if (any(names(x) == "year")) {
    years <- x$year
  } else {
    years <- x[,1]
  }
  n <- length(unique(years))
  no.vars <- dim(x)[2] - 2
  .names <- names(x)[-c(1:2)]

  ## create a full matrix from -1 to 12 for every parameter, names are
  ## not needed for this step

  m <- list()
  
  if (pad) {
    npad <- n + 1
    padyears <- c(min(years) - 1, unique(years))
    x <- rbind(x[1:12,], x)
    x[1:12, 1] <- x[1:12, 1] - 1
    x[1:12, 3:(2 + no.vars)] <- NA
    years <- x[,1]
  } else {
    npad <- n 
    padyears <- unique(years)
  }

  for (k in 1:no.vars) {
    m_k <- matrix(NA, nrow = 24, ncol = npad - 1)
    colnames(m_k) <- padyears[-1]
    for (i in 2:npad) {
      start_with <- which(years == padyears[i - 1])[1] # previous january
      for (j in 1:24) {                 # loop through months
        m_k[j, (i - 1)] <- x[(start_with + j - 1), 2 + k]
      }
    }
    m[[k]] <- t(m_k)
  }

  names(m) <- .names
  class(m) <- list("ctpmat", "list")
  m
}
