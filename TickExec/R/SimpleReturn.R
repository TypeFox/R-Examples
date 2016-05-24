#### Calculate return of two given price matrices ####
#### priceTo should be later in time than priceFrom ####

SimpleReturn <- function (priceFrom, priceTo, diff = 0) {
  ## formalize arguments ##
  priceFrom <- as.matrix(priceFrom)
  priceTo   <- as.matrix(priceTo)
  
  ## check dimensions ##
  l = nrow(priceFrom)
  if (l != nrow(priceTo)) {
    stop('Unequal number of rows.')
  }
  if (ncol(priceFrom) != ncol(priceTo)) {
    stop('Unequal number of columns.')
  }
  
  ## calculation ##
  if (diff > 0) {
    if (l <= diff ) {
      stop('not enough rows.')
    }
    mTo <- priceTo[-diff, ]
    mFrom <- priceFrom[1:(l-diff), ]
    prefix <- matrix(ncol = ncol(priceFrom), nrow = diff, NA)
    return (rbind(prefix, mTo / mFrom - 1))
  }else if (diff == 0) {
    return (priceTo / priceFrom - 1)
  }
}