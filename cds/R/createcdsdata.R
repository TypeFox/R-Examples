#' Create a cdsdata Object
#' 
#' Create a cdsdata object from a data frame or matrix.
#' 
#' @param x A data frame or matrix containing the data.
#' @param q Optional; the maximum rating category, so that the rating scale used for all 
#' items are \code{1:q}.
#' @keywords multivariate
#' @export createcdsdata
createcdsdata <- function(x, q = NULL){
  # Capture call
  cll <- match.call()
  
  # Type-check and coerce to matrix
  stopifnot(is.data.frame(x) | is.matrix(x))
  if(is.data.frame(x)) x <- as.matrix.data.frame(x)
  
  # Check for missings
  if(any(is.na(x))) stop("Missings not allowed")
  
  # Obtain dimensions and q
  n <- nrow(x)
  m <- ncol(x)
  if(is.null(q)){
    q <- max(x, na.rm = TRUE)
    message("Setting q = max(x, na.rm = TRUE) since none supplied")
  }
  scales <- 1:q
  
  # Calculate Fr.cent.rs and Fr.rs
  bounds <- matrix(scales[-q] + 0.5, byrow = TRUE, ncol = q  - 
                     1, nrow = n)
  colnames(bounds) <- paste0("b", 1:(q - 1))
  Tmat.rs <- t(apply(cbind(x, bounds), 1, rank, na.last = "keep")) - 1
  Fr.rs <- rbind(Tmat.rs, m + q - 2 - Tmat.rs)
  Fr.cent.rs <- Fr.rs - 0.5*(m + q - 2)*matrix(1, nrow = 2*n, ncol = m + q - 1)
  
  # Create output object and set class
  out <- list(postrs = x, postbl = x, Fr.cent.rs = Fr.cent.rs, Fr.cent.bl = Fr.cent.rs, 
              Fr.rs = Fr.rs, Fr.bl = Fr.rs, m0 = m, munique = 0,
              scales = scales, m = m, call = cll)
  attr(out, "class") <- c("cdsdata", "icdsdata", "list")
  return(out)
}
