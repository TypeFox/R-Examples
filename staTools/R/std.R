#' Unity-based Normalization
#'
#' Unity-based normalization of a vector.
#' @param x A vector to normalize.
#' @keywords plot fit powerlaw distributions
#' @export std
#' @examples
#' x = moby
#' z = std(x)

std = function(x)
{
  z = (x-min(x)) / (max(x)-min(x))
  return(z)
}