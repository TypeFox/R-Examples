#' Select and return valid dummy outlier masks
#'
#' Valid dummy outlier masks are integers whose bitwise AND with the given \code{invalid_mask} is zero. This function returns the subset of integers of the given vector that does not share any bits with the given \code{invalid_mask}.
#' @param all_outlier_masks A vector of possible outlier masks (integers).
#' @param invalid_mask An integer encoding the invalid columns.
#' @return The valid outlier_masks as a vector of integers.
#' @examples
#' all_outlier_masks <- c(0, 1, 2, 3, 4, 5, 6, 7)
#' invalid_mask <- 1
#' autovarCore:::select_valid_masks(all_outlier_masks, invalid_mask)
#' @export
select_valid_masks <- function(all_outlier_masks, invalid_mask) {
  if (invalid_mask == 0)
    return(all_outlier_masks)
  result <- NULL
  for (outlier_mask in all_outlier_masks)
    if (bitwAnd(outlier_mask, invalid_mask) == 0)
      result <- c(result, outlier_mask)
  result
}
