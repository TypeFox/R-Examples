#' Get dimensions of all data sets
#'
#' Get the size of the data sets behind each of the GUID's offered by the base
#' URL. This function will iterate through all the GUID's available and present
#' the results as a data frame with the GUID, the number of rows, number of
#' columns and dimension (as total number of cells) for each data set.
#'
#' @param base_url The base URL of the Junar service
#' @param api_key The user's API key for the Junar service
#' @keywords GUID
#' @export

get_dimensions <- function(base_url, api_key) {
  if (missing(base_url)) {
    warning("Please add a valid base URL")
  }
  if (missing(api_key)) {
    warning("Please add a valid API key for the base URL you are trying to access")
  }
  try({
    guid_list <- list_guid(base_url, api_key)

    for (guid in guid_list) {
     current_set <- get_data(base_url, api_key, guid)
     result_row <- c("GUID" = guid, "NROW" = nrow(current_set), "NCOL" = ncol(current_set),
                     "DIM" = nrow(current_set)*ncol(current_set))
     if (!exists("result_df")) {
       # TODO: There must be a better way to do this
       result_df <- data.frame(GUID = 0, NROW = 0, NCOL = 0, DIM = 0)
       result_df <- rbind(result_df, result_row)
       result_df <- result_df[-1,]
      } else {
       result_df <- rbind(result_df, result_row)
      }
    }
    result_df$NROW <- as.numeric(result_df$NROW)
    result_df$NCOL <- as.numeric(result_df$NCOL)
    result_df$DIM  <- as.numeric(result_df$DIM)
    return(result_df)
  })
}
