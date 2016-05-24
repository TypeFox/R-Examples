
#' Format a df of addresses for upload
#'
#' Format a df of addresses for upload
#' @param addresses_df A df with two columns, a unique id and address character strings to be geocoded. Only results with valid lat, lng are returned, so the unique id is used to match back to the data later on.
#' @return a long string consisting of all of the addresses to be geocoded formatted for the POST request in geocodeHERE_batch_upload()
#' @keywords geocode batch
#' @export
#' @examples
#' addresses <- chicago_landmarks[,"Address"]
#' addresses <- paste(addresses, "chicago IL")
#' addresses_df <- data.frame(id=1:length(addresses), addresses=addresses)
#' address_str <- df_to_string(addresses_df)
df_to_string <- function(addresses_df){
  if(nrow(addresses_df) > 9999){stop("'addresses_df' must be less than 10,000 rows")}

  header <- "recID|searchText"
  therest <- paste(paste(addresses_df[,1], addresses_df[,2], sep="|"), collapse="\n")
  final <- paste(header, therest, sep="\n")

  return(final)
}



