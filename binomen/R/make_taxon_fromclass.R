#' Make taxon from class
#'
#' @export
#' @param data A data.frame
#' @param authority Taxonomic authority
#' @examples \dontrun{
#' df <- data.frame(rank=c('family','tribe','subtribe','genus','subgenus','species'),
#'                  name=c('Helianthi','Helianthi','Helianthi','Poa','Festuci','Poa annua'),
#'                  id=c(1,2,3,4,5,6),
#'                  stringsAsFactors = FALSE)
#' apply(df, 1, make_taxon_fromclass)
#' }
make_taxon_fromclass <- function(data, authority="none"){
  rank <- data[['rank']]
  name <- data[['name']]
  id <- as.numeric(data[['id']])
  res <- list(taxonref(rank = rank, name = name, id = id))
  names(res) <- rank
  res[1]
}
