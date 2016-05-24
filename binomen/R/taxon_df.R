#' Taxon data.frame
#'
#' @export
#' @param x A data.frame of taxa
#' @examples
#' # subset data.frame using taxonomy
#' df <- data.frame(family=c('Asteraceae','Asteraceae','Asteraceae','Poaceae','Poaceae','Poaceae'),
#'                  tribe=c('Helianthi','Helianthi','Helianthi','Poaeae','Festuci','Poaeae'),
#'                  genus=c('Helianthus','Helianthus','Madia','Poa','Festuca','Holodiscus'),
#'                  stringsAsFactors = FALSE)
#' df2 <- taxon_df(df)
#' df2 %>% pick(family)
#' df2 %>% pick(genus, tribe)

taxon_df <- function(x){
  structure(x, class = c('taxondf', 'data.frame'))
}
