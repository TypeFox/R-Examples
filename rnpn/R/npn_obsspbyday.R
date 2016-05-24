#' Get observations by day for a particular species or set of species.
#'
#' @export
#' @param speciesid Required. Species id numbers, from 1 to infinity, potentially,
#'     use e.g., c(52, 53, etc.) if more than one species desired (numeric)
#' @param startdate start date of data period desired, see format in examples (character)
#' @param enddate end date of data period desired, see format in examples (character)
#' @template curl
#' @return Number of observations by day, in a list
#' @examples \dontrun{
#' out <- npn_obsspbyday(speciesid=357, startdate='2010-04-01', enddate='2012-01-05')
#' head(out[[1]])
#'
#' # Lookup names
#' temp <- lookup_names(name='bird', type='common')
#' comnames <- temp[temp$species_id %in% c(357, 359, 1108), 'common_name']
#'
#' out <- npn_obsspbyday(speciesid=c(357, 359, 1108), startdate='2010-04-01', enddate='2013-09-31')
#' names(out) <- comnames
#' library("plyr")
#' df <- ldply(out)
#' df$date <- as.Date(df$date)
#'
#' library('ggplot2')
#' ggplot(df, aes(date, count)) +
#'  geom_line() +
#'  theme_grey(base_size=20) +
#'  facet_grid(.id ~.)
#' }

npn_obsspbyday <- function(speciesid = NULL, startdate = NULL, enddate = NULL, ...) {
  args <- npnc(list(start_date = startdate, end_date = enddate))
  for (i in seq_along(speciesid)) {
    args[paste('species_id[',i,']', sep = '')] <- speciesid[i]
  }
  tt <- npn_GET(paste0(base(), 'observations/getObservationsForSpeciesByDay.json'), args, ...)
  df_list <- lapply(tt$all_species$species, function(x) {
    setDF(rbindlist(x[2]$count_list, fill = TRUE, use.names = TRUE))
  })
  df_list <- lapply(df_list, function(x){
    x$date <- gsub("\\s.+", "", x$date)
    x
  })
  structure(df_list, .Names = speciesid)
}
