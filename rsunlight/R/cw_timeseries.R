#' Find the popularity of a phrase over a period of time.
#'
#' @export
#' @template cw
#' @template cw_timeseries
#' @return A data.frame with number of times (and percentages of the text result,
#' if selected).
#' @examples \dontrun{
#' # Search for a phrase to get a time series of
#' cw_timeseries(phrase='climate change')
#'
#' # Get percentages back, which are not returned by default
#' cw_timeseries(phrase='climate change', percentages='true')
#'
#' # Get a list of how many times the phrase "united states" was said by
#' # legislators from Virginia on each day of the most recent Congress:
#' cw_timeseries(phrase='united states', entity_type='state',
#' entity_value='VA')
#'
#' # Get a list of how many times the phrase "united States" appears in the
#' # Congressional Record on each day between Jan. 1, 2010, and June 1, 2010:
#' cw_timeseries(phrase='united states', start_date='2009-01-01',
#' end_date='2009-06-01')
#'
#' # Get a list of how many times the phrase "united states" appears in the
#' # Congressional Record in each month between January and June, 2010:
#' cw_timeseries(phrase='united states', start_date='2009-01-01',
#' end_date='2009-04-30', granularity='month')
#'
#' # Plot data
#' library('ggplot2')
#' dat <- cw_timeseries(phrase='climate change')
#' ggplot(dat, aes(day, count)) +
#'    geom_line() +
#'    theme_grey(base_size=20)
#'
#' dat_d <- cw_timeseries(phrase = 'climate change', party = "D")
#' dat_d$party <- rep("D", nrow(dat_d))
#' dat_r <- cw_timeseries(phrase = 'climate change', party = "R")
#' dat_r$party <- rep("R", nrow(dat_r))
#' dat_both <- rbind(dat_d, dat_r)
#' ggplot(dat_both, aes(day, count, colour = party)) +
#'    geom_line() +
#'    theme_grey(base_size = 20) +
#'    scale_colour_manual(values=c("blue", "red"))
#' }

cw_timeseries <- function(phrase=NULL, date = NULL, start_date=NULL, end_date=NULL,
  chamber=NULL, state=NULL, party=NULL, bioguide_id=NULL, mincount=NULL,
  percentages=NULL, granularity='day', entity_type=NULL, entity_value=NULL, as = 'table',
  key = NULL, ...) {

  key <- check_key(key)
  splitt <- function(x) paste(str_sub(x, 1, 4), "-", str_sub(x, 5, 6), sep = "")
  args <- sc(list(apikey=key, phrase=phrase, start_date=start_date,
                       date = date, end_date=end_date, chamber=chamber, state=state,
                       party=party, bioguide_id=bioguide_id, mincount=mincount,
                       percentages=percentages, granularity=granularity,
                       entity_type=entity_type, entity_value=entity_value))

  tmp <- give_cg(as, cwurl(), "/dates.json", args, ...)
  if (as %in% c('response', 'list')) {
    tmp
  } else {
    if (granularity == 'day') {
      tmp$day <- as.Date(tmp$day)
    } else if (granularity == 'month') {
      tmp$month <- as.Date(sprintf("%s-01", sapply(tmp$month, splitt)))
    } else if (granularity == 'year') {
      tmp$year <- as.Date(sprintf("%s-01-01", tmp$year))
    }
    tmp
  }
}
