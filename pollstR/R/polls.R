#'@include pollstr-package.R
NULL

# Create URL for the charts API method
pollstr_polls_url <- function(page, chart, state, topic, before, after,
                              sort, showall) {
  query <- list()
  if (! is.null(page)) {
    query[["page"]] <- as.character(page)[1]
  }
  if (! is.null(chart)) {
    query[["chart"]] <- as.character(chart)[1]
  }
  if (! is.null(state)) {
    query[["state"]] <- as.character(state)[1]
  }
  if (! is.null(topic)) {
    query[["topic"]] <- as.character(topic)[1]
  }
  if (! is.null(before)) {
    before <- before[1]
    if (inherits(before, "Date")) before <- format(before, "%Y-%m-%d")
    query[["before"]] <- as.character(before)[1]
  }
  if (! is.null(after)) {
    after <- after[1]
    if (inherits(after, "Date")) after <- format(after, "%Y-%m-%d")
    query[["after"]] <- as.character(after)[1]
  }
  if (sort) {
    query[["sort"]] <- "updated"
  }
  if (! is.null(showall)) {
    query[["showall"]] <- if (showall) "true" else "false"
  }
  if (! length(query)) {
    query <- NULL
  }
  modify_url(paste(.POLLSTR_API_URL, "polls", sep="/"), query = query)
}

polls2df <- function(.data) {

  # Remove polls without questions or an error occurs in the conversion
  .data <- .data[lapply(.data, function(x) length(x[["questions"]])) > 0]

  polls <- ldply(.data,
                 function(x) {
                   y <- convert_df(x[setdiff(names(x),
                                             c("questions", "survey_houses",
                                               "sponsors"))])
                   y[["start_date"]] <- as.Date(y[["start_date"]])
                   y[["end_date"]] <- as.Date(y[["end_date"]])
                   y[["last_updated"]] <- as.POSIXct(y[["last_updated"]],
                                                     "%Y-%m-%dT%H:%M:%OSZ",
                                                     tz = "GMT")
                   y
                 })
  
  # Convert polls
  for (i in c("id")) {
    polls[[i]] <- as.integer(polls[[i]])
  }
  
  clean_subpopulations <- function(x) {
    merge(convert_df(x[c("name", "observations", "margin_of_error")]),
          ldply(x[["responses"]], convert_df))
  }
  
  clean_questions <- function(x) {
    subpops <- ldply(x[["subpopulations"]], clean_subpopulations)
    subpops <- rename(subpops, c(name = "subpopulation"))
    merge(convert_df(x[c("name", "chart", "topic", "state")]),
          subpops)
  }
  
  questions <-
    ldply(.data,
          function(x) {
            ques <- rename(ldply(x[["questions"]], clean_questions),
                           c(name = "question"))
            ques[["id"]] <- x[["id"]]
            ques
          })
  # convert
  for (i in c("observations", "id")) {
    questions[[i]] <- as.integer(questions[[i]])
  }
  
  clean_sponsors <- function(x) {
    sponsors <- x[["sponsors"]]
    if (length(sponsors)) {
      sponsors <- ldply(sponsors, convert_df)
      sponsors[["id"]] <- x[["id"]]
      sponsors
    } else {
      NULL
    }
  }   
  sponsors <- ldply(.data, clean_sponsors)
  
  clean_survey_houses <- function(x) {
    survey_houses <- x[["survey_houses"]]
    if (length(survey_houses)) {
      survey_houses <- ldply(survey_houses, convert_df)
      survey_houses[["id"]] <- x[["id"]]
      survey_houses
    } else {
      NULL
    }
  }
  survey_houses <- ldply(.data, clean_survey_houses)
  
  
  structure(list(polls = polls,
                 questions = questions,
                 survey_houses = survey_houses,
                 sponsors = sponsors),
            class = "pollstr_polls")
}

get_poll <- function(page, chart, state, topic, before, after, sort, showall,
                     as = "parsed") {
  url <- pollstr_polls_url(page, chart, state, topic, before, after, sort,
                           showall)
  get_url(url, as = as)
}

#' Get a list of polls
#'
#' @param page Return page number
#' @param chart List polls related to the specified chart. Chart names are the \code{slug} returned by \code{pollstr_charts}.
#' @param state Only include charts from a single state. Use 2-letter pstate abbreviations. "US" will return all national charts.
#' @param topic Only include charts related to a specific topic. See the \url{http://elections.huffingtonpost.com/pollster/api} for examples.
#' @param before Only list polls that ended on or bfore the specified date.
#' @param after Only list polls that ended on or bfore the specified date.
#' @param sort If \code{TRUE}, then sort polls by the last updated time.
#' @param showall Include polls for races that were once possible but didn't happen (e.g. Gingrich vs. Obama 2012)
#' @param max_pages Maximum number of pages to get.
#' @param convert Rearrange the data returned by the API into easier to use data frames.
#'
#' @references \url{http://elections.huffingtonpost.com/pollster/api}
#' @return If \code{convert=TRUE}, a \code{"pollstr_polls"} object with elements
#' \describe{
#' \item{\code{polls}}{A \code{data.frame} with entries for each poll.}
#' \item{\code{questions}}{A \code{data.frame} with entries for each question asked in the polls.}
#' \item{\code{survey_houses}}{A \code{data.frame} with the survey houses of the polls. There can be multiple survey houses for a poll.}
#' \item{\code{sponsors}}{A \code{data.frame} with the sponsors of the polls. Not all polls have sponsors.}
#' }
#' Otherwise, a \code{"list"} in the original structure of the json returned by the API.
#' @examples
#' \dontrun{
#' # Get polls related to a chart pulled programmatically with
#' # pollstr_charts()
#' all_charts <- pollstr_charts()
#' pollstr_polls(chart=all_charts$slug[1])
#' # Lookup polls related to a specific topic
#' pollstr_polls(topic='2016-president')
#' }
#' @export
pollstr_polls <- function(page = 1, chart = NULL, state = NULL,
                          topic = NULL, before = NULL, after = NULL,
                          sort = FALSE, showall = NULL, max_pages = 1,
                          convert = TRUE) {
  .data <- list()
  i <- 0L
  while (i < max_pages) {
    newdata <- get_poll(page + i, chart, state, topic, before, after, sort,
                        showall)
    if (length(newdata)) {
      .data <- append(.data, newdata)
    } else {
      break
    }
    i <- i + 1L
  }
  if (convert) .data <- polls2df(.data)
  .data
}
