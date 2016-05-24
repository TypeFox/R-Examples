#' Retrieve filtered votes from a database
#'
#' Function \code{get_filtered_votes} reads filtered votes from a database.
#'
#' @details
#' Function \code{get_filtered_votes} reads filtered votes from a database.
#' The result of this function is an invisible data frame with statements' data.
#'
#' Possible filters:
#' \enumerate{
#' \item clubs - names of clubs. This filter is a character vector with elements
#' like for example: 'PO', 'PiS', 'SLD'. It is possible to choose more than one club.
#' \item dates - period of time. This filter is a character vector with two elements
#' in date format 'YYYY-MM-DD', where the first describes left boundary of period and
#' the second right boundary. It is possible to choose only one day, just try the same
#' date as first and second element of vector.
#' \item terms_of_office - range of terms of office's numbers. This filter is a integer
#' vector with two elements, where the first describes a left boundary of range
#' and the second a right boundary. It is possible to choose only one term of office,
#' just try the same number as first and second element of vector.
#' \item meetings - range of meetings' numbers. This filter is a integer vector with two
#' elements, where the first describes a left boundary of range and the second a right
#' boundary. It is possible to choose only one meeting, just try the same number
#' as first and second element of vector.
#' \item votings - range of votings' numbers. This filter is a integer vector with two
#' elements, where the first describes a left boundary of range and the second a right
#' boundary. It is possible to choose only one voting, just try the same number
#' as first and second element of vector.
#' \item deputies - full names of deputies. This filter is a character vector with full
#' names of deputies in format: 'surname first_name second_name'. If you are not sure
#' if the deputy you were thinking about has second name, try 'surname first_name' or
#' just 'surname'. There is high probability that proper deputy will be chosen.
#' It is possible to choose more than one deputy.
#' \item topics - text patterns. This filter is a character vector with text patterns of
#' topics that you are interested about. Note that the votings' topics are written like
#' sentences, so remember about case inflection of nouns and adjectives and use stems of
#' words as patterns. For example if you want to find votings about education (in Polish:
#' szkolnictwo) try 'szkolnictw'. It is possible to choose more than one pattern.}
#'
#' If you did not choose any filter, the whole database will be downloaded.
#' Note that, due to data size (<= ~150 MB) it may take few seconds / minutes
#' to download all votes.
#'
#' Because of encoding issue on Windows operation system, you also need to select
#' if you use Windows.
#'
#' @usage get_filtered_votes(dbname = 'sejmrp', user = 'reader',
#'   password = 'qux94874', host = 'services.mini.pw.edu.pl',
#'   windows = .Platform$OS.type == 'windows', clubs = character(0),
#'   dates = character(0), terms_of_office = integer(0), 
#'   meetings = integer(0), votings = integer(0),
#'   deputies = character(0), topics = character(0))
#'
#' @param dbname name of database; default: 'sejmrp'
#' @param user name of user; default: 'reader'
#' @param password password of database; default: 'qux94874'
#' @param host name of host; default: 'services.mini.pw.edu.pl'
#' @param windows information of used operation system; default: .Platform$OS.type == 'windows'
#' @param clubs names of clubs that will be taken to filter data from database;
#' default: character(0)
#' @param dates period of time that will be taken to filter data from database;
#' default: character(0)
#' @param terms_of_office range of terms of office's numbers that will be taken to filter data
#' from database; default: integer(0)
#' @param meetings range of meetings' numbers that will be taken to filter data from database;
#' default: integer(0)
#' @param votings range of votings' numbers that will be taken to filter data from database;
#' default: integer(0)
#' @param deputies full names of deputies that will be taken to filter data from database;
#' default: character(0)
#' @param topics text patterns that will be taken to filter data from database;
#' default: character(0)
#'
#' @return data frame with NULL
#'
#' @examples
#' \dontrun{
#' filtered_votes <- get_filtered_votes()
#' dim(filtered_votes)
#' # [1] 2826483       9
#' names(filtered_votes)
#' [1] 'surname_name' 'nr_term_of_office' 'club' 'vote' 'id_voting'
#' [6] 'nr_meeting' 'nr_voting' 'date_meeting' 'topic_voting'
#' object.size(filtered_votes)
#' # 148694336 bytes}
#'
#' @note
#' Default parameters use privilages of 'reader'. It can only SELECT data from database.
#'
#' All information is stored in PostgreSQL database.
#'
#' @author Piotr Smuda
#'
#' @export
#'
#' @importFrom dplyr src_postgres
#' @importFrom dplyr tbl
#' @importFrom dplyr sql
#' @importFrom dplyr filter
#' @importFrom dplyr between
#' @importFrom dplyr mutate
#' @importFrom dplyr collect
#'

get_filtered_votes <- function(dbname = "sejmrp", user = "reader", password = "qux94874", host = "services.mini.pw.edu.pl",
                               windows = .Platform$OS.type == "windows", clubs = character(0), dates = character(0),
                               terms_of_office = integer(0), meetings = integer(0), votings = integer(0), deputies = character(0),
                               topics = character(0)) {
    stopifnot(is.character(dbname), is.character(user), is.character(password), is.character(host), is.logical(windows), is.character(clubs),
        is.character(dates), is.numeric(terms_of_office), is.numeric(meetings), is.numeric(votings), is.character(deputies), is.character(topics),
        all(c(terms_of_office, meetings, votings)%%1 == 0))
    length_clubs <- length(clubs)
    length_dates <- length(dates)
    length_terms_of_office <- length(terms_of_office)
    length_meetings <- length(meetings)
    length_votings <- length(votings)
    length_deputies <- length(deputies)
    length_topics <- length(topics)
    stopifnot(length_clubs >= 0, length_dates == 0 | length_dates == 2, length_terms_of_office == 0 | length_terms_of_office == 2,
              length_meetings == 0 | length_meetings == 2, length_votings == 0 | length_votings == 2, length_deputies >= 0, length_topics >= 0)

    # connecting to database to add information about new SELECT to the counter table
    drv <- dbDriver("PostgreSQL")
    database_diet <- dbConnect(drv, dbname = dbname, user = user, password = password, host = host)

    # add information about new SELECT to the counter table
    dbSendQuery(database_diet, paste0("INSERT INTO counter (what, date) VALUES ('filt_votes','", Sys.Date(), "')"))

    suppressWarnings(dbDisconnect(database_diet))

    # fake variables in order to pass CRAN CHECK
    club <- NULL
    date_meeting <- NULL
    nr_term_of_office <- NULL
    nr_meeting <- NULL
    nr_voting <- NULL
    surname_name <- NULL
    topic_voting <- NULL
    `%SIMILAR TO%` <- NULL

    # connecting to database with dplyr to get votes
    database_diet <- src_postgres(dbname = dbname, user = user, password = password, host = host)

    # read data dodac potem
    votes <- tbl(database_diet, sql("SELECT d.surname_name, d.nr_term_of_office, v.club, v.vote, vv.id_voting,
                                    vv.nr_meeting, vv.nr_voting, vv.date_meeting, vv.topic_voting
                                    FROM votes v, votings vv, deputies d
                                    WHERE v.id_voting = vv.id_voting AND v.nr_term_of_office = vv.nr_term_of_office
                                    AND d.id_deputy = v.id_deputy AND d.nr_term_of_office = vv.nr_term_of_office"))

    # clubs filter
    if (length_clubs > 0) {
        clubs <- paste0("(%", clubs, "%)")
        clubs <- paste0(clubs, collapse = "|")
        votes <- filter(votes, club %SIMILAR TO% clubs)
    }

    # dates filter
    if (length_dates == 2) {
        votes <- filter(votes, between(date_meeting, dates[1], dates[2]))
    }

    # terms_of_office filter
    if (length_terms_of_office == 2) {
      votes <- filter(votes, between(nr_term_of_office, terms_of_office[1], terms_of_office[2]))
    }
    
    # meetings filter
    if (length_meetings == 2) {
        votes <- filter(votes, between(nr_meeting, meetings[1], meetings[2]))
    }

    # votings filter
    if (length_votings == 2) {
        votes <- filter(votes, between(nr_voting, votings[1], votings[2]))
    }

    # deputies filter
    if (length_deputies > 0) {
        # changing polish characters for any character
        deputies <- stri_replace_all_regex(deputies, "[^a-zA-Z %]", "_")
        deputies <- paste0("(%", deputies, "%)")
        deputies <- paste0(deputies, collapse = "|")
        votes <- filter(votes, surname_name %SIMILAR TO% deputies)
    }

    # topics filter
    if (length_topics > 0) {
        # changing polish characters for any character
        topics <- stri_replace_all_regex(topics, "[^a-zA-Z %]", "_")
        topics <- paste0("(%", topics, "%)")
        topics <- paste0(topics, collapse = "|")
        votes <- filter(votes, topic_voting %SIMILAR TO% topics)
    }

    # reading data
    votes <- as.data.frame(collect(votes, stringsAsFactors = FALSE))

    # if empty result of query
    if (nrow(votes) == 0) {
        suppressWarnings(dbDisconnect(database_diet$con))
        return(votes)
    }

    # encoding for windows
    if (windows) {
        votes[, 1] <- iconv(votes[, 1], from = "UTF-8", to = "Windows-1250")
        votes[, 3] <- iconv(votes[, 3], from = "UTF-8", to = "Windows-1250")
        votes[, 4] <- iconv(votes[, 4], from = "UTF-8", to = "Windows-1250")
        votes[, 9] <- iconv(votes[, 9], from = "UTF-8", to = "Windows-1250")
    }

    suppressWarnings(dbDisconnect(database_diet$con))
    return(invisible(votes))
}
