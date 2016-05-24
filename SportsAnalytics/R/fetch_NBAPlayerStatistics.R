

#' Fetch player statistics for NBA season statistics
#'
#' @param season Season to fetch (\code{"88-89"}, \code{"89-90"},
#    \code{...})
#' @param what Complete, home or away statistics
#'
#' @return Data frame; see \code{\link{NBAPlayerStatistics0910}} for a
#'   description.
#'
#' @examples
#'   \dontrun{
#'     fetch_NBAPlayerStatistics("07-08")
#'     fetch_NBAPlayerStatistics("07-08", what = ".Home")
#'     fetch_NBAPlayerStatistics("07-08", what = ".Away")
#'   }
#'
#' @export
fetch_NBAPlayerStatistics <- function(season = "09-10", what = c("", ".Home", ".Away")) {
  capwords <- function(s, strict = FALSE) {
    cap <- function(s)
      paste(toupper(substring(s,1,1)),
        {s <- substring(s,2); if(strict) tolower(s) else s}, sep = "", collapse = " " )

    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
  }

  what <- match.arg(what)

  url <- sprintf("http://dougstats.com/%s%sRD.txt", season, what)

  stats <- read.table(url, header = TRUE, stringsAsFactors = FALSE)

  stats$Team <- toupper(stats$Team)
  stats$Team <- factor(stats$Team)

  stats$Player <- local({
    name <- strsplit(stats$Player, ",")
    name <- lapply(name, capwords)
    name <- sapply(name, function(x) paste(rev(x), collapse = " "))

    name
  })

  stats$PS[stats$PS == "??"] <- NA
  stats$PS <- factor(stats$PS)

  stats$League <- "NBA"
  stats$League <- factor(stats$League)

  cn <- .colname.NBAPlayerStatistics()
  stats <- stats[, names(cn)]
  names(stats) <- unname(cn)

  stats
}



.colname.NBAPlayerStatistics <- function() {
  ## http://dougstats.com/ExelNotes.html
  c(League = "League",
    Player = "Name",
    Team   = "Team",
    PS     = "Position",
    GP     = "GamesPlayed",
    Min    = "TotalMinutesPlayed",
    FGM    = "FieldGoalsMade",
    FGA    = "FieldGoalsAttempted",
    X3M    = "ThreesMade",
    X3A    = "ThreesAttempted",
    FTM    = "FreeThrowsMade",
    FTA    = "FreeThrowsAttempted",
    OR     = "OffensiveRebounds",
    TR     = "TotalRebounds",
    AS     = "Assists",
    ST     = "Steals",
    TO     = "Turnovers",
    BK     = "Blocks",
    PF     = "PersonalFouls",
    DQ     = "Disqualifications",
    PTS    = "TotalPoints",
    TC     = "Technicals",
    EJ     = "Ejections",
    FF     = "FlagrantFouls",
    Sta    = "GamesStarted")
}
