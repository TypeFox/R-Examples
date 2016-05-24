#'
#' Partial parser for game-log files
#'
#' Instead of returning the entire file, this function allows the user
#' to choose the columns and date for game-log data.
#'
#' @param year A single four-digit year.
#' @param glFields character.  The desired game-log columns.  This should be a
#' subset of \code{gamelogFields}, and \strong{not} the entire vector.
#' @param date One of either NULL (the default), or a single four-digit
#' character string identifying the date 'mmdd'
#' @param ... further arguments passed to \code{\link[utils]{download.file}}
#'
#' @importFrom data.table fread
#' @importFrom data.table setnames
#'
#' @export
#'
#' @return
#' \itemize{
#' \item{\code{getPartialGamelog}}{ - A data table with dimensions \code{length(date)} x \code{length(glFields)} if
#' \code{date} is not NULL, otherwise the row dimension is the nuber of games for the given year.}
#' \item{\code{gamelogFields}}{ - A character vector of possible values to choose from for the
#' \code{glFlields} argument in \code{getPartialGamelog}.}
#' }
#'
#' @examples ## Get Homerun and RBI info for the 2012 season, with park ID
#' f <- grep("HR|RBI|Park", gamelogFields, value = TRUE)
#' getPartialGamelog(2012, glFields = f)
#'
#' ## Get Homerun and RBI info for August 25, 2012 - with park ID
#' getPartialGamelog(glFields=f, date = "20120825")
#'
getPartialGamelog <- function(year, glFields, date = NULL, ...) {

    ## check 'glFields' against package variable 'retrosheetFields$gamelog'
    if(identical(glFields, retrosheetFields$gamelog)) {
        stop(shQuote("getPartialGamelog"), " is for efficiently returning a small subset of the entire file. For the full table, use ", shQuote("getRetrosheet(\"game\", year)"))
    }

    if(missing(year)) {
        year <- substr(date, 1L, 4L)
    }

    ## define the url
    u <- "http://www.retrosheet.org"
    full <- sprintf("%s/gamelogs/gl%s.zip", u, year)

    ## download the file
    tmp <- tempfile()
    on.exit(unlink(tmp))
    download.file(full, destfile = tmp, ...)

    ## extract the text file
    fname <- unzip(tmp, list = TRUE)$Name
    unzip(tmp, files = fname)
    on.exit(unlink(fname), add = TRUE)

    ## match 'glFields' against the internal name vector
    sel <- union(1L, sort(match(glFields, retrosheetFields$gamelog)))

    ## read the data
    out <- if(is.null(date)) {

        fread(fname, select = sel, header = FALSE)

    } else if(is.character(date)) {

        ## get the first column only - this is the 'Date' column
        sc <- scan(fname, sep = ",", flush = TRUE, what = "", quote = "\"", quiet = TRUE)

        ## find rows of matched dates
        if(!length(wh <- which(sc == date))) {
            stop("invalid 'date' given")
        }

        ## read the file - selecting specified date and columns
        suppressWarnings(fread(fname, select = sel, header = FALSE,
            skip = min(wh)-1L, nrows = diff(range(wh))))

    }

    ## set the names
    setnames(out, retrosheetFields$gamelog[sel])

    ## return the table
    out
}

#' @rdname getPartialGamelog
#'
#' @name gamelogFields
#'
#' @export
#'
gamelogFields <- c("Date", "DblHdr", "Day", "VisTm", "VisTmLg", "VisTmGNum", "HmTm",
    "HmTmLg", "HmTmGNum", "VisRuns", "HmRuns", "NumOuts", "DayNight",
    "Completion", "Forfeit", "Protest", "ParkID", "Attendance", "Duration",
    "VisLine", "HmLine", "VisAB", "VisH", "VisD", "VisT", "VisHR",
    "VisRBI", "VisSH", "VisSF", "VisHBP", "VisBB", "VisIBB", "VisK",
    "VisSB", "VisCS", "VisGDP", "VisCI", "VisLOB", "VisPs", "VisER",
    "VisTER", "VisWP", "VisBalks", "VisPO", "VisA", "VisE", "VisPassed",
    "VisDB", "VisTP", "HmAB", "HmH", "HmD", "HmT", "HmHR", "HmRBI",
    "HmSH", "HmSF", "HmHBP", "HmBB", "HmIBB", "HmK", "HmSB", "HmCS",
    "HmGDP", "HmCI", "HmLOB", "HmPs", "HmER", "HmTER", "HmWP", "HmBalks",
    "HmPO", "HmA", "HmE", "HmPass", "HmDB", "HmTP", "UmpHID", "UmpHNm",
    "Ump1BID", "Ump1BNm", "Ump2BID", "Ump2BNm", "Ump3BID", "Ump3BNm",
    "UmpLFID", "UmpLFNm", "UmpRFID", "UmpRFNm", "VisMgrID", "VisMgrNm",
    "HmMgrID", "HmMgrNm", "WinPID", "WinPNm", "PID", "PNAme", "SavePID",
    "SavePNm", "GWinRBIID", "GWinRBINm", "VisStPchID", "VisStPchNm",
    "HmStPchID", "HmStPchNm", "VisBat1ID", "VisBat1Nm", "VisBat1Pos",
    "VisBat2ID", "VisBat2Nm", "VisBat2Pos", "VisBat3ID", "VisBat3Nm",
    "VisBat3Pos", "VisBat4ID", "VisBat4Nm", "VisBat4Pos", "VisBat5ID",
    "VisBat5Nm", "VisBat5Pos", "VisBat6ID", "VisBat6Nm", "VisBat6Pos",
    "VisBat7ID", "VisBat7Nm", "VisBat7Pos", "VisBat8ID", "VisBat8Nm",
    "VisBat8Pos", "VisBat9ID", "VisBat9Nm", "VisBat9Pos", "HmBat1ID",
    "HmBat1Nm", "HmBat1Pos", "HmBat2ID", "HmBat2Nm", "HmBat2Pos",
    "HmBat3ID", "HmBat3Nm", "HmBat3Pos", "HmBat4ID", "HmBat4Nm",
    "HmBat4Pos", "HmBat5ID", "HmBat5Nm", "HmBat5Pos", "HmBat6ID",
    "HmBat6Nm", "HmBat6Pos", "HmBat7ID", "HmBat7Nm", "HmBat7Pos",
    "HmBat8ID", "HmBat8Nm", "HmBat8Pos", "HmBat9ID", "HmBat9Nm",
    "HmBat9Pos", "Additional", "Acquisition")

