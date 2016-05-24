#' @title Score the FACT-GP
#' 
#' @description
#' Generates all of the scores of the Functional Assessment of Cancer 
#' Therapy-General (FACT-GP) from item responses.
#' 
#' @details
#' Given a data frame that includes all of the FACT-GP items as variables,
#' appropriately named, this function generates all of the FACT-GP scale
#' scores.  It is crucial that the item variables in the supplied data frame
#' are named according to FACT conventions.  For example, the first physical
#' well-being item should be named GP1, the second GP2, and so on.  Please
#' refer to the materials provided by \url{http://www.facit.org} for the
#' particular questionnaire you are using.  In particular, refer to the left
#' margin of the official questionnaire (i.e., from facit.org) for the
#' appropriate item variable names.
#' 
#' @section Note:
#' Keep in mind that this function (and R in general) is case-sensitive.
#'
#' All variables should be in numeric or integer format.
#'
#' This scoring function expects missing item responses to be coded as NA,
#' 8, or 9, and valid item responses to be coded as 0, 1, 2, 3, or 4.  Any
#' other value for any of the items will result in an error message and no
#' scores.
#'
#' Some item variables are reverse coded for the purpose of generating the
#' scale scores.  The official (e.g., from \url{http://www.facit.org}) SAS
#' and SPSS scoring algorithms for this questionnaire automatically replace
#' the original items with their reverse-coded versions.  This can be
#' confusing if you accidentally run the algorithm more than once on your
#' data.  As its default, \code{scoreFACT_GP} DOES NOT replace any of your
#' original item variables with the reverse coded versions.  However, for
#' consistentcy with the behavior of the other versions on
#' \url{http://www.facit.org}, the \code{updateItems} argument is
#' provided.  If set to \code{TRUE}, any item that is supposed to be
#' reverse coded will be replaced with its reversed version in the data
#' frame returned by \code{scoreFACT_GP}.
#' 
#' 
#' @param df A data frame with the FACT-GP items, appropriately-named.
#' @param updateItems Logical, if \code{TRUE} any original item that is
#' reverse coded for scoring will be replaced by its reverse coded version
#' in the returned data frame, and any values of 8 or 9 will be replaced
#' with NA.  The default, \code{FALSE}, returns the original items
#' unmodified.  
#' @param keepNvalid Logical, if \code{TRUE} the function
#' returns an additional variable for most of the returned scale scores
#' containing the number of valid, non-missing responses from each
#' respondent to the items on the given scale.  If \code{FALSE} (the
#' default), these variables are omitted from the returned data frame.
#' 
#' 
#' @return The original data frame is returned (optionally with modified
#' items if \code{updateItems = TRUE}) with new variables corresponding to
#' the scored scales. If \code{keepNvalid = TRUE}, for each scored scale an
#' additional variable is returned that contains the number of valid
#' responses each respondent made to the items making up the given scale.
#' These optional variables have names of the format \code{SCALENAME_N}.
#' The following scale scores are returned:
#' 
#' \describe{
#'   \item{PWB}{Physical Well-Being subscale}
#'   \item{SWB}{Social/Family Well-Being subscale}
#'   \item{EWB}{Emotional Well-Being subscale}
#'   \item{FWB}{Physical Well-Being subscale}
#'   \item{FACT_GP}{FACT-GP Total Score (i.e., PWB+SWB+EWB+FWB)}
#' }
#' 
#' @references FACT-GP Scoring Guidelines, available at \url{http://www.facit.org}
#' 
#' @export
#' 
#' @examples
#' set.seed(6375309)
#' ## Generating random item responses for 8 fake respondents
#' exampleDat <- t(replicate(8, sample(0:4, size = 27, replace = TRUE)))
#' ## Making half of respondents missing about 10% of items, 
#' ## half missing about 50%.
#' miss10 <- t(replicate(4, sample(c(0, 9), prob = c(0.9, 0.1), 
#'     size = 27, replace = TRUE)))
#' miss50 <- t(replicate(4, sample(c(0, 9), prob = c(0.5, 0.5), 
#'     size = 27, replace = TRUE)))
#' missMtx <- rbind(miss10, miss50)
#' ## Using 9 as the code for missing responses
#' exampleDat[missMtx == 9] <- 9
#' exampleDat <- as.data.frame(cbind(ID = paste0('ID', 1:8), 
#'     as.data.frame(exampleDat)))
#' names(exampleDat) <- c('ID', 'GP1', 'GP2', 'GP3', 'GP4', 'GP5',
#'     'GP6', 'GP7', 'GS1', 'GS2', 'GS3', 'GS4', 'GS5', 'GS6', 
#'     'GS7', 'GE1', 'GE2', 'GE3', 'GE4', 'GE5', 'GE6', 'GF1', 
#'     'GF2', 'GF3', 'GF4', 'GF5', 'GF6', 'GF7')
#' 
#' ## Returns data frame with scale scores and with original items untouched
#' scoredDat <- scoreFACT_GP(exampleDat)
#' names(scoredDat)
#' scoredDat
#' ## Returns data frame with scale scores, with the appropriate items
#' ## reverse scored, and with item values of 8 and 9 replaced with NA.
#' ## Also illustrates the effect of setting keepNvalid = TRUE.
#' scoredDat <- scoreFACT_GP(exampleDat, updateItems = TRUE, keepNvalid = TRUE)
#' names(scoredDat)
#' ## Descriptives of scored scales
#' summary(scoredDat[, c('PWB', 'SWB', 'EWB', 'FWB', 'FACT_GP')])
scoreFACT_GP <- function(df, updateItems = FALSE, keepNvalid = FALSE) {
    dfup <- df
    names(dfup) <- toupper(names(df))
    items <- dfup[, c("GP1", "GP2", "GP3", "GP4", "GP6", "GP7", "GS1", 
        "GS2", "GS3", "GS6", "GS7", "GE1", "GE4", "GE5", "GE6", "GF1", 
        "GF2", "GF3", "GF5", "GF6", "GF7")]
    if (any(!(as.matrix(items) %in% c(0:4, 8, 9, NA)))) {
        stop("At least 1 response is out of range (i.e., not 0-4, 8, 9, or NA)")
        break
    }
    makeMiss <- function(x) {
        x[x %in% c(8, 9)] <- NA
        return(x)
    }
    items <- as.data.frame(lapply(items, makeMiss))
    PWB_names <- c("GP1", "GP2", "GP3", "GP4", "GP6", "GP7")
    SWB_names <- c("GS1", "GS2", "GS3", "GS6", "GS7")
    EWB_names <- c("GE1", "GE4", "GE5", "GE6")
    FWB_names <- c("GF1", "GF2", "GF3", "GF5", "GF6", "GF7")
    revNames <- c("GP1", "GP2", "GP3", "GP4", "GP6", "GP7", "GE1", "GE4", 
        "GE5", "GE6")
    revHelper <- function(x) {
        return(4 - x)
    }
    items[, revNames] <- lapply(items[, revNames], revHelper)
    valid_N <- as.data.frame(lapply(items, function(x) as.numeric(!is.na(x))))
    PWB_N <- rowSums(valid_N[, PWB_names])
    SWB_N <- rowSums(valid_N[, SWB_names])
    EWB_N <- rowSums(valid_N[, EWB_names])
    FWB_N <- rowSums(valid_N[, FWB_names])
    FACT_GP_N <- PWB_N + SWB_N + EWB_N + FWB_N
    PWB <- round(rowMeans(items[, PWB_names], na.rm = TRUE) * length(PWB_names), 
        3)
    PWB[PWB_N/length(PWB_names) <= 0.5] <- NA
    SWB <- round(rowMeans(items[, SWB_names], na.rm = TRUE) * length(SWB_names), 
        3)
    SWB[SWB_N/length(SWB_names) <= 0.5] <- NA
    EWB <- round(rowMeans(items[, EWB_names], na.rm = TRUE) * length(EWB_names), 
        3)
    EWB[EWB_N/length(EWB_names) <= 0.5] <- NA
    FWB <- round(rowMeans(items[, FWB_names], na.rm = TRUE) * length(FWB_names), 
        3)
    FWB[FWB_N/length(FWB_names) <= 0.5] <- NA
    FACT_GP <- round(PWB + SWB + EWB + FWB, 3)
    FACT_GP[FACT_GP_N/ncol(items) <= 0.8] <- NA
    if (updateItems) {
        dfItemPos <- unlist(sapply(names(items), function(x) grep(x, names(df), 
            ignore.case = TRUE, value = FALSE)))
        names(df)[dfItemPos] <- toupper(names(df)[dfItemPos])
        df[, names(items)] <- items
    }
    if (keepNvalid) {
        dfOut <- as.data.frame(cbind(df, PWB_N, SWB_N, EWB_N, FWB_N, FACT_GP_N, 
            PWB, SWB, EWB, FWB, FACT_GP))
    } else {
        dfOut <- as.data.frame(cbind(df, PWB, SWB, EWB, FWB, FACT_GP))
    }
    return(dfOut)
} 
