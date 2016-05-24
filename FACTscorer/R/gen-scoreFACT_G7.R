#' @title Score the FACT-G7
#' 
#' @description
#' Generates all of the scores of the 7-item short-form of the Functional 
#' Assessment of Cancer Therapy-General (FACT-G7) from item responses.
#' 
#' @details
#' Given a data frame that includes all of the FACT-G7 items as variables,
#' appropriately named, this function generates all of the FACT-G7 scale
#' scores.  It is crucial that the item variables in the supplied data frame
#' are named according to FACT conventions.  For example, the first item should be 
#' named GP1, the second GP2, and so on.  Please refer to the materials provided 
#' by \url{http://www.facit.org} for the particular questionnaire you are using.  
#' In particular, refer to the left margin of the official questionnaire (i.e., 
#' from facit.org) for the appropriate item variable names.
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
#' data.  As its default, \code{scoreFACTG7} DOES NOT replace any of your
#' original item variables with the reverse coded versions.  However, for
#' consistentcy with the behavior of the other versions on
#' \url{http://www.facit.org}, the \code{updateItems} argument is
#' provided.  If set to \code{TRUE}, any item that is supposed to be
#' reverse coded will be replaced with its reversed version in the data
#' frame returned by \code{scoreFACTG7}.
#' 
#' 
#' @param df A data frame with the FACT-G7 items, appropriately-named.
#' @param updateItems Logical, if \code{TRUE} any original item that is
#' reverse coded for scoring will be replaced by its reverse coded version
#' in the returned data frame, and any values of 8 or 9 will be replaced
#' with NA.  The default, \code{FALSE}, returns the original items
#' unmodified.  
#' @param keepNvalid Logical, if \code{TRUE} the function
#' returns an additional variable for the returned scale scores
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
#'   \item{FACTG7}{FACT-G7 Total Score}
#' }
#' 
#' @references FACT-G7 Scoring Guidelines, available at \url{http://www.facit.org}
#' 
#' @export
#' 
#' @examples
#' ## Setting up item names for fake data
#' itemNames <- c('GP1', 'GP4', 'GP2', 'GE6', 'GF5', 'GF3', 'GF7')
#' set.seed(6375309)
#' ## Generating random item responses for 8 fake respondents
#' exampleDat <- t(replicate(8, sample(0:4, size = length(itemNames), replace = TRUE)))
#' ## Making half of respondents missing about 10% of items, 
#' ## half missing about 50%.
#' miss10 <- t(replicate(4, sample(c(0, 9), prob = c(0.9, 0.1), 
#'     size = length(itemNames), replace = TRUE)))
#' miss50 <- t(replicate(4, sample(c(0, 9), prob = c(0.5, 0.5), 
#'     size = length(itemNames), replace = TRUE)))
#' missMtx <- rbind(miss10, miss50)
#' ## Using 9 as the code for missing responses
#' exampleDat[missMtx == 9] <- 9
#' exampleDat <- as.data.frame(cbind(ID = paste0('ID', 1:8), 
#'     as.data.frame(exampleDat)))
#' names(exampleDat) <- c('ID', itemNames)
#' 
#' ## Returns data frame with scale scores and with original items untouched
#' scoredDat <- scoreFACTG7(exampleDat)
#' names(scoredDat)
#' scoredDat
#' ## Returns data frame with scale scores, with the appropriate items
#' ## reverse scored, and with item values of 8 and 9 replaced with NA.
#' ## Also illustrates the effect of setting keepNvalid = TRUE.
#' scoredDat <- scoreFACTG7(exampleDat, updateItems = TRUE, keepNvalid = TRUE)
#' names(scoredDat)
#' ## Descriptives of scored scales
#' summary(scoredDat[, c('FACTG7')])
scoreFACTG7 <- function(df, updateItems = FALSE, keepNvalid = FALSE) {
    dfup <- df
    names(dfup) <- toupper(names(df))
    items <- dfup[, c("GP1", "GP4", "GP2", "GE6", "GF5", "GF3", "GF7")]
    if (any(!(as.matrix(items) %in% c(0:4, 8, 9, NA)))) {
        stop("At least 1 response is out of range (i.e., not 0-4, 8, 9, or NA)")
        break
    }
    makeMiss <- function(x) {
        x[x %in% c(8, 9)] <- NA
        return(x)
    }
    items <- as.data.frame(lapply(items, makeMiss))
    PWB_names <- c("GP1", "GP4", "GP2", "GE6", "GF5", "GF3", "GF7")
    revNames <- c("GP1", "GP4", "GP2", "GE6")
    revHelper <- function(x) {
        return(4 - x)
    }
    items[, revNames] <- lapply(items[, revNames], revHelper)
    valid_N <- as.data.frame(lapply(items, function(x) as.numeric(!is.na(x))))
    PWB_N <- rowSums(valid_N[, PWB_names])
    FACTG7_N <- PWB_N
    PWB <- round(rowMeans(items[, PWB_names], na.rm = TRUE) * length(PWB_names), 
        3)
    PWB[PWB_N/length(PWB_names) <= 0.5] <- NA
    FACTG7 <- round(PWB, 3)
    FACTG7[FACTG7_N/ncol(items) <= 0.8] <- NA
    if (updateItems) {
        dfItemPos <- unlist(sapply(names(items), function(x) grep(x, names(df), 
            ignore.case = TRUE, value = FALSE)))
        names(df)[dfItemPos] <- toupper(names(df)[dfItemPos])
        df[, names(items)] <- items
    }
    if (keepNvalid) {
        dfOut <- as.data.frame(cbind(df, FACTG7_N, FACTG7))
    } else {
        dfOut <- as.data.frame(cbind(df, FACTG7))
    }
    return(dfOut)
} 
