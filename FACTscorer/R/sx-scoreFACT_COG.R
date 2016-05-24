#' @title Score the FACT-COG
#'
#' @description
#' Generates all of the scores of the Functional Assessment of Cancer Therapy -
#' Cognitive function issues - Patient Satisfaction (FACT-COG, v4) from item responses.
#'
#' @details
#' Given a data frame that includes all of the FACT-COG (Version 4) items as
#' variables, appropriately named, this function generates all of the FACT-COG
#' scale scores.  It is crucial that the item variables in the supplied data
#' frame are named according to FACT conventions.  For example, the first
#' Perceived Cognitive Impairments item should be named CogA1, the second CogA3,
#' and so on. Please refer to the materials provided by \url{http://www.facit.org}
#' for the particular questionnaire you are using.  In particular, refer to the
#' left margin of the official questionnaire (i.e., from facit.org) for the
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
#' data.  As its default, \code{scoreFACT_COG} DOES NOT replace any of your
#' original item variables with the reverse coded versions.  However, for
#' consistentcy with the behavior of the other versions on
#' \url{http://www.facit.org}, the \code{updateItems} argument is
#' provided.  If set to \code{TRUE}, any item that is supposed to be
#' reverse coded will be replaced with its reversed version in the data
#' frame returned by \code{scoreFACT_COG}.
#'
#'
#' @param df A data frame with the FACT-COG items, appropriately-named.
#' @param updateItems Logical, if \code{TRUE} any original item that is
#' reverse coded for scoring will be replaced by its reverse coded version
#' in the returned data frame, and any values of 8 or 9 will be replaced
#' with NA.  The default, \code{FALSE}, returns the original items
#' unmodified.
#' @param keepNvalid Logical, if \code{TRUE} the function
#' returns an additional variable for each of the returned scale scores
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
#'   \item{CogPCI}{Perceived Cognitive Impairments subscale}
#'  \item{CogQOL}{Impact of perceived cognitive impairments on quality of life subscale}
#'   \item{CogOth}{Comments from Others subscale}
#'   \item{CogPCA}{Perceived Cognitive Abilities subscale}
#' }
#'
#' @references FACT-COG Scoring Guidelines, available at \url{http://www.facit.org}
#'
#' @export
#'
#' @examples
#' ## Setting up item names for fake data
#' AC_names1 <- c('CogA1', 'CogA3', 'CogC7', 'CogM9', 'CogM10', 'CogM12', 'CogV13',
#'                'CogV15', 'CogV16', 'CogV17b', 'CogF19', 'CogF23', 'CogF24', 'CogF25',
#'                'CogC31', 'CogC32', 'CogC33a', 'CogC33c')
#' AC_names2 <- c('CogQ35', 'CogQ37', 'CogQ38', 'CogQ41')
#' AC_names3 <- c('CogO1', 'CogO2', 'CogO3', 'CogO4')
#' AC_names4 <- c('CogPC1', 'CogPV1', 'CogPM1','CogPM2', 'CogPF1', 'CogPCh1', 'CogPCh2')
#' AC_names <- c(AC_names1, AC_names2, AC_names3, AC_names4)
#' itemNames <- AC_names
#' ## Generating random item responses for 8 fake respondents
#' set.seed(6375309)
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
#' scoredDat <- scoreFACT_COG(exampleDat)
#' names(scoredDat)
#' scoredDat
#' ## Returns data frame with scale scores, with the appropriate items
#' ## reverse scored, and with item values of 8 and 9 replaced with NA.
#' ## Also illustrates the effect of setting keepNvalid = TRUE.
#' scoredDat <- scoreFACT_COG(exampleDat, updateItems = TRUE, keepNvalid = TRUE)
#' names(scoredDat)
#' ## Descriptives of scored scales
#' summary(scoredDat[, c('CogPCI', 'CogQOL', 'CogOth', 'CogPCA')])
scoreFACT_COG <- function(df, updateItems = FALSE, keepNvalid = FALSE) {
    dfup <- df
    names(dfup) <- toupper(names(dfup))
    AC_names1 <- toupper(c("CogA1", "CogA3", "CogC7", "CogM9", "CogM10",
        "CogM12", "CogV13", "CogV15", "CogV16", "CogV17b", "CogF19", "CogF23",
        "CogF24", "CogF25", "CogC31", "CogC32", "CogC33a", "CogC33c"))
    AC_names2 <- toupper(c("CogQ35", "CogQ37", "CogQ38", "CogQ41"))
    AC_names3 <- toupper(c("CogO1", "CogO2", "CogO3", "CogO4"))
    AC_names4 <- toupper(c("CogPC1", "CogPV1", "CogPM1", "CogPM2", "CogPF1",
        "CogPCh1", "CogPCh2"))
    AC_names <- c(AC_names1, AC_names2, AC_names3, AC_names4)
    revNames <- toupper(c("CogA1", "CogA3", "CogC7", "CogM9", "CogM10",
        "CogM12", "CogV13", "CogV15", "CogV16", "CogV17b", "CogF19", "CogF23",
        "CogF24", "CogF25", "CogC31", "CogC32", "CogC33a", "CogC33c", "CogQ35",
        "CogQ37", "CogQ38", "CogQ41", "CogO1", "CogO2", "CogO3", "CogO4"))
    AC_items <- dfup[, AC_names]
    if (any(!(as.matrix(AC_items) %in% c(0:4, 8, 9, NA)))) {
        stop("At least 1 response is out of range (i.e., not 0-4, 8, 9, or NA)")
        break
    }
    makeMiss <- function(x) {
        x[x %in% c(8, 9)] <- NA
        return(x)
    }
    AC_items <- as.data.frame(lapply(AC_items, makeMiss))
    revHelper <- function(x) {
        return(4 - x)
    }
    AC_items[, revNames] <- lapply(AC_items[, revNames], revHelper)
    valid_N <- as.data.frame(lapply(AC_items, function(x) as.numeric(!is.na(x))))
    AC_N1 <- rowSums(valid_N[, AC_names1])
    AC_N2 <- rowSums(valid_N[, AC_names2])
    AC_N3 <- rowSums(valid_N[, AC_names3])
    AC_N4 <- rowSums(valid_N[, AC_names4])
    TOTAL_N <- AC_N1 + AC_N2 + AC_N3 + AC_N4
    AC1 <- round(rowMeans(AC_items[, AC_names1], na.rm = TRUE) * length(AC_names1),
        3)
    AC1[AC_N1/length(AC_names1) <= 0.5] <- NA
    AC2 <- round(rowMeans(AC_items[, AC_names2], na.rm = TRUE) * length(AC_names2),
        3)
    AC2[AC_N2/length(AC_names2) <= 0.5] <- NA
    AC3 <- round(rowMeans(AC_items[, AC_names3], na.rm = TRUE) * length(AC_names3),
        3)
    AC3[AC_N3/length(AC_names3) <= 0.5] <- NA
    AC4 <- round(rowMeans(AC_items[, AC_names4], na.rm = TRUE) * length(AC_names4),
        3)
    AC4[AC_N4/length(AC_names4) <= 0.5] <- NA
    CogPCI_N <- AC_N1
    CogQOL_N <- AC_N2
    CogOth_N <- AC_N3
    CogPCA_N <- AC_N4
    CogPCI <- AC1
    CogQOL <- AC2
    CogOth <- AC3
    CogPCA <- AC4
    if (updateItems) {
        dfItemPos <- unlist(sapply(AC_names, function(x) grep(x, names(df),
            ignore.case = TRUE, value = FALSE)))
        names(df)[dfItemPos] <- toupper(names(df)[dfItemPos])
        df[, AC_names] <- AC_items
    }
    if (keepNvalid) {
        dfOut <- as.data.frame(cbind(df, CogPCI_N, CogQOL_N, CogOth_N,
            CogPCA_N, CogPCI, CogQOL, CogOth, CogPCA))
    } else {
        dfOut <- as.data.frame(cbind(df, CogPCI, CogQOL, CogOth, CogPCA))
    }
    return(dfOut)
}
