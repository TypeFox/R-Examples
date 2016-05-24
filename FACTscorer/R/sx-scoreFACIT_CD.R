#' @title Score the FACIT-CD
#'
#' @description
#' Generates all of the scores of the Functional Assessment of Chronic Illness
#' Therapy - Cervical Dysplasia (FACIT_CD, v4) from item responses.
#'
#' @details
#' Given a data frame that includes all of the FACIT-CD (Version 4) items as
#' variables, appropriately named, this function generates all of the FACIT-CD
#' scale scores.  It is crucial that the item variables in the supplied data
#' frame are named according to FACT conventions.  For example, the first
#' physical well-being item should be named GP1, the second GP2, and so on.
#' Please refer to the materials provided by \url{http://www.facit.org} for the
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
#' data.  As its default, \code{scoreFACIT_CD} DOES NOT replace any of your
#' original item variables with the reverse coded versions.  However, for
#' consistentcy with the behavior of the other versions on
#' \url{http://www.facit.org}, the \code{updateItems} argument is
#' provided.  If set to \code{TRUE}, any item that is supposed to be
#' reverse coded will be replaced with its reversed version in the data
#' frame returned by \code{scoreFACIT_CD}.
#'
#'
#' @param df A data frame with the FACIT-CD items, appropriately-named.
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
#'   \item{PWB}{Physical Well-Being subscale}
#'  \item{TS}{Treatmnet Satisfation subscale}
#'   \item{GP}{General Perceptions subscale}
#'   \item{EWB}{Emotional Well-Being subscale}
#'   \item{REL}{Relationships subscale}
#'   \item{FACIT_CD_TOTAL}{FACIT-CD Total Score (i.e., PWB+TS+GP+EWB+REL)}
#' }
#'
#' @references FACIT-CD Scoring Guidelines, available at \url{http://www.facit.org}
#'
#' @export
#'
#' @examples
#' ## Setting up item names for fake data
#' AC_names1 <- c('CD1', 'CD2', 'CD3', 'Cx1', 'GP5', 'ES8', 'CD4', 'CD5')
#' AC_names2 <- c('GR1', 'CD6', 'CD7', 'CD8')
#' AC_names3 <- c('GF1', 'GF3', 'HI11', 'Sp9', 'GF7', 'CD9', 'CD10')
#' AC_names4 <- c('CD11', 'CD12', 'CD13', 'BMT18', 'CD14', 'CD15', 'CD16', 'CD17',
#'                'CD18', 'CD19', 'CD20')
#' AC_names5 <- c('CD21', 'CD22', 'GS1', 'HI3')
#' AC_names <- c(AC_names1, AC_names2, AC_names3, AC_names4, AC_names5)
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
#' scoredDat <- scoreFACIT_CD(exampleDat)
#' names(scoredDat)
#' scoredDat
#' ## Returns data frame with scale scores, with the appropriate items
#' ## reverse scored, and with item values of 8 and 9 replaced with NA.
#' ## Also illustrates the effect of setting keepNvalid = TRUE.
#' scoredDat <- scoreFACIT_CD(exampleDat, updateItems = TRUE, keepNvalid = TRUE)
#' names(scoredDat)
#' ## Descriptives of scored scales
#' summary(scoredDat[, c('PWB', 'TS', 'GP', 'EWB', 'REL', 'FACIT_CD_TOTAL')])
scoreFACIT_CD <- function(df, updateItems = FALSE, keepNvalid = FALSE) {
    dfup <- df
    names(dfup) <- toupper(names(df))
    AC_names1 <- toupper(c("CD1", "CD2", "CD3", "Cx1", "GP5", "ES8", "CD4",
        "CD5"))
    AC_names2 <- toupper(c("GR1", "CD6", "CD7", "CD8"))
    AC_names3 <- toupper(c("GF1", "GF3", "HI11", "Sp9", "GF7", "CD9", "CD10"))
    AC_names4 <- toupper(c("CD11", "CD12", "CD13", "BMT18", "CD14", "CD15",
        "CD16", "CD17", "CD18", "CD19", "CD20"))
    AC_names5 <- toupper(c("CD21", "CD22", "GS1", "HI3"))
    AC_names <- c(AC_names1, AC_names2, AC_names3, AC_names4, AC_names5)
    revNames <- toupper(c("CD1", "CD2", "CD3", "Cx1", "GP5", "ES8", "CD4",
        "CD5", "CD11", "CD12", "CD13", "BMT18", "CD14", "CD15", "CD16",
        "CD17", "CD18", "CD19", "CD20"))
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
    AC_N5 <- rowSums(valid_N[, AC_names5])
    TOTAL_N <- AC_N1 + AC_N2 + AC_N3 + AC_N4 + AC_N5
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
    AC5 <- round(rowMeans(AC_items[, AC_names5], na.rm = TRUE) * length(AC_names5),
        3)
    AC5[AC_N5/length(AC_names5) <= 0.5] <- NA
    TOTAL <- AC1 + AC2 + AC3 + AC4 + AC5
    TOTAL[TOTAL_N/length(AC_names) <= 0.8] <- NA
    FACIT_CD_TOTAL_N <- TOTAL_N
    PWB_N <- AC_N1
    TS_N <- AC_N2
    GP_N <- AC_N3
    EWB_N <- AC_N4
    REL_N <- AC_N5
    PWB <- AC1
    TS <- AC2
    GP <- AC3
    EWB <- AC4
    REL <- AC5
    FACIT_CD_TOTAL <- TOTAL
    if (updateItems) {
        dfItemPos <- unlist(sapply(AC_names, function(x) grep(x, names(df),
            ignore.case = TRUE, value = FALSE)))
        names(df)[dfItemPos] <- toupper(names(df)[dfItemPos])
        df[, AC_names] <- AC_items
    }
    if (keepNvalid) {
        dfOut <- as.data.frame(cbind(df, PWB_N, TS_N, GP_N, EWB_N, REL_N,
            FACIT_CD_TOTAL_N, PWB, TS, GP, EWB, REL, FACIT_CD_TOTAL))
    } else {
        dfOut <- as.data.frame(cbind(df, PWB, TS, GP, EWB, REL, FACIT_CD_TOTAL))
    }
    return(dfOut)
}
