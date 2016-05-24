# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

# TODO: error handling
# TODO: account for inflation-adjustment

#' Equivalized disposable income
#' 
#' Compute the equivalized disposable income from household and personal income
#' variables.
#' 
#' All income components should already be imputed, otherwise \code{NA}s are
#' simply removed before the calculations.
#' 
#' @param hid if \code{data=NULL}, a vector containing the household ID.
#' Otherwise a character string specifying the column of \code{data} that
#' contains the household ID.
#' @param hplus if \code{data=NULL}, a \code{data.frame} containing the
#' household income components that have to be added.  Otherwise a character
#' vector specifying the columns of \code{data} that contain these income
#' components.
#' @param hminus if \code{data=NULL}, a \code{data.frame} containing the
#' household income components that have to be subtracted.  Otherwise a
#' character vector specifying the columns of \code{data} that contain these
#' income components.
#' @param pplus if \code{data=NULL}, a \code{data.frame} containing the personal
#' income components that have to be added.  Otherwise a character vector
#' specifying the columns of \code{data} that contain these income components.
#' @param pminus if \code{data=NULL}, a \code{data.frame} containing the
#' personal income components that have to be subtracted.  Otherwise a character
#' vector specifying the columns of \code{data} that contain these income
#' components.
#' @param eqSS if \code{data=NULL}, a vector containing the equivalized
#' household size.  Otherwise a character string specifying the column of
#' \code{data} that contains the equivalized household size.  See
#' \code{\link{eqSS}} for more details.
#' @param year if \code{data=NULL}, a vector containing the year of the survey.
#' Otherwise a character string specifying the column of \code{data} that
#' contains the year.
#' @param data a \code{data.frame} containing EU-SILC survey data, or
#' \code{NULL}.
#' 
#' @return A numeric vector containing the equivalized disposable income for
#' every individual in \code{data}.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{eqSS}}
#' 
#' @references Working group on Statistics on Income and Living Conditions
#' (2004) Common cross-sectional EU indicators based on EU-SILC; the gender pay
#' gap.  \emph{EU-SILC 131-rev/04}, Eurostat.
#' @keywords survey
#' 
#' @examples
#' data(eusilc)
#' 
#' # compute a simplified version of the equivalized disposable income
#' # (not all income components are available in the synthetic data)
#' hplus <- c("hy040n", "hy050n", "hy070n", "hy080n", "hy090n", "hy110n")
#' hminus <- c("hy130n", "hy145n")
#' pplus <- c("py010n", "py050n", "py090n", "py100n", 
#'     "py110n", "py120n", "py130n", "py140n")
#' eqIncome <- eqInc("db030", hplus, hminus, 
#'     pplus, character(), "eqSS", data=eusilc)
#' 
#' # combine with household ID and equivalized household size
#' tmp <- cbind(eusilc[, c("db030", "eqSS")], eqIncome)
#' 
#' # show the first 8 rows
#' head(tmp, 8)
#' 
#' @export

eqInc <- function(hid, hplus, hminus, pplus, pminus, 
        eqSS, year = NULL, data = NULL) {
    ## initializations
    if(is.null(data)) {
        data <- data.frame(hid=hid)
        hid <- "hid"
        if(!is.null(year)) {
            data <- cbind(year=year, data)
            year <- "year"
        }
        npplus <- names(pplus)
        npminus <- names(pminus)
    } else {
        hplus <- data[, hplus]
        hminus <- data[, hminus]
        npplus <- pplus
        pplus <- data[, npplus]
        npminus <- pminus
        pminus <- data[, npminus]
        eqSS <- data[, eqSS]
        data <- data[, c(year, hid), drop=FALSE]
    }
    ## calculations
    hy020h <- rowSums(hplus, na.rm=TRUE) - rowSums(hminus, na.rm=TRUE)
    tmp <- aggregate(data.frame(pplus,pminus), data, sum, na.rm=TRUE)
    hy020p <- rowSums(tmp[,npplus], na.rm=TRUE) - 
        rowSums(tmp[,npminus], na.rm=TRUE)
    if(is.null(year)) {
        names(hy020p) <- tmp[, hid]
        hy020p <- unname(hy020p[as.character(data[, hid])])
    } else {
        tmp <- cbind(tmp[, c(year, hid), drop=FALSE], .hy020p=hy020p)
        data <- cbind(data, .ID=1:nrow(data))  # add ID to original data
        data <- merge(data, tmp, sort=FALSE)  # merge with original data set
        ## order according to original data and extract hy020p
        hy020p <- data$.hy020p[order(data$.ID)]
    }
    ## return result
    (hy020h + hy020p) / eqSS
}
