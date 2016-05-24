# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

# TODO: error handling

#' Equivalized household size
#' 
#' Compute the equivalized household size according to the modified OECD scale
#' adopted in 1994.
#' 
#' @param hid if \code{data=NULL}, a vector containing the household ID.
#' Otherwise a character string specifying the column of \code{data} that
#' contains the household ID.
#' @param age if \code{data=NULL}, a vector containing the age of the
#' individuals.  Otherwise a character string specifying the column of
#' \code{data} that contains the age.
#' @param year if \code{data=NULL}, a vector containing the year of the survey.
#' Otherwise a character string specifying the column of \code{data} that
#' contains the year.
#' @param data a \code{data.frame} containing EU-SILC survey data, or
#' \code{NULL}.
#' 
#' @return A numeric vector containing the equivalized household size for every
#' observation in \code{data}.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{eqInc}}
#' 
#' @references Working group on Statistics on Income and Living Conditions
#' (2004) Common cross-sectional EU indicators based on EU-SILC; the gender pay
#' gap.  \emph{EU-SILC 131-rev/04}, Eurostat.
#' 
#' @keywords survey
#' 
#' @examples
#' data(eusilc)
#' 
#' # calculate equivalized household size
#' eqSS <- eqSS("db030", "age", data=eusilc)
#' 
#' # combine with household ID and household size
#' tmp <- cbind(eusilc[, c("db030", "hsize")], eqSS)
#' 
#' # show the first 8 rows
#' head(tmp, 8)
#' 
#' @export

eqSS <- function(hid, age, year = NULL, data = NULL) {
    ## initializations
    if(is.null(data)) {
        data <- data.frame(hid=hid)
        hid <- "hid"
        if(!is.null(year)) {
            data <- cbind(year=year, data)
            year <- "year"
        }
    } else {
        age <- data[, age]
        data <- data[, c(year, hid), drop=FALSE]
    }
    ## calculations
    i <- if(is.null(year)) 2 else 3
    tmp <- as.data.frame(table(data))  # number of household members
    hm14p <- as.data.frame(table(data[age >= 14,]))[, i]  # at least 14 years
    hm13m <- tmp[, i] - hm14p  # younger than 14
    tmp[, i] <- 1 + 0.5*(hm14p-1) + 0.3*hm13m  # eqSS for househoulds
    names(tmp) <- c(year, hid, ".eqSS")
    data <- cbind(data, .ID=1:nrow(data))  # add ID to original data
    data <- merge(data, tmp, sort=FALSE)  # merge with original data set
    ## order according to original data and extract eqSS
    data$.eqSS[order(data$.ID)]
}
