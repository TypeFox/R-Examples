#'
#' @title Functions for Working with the 'www.estimize.com' Web Services
#'
#' @description Provides the user with functions to develop their trading
#' strategy, uncover actionable trading ideas, and monitor consensus shifts with
#' crowdsourced earnings and economic estimate data directly from
#' <www.estimize.com>. Further information regarding the web services this
#' package invokes can be found at <www.estimize.com/api>.
#'
#' @details Provides a complete implementation of all web methods available at
#' at \href{http://api.estimize.com}{api.estimize.com} -- functions available
#' include:
#'     \itemize{
#'     \item \code{\link{GetCompanies}}
#'     \item \code{\link{GetCompany}}
#'     \item \code{\link{GetCompanyEstimates}}
#'     \item \code{\link{GetCompanyReleases}}
#'     \item \code{\link{GetEstimates}}
#'     \item \code{\link{GetReleaseConsensus}}
#' }
#'
#' Note that the \code{\link{SetKey}} function must be called \emph{once} prior
#' to invoking any of the other functions in this package.
#'
#' @import RJSONIO
#' @import RCurl
#'
#' @docType package
#'
#' @name restimizeapi
#'
NULL

.rEstimizeAPI.env <- new.env()

#' The api.estimize.com url, which is used elsewhere in this package
#' when invoking various web methods.
#'
.baseURL <- "http://api.estimize.com/"

NAME <- "name"

TICKER <- "ticker"
TICKER_KEY <- TICKER

COMPANY_NAME <- "Company Name"

FISCAL_YEAR <- "Fiscal Year"
FISCAL_YEAR_KEY <- "fiscal_year"

FISCAL_QUARTER <- "Fiscal Quarter"
FISCAL_QUARTER_KEY <- "fiscal_quarter"

EPS <- "EPS"
EPS_KEY <- "eps"

EPS_REVISION = "EPS Revision"

EPS_REVENUE <- "EPS Revenue"
REVENUE_KEY <- "revenue"

CONSENSUS_EPS_ESTIMATE <- "Consensus EPS Estimate"
CONSENSUS_EPS_ESTIMATE_KEY <- "consensus_eps_estimate"

CONSENSUS_REVENUE_ESTIMATE <- "Consensus Revenue Estimate"
CONSENSUS_REVENUE_ESTIMATE_KEY <- "consensus_revenue_estimate"

WALLSTREET_REVENUE_ESTIMATE <- "Wall Street Revenue Estimate"
WALLSTREET_REVENUE_ESTIMATE_KEY <- "wallstreet_revenue_estimate"

WALLSTREET_EPS_ESTIMATE <- "Wall Street EPS Estimate"
WALLSTREET_EPS_ESTIMATE_KEY <- "wallstreet_eps_estimate"

RELEASE_DATE <- "Release Date"
RELEASE_DATE_KEY <- "release_date"

ID <- "Id"
ID_KEY <- "id"

USERNAME <- "User Name"
USERNAME_KEY <- "username"

CREATED_AT <- "Created At"
CREATED_AT_KEY <- "created_at"

ANALYST_ID = "Analyst Id"
ANALYST_ID_KEY <- "analyst_id"

MEAN_KEY <- "mean"
MEAN <- "Mean"

HIGH_KEY <- "high"
HIGH <- "High"

LOW_KEY <- "low"
LOW <- "Low"

STANDARD_DEVIATION_KEY <- "standard_deviation"
STANDARD_DEVIATION <- "Standard Deviation"

COUNT_KEY <- "count"
COUNT <- "Count"

UPDATED_AT_KEY <- "updated_at"
UPDATED_AT <- "Updated At"

WALL_STREET <- "Wall Street"
ESTIMIZE <- "Estimize"
EPS <- "EPS"
REVENUE <- "Revenue"
REVENUE_REVISION <- "Revenue Revision"

#' Function halts the execution of the script if the paramValue is null.
#'
#' @param paramName The parameter name that is being checked for a null value.
#'
#' @param paramValue the parameter value which cannot be null.
#'
.assertNotNull <- function (paramName, paramValue) {
    if (is.null (paramValue)) {
        stop (paste ("The value of the parameter ", paramName,
            " cannot be  null."))
    }
}

#' Function halts the script if the ticker, year, or quarter is null.
#'
#' @param ticker The ticker should always be a non-NULL value.
#' @param The year when non-NULL should also be a numeric.
#' @param The quarter when non-NULL should also be a numeric between one and
#' four inclusive.
#'
.CheckParams <- function (ticker, year, quarter) {
    
    if (!is.null (year) && !is.numeric (year))
        stop ("The year must be a number")
    
    if (!is.null (quarter)
        && !is.numeric (quarter))
        stop ("The quarter must be a number between 1 and 4 inclusive.")
    else if (!is.null (quarter)
        && is.numeric (quarter)
        && !((1 <= quarter) & (quarter <= 4)))
        stop ("The quarter must be a number between 1 and 4 inclusive.")
}

#' An action which does nothing and simply returns the data.
#'
#' @param data The data.
#'
#' @return The data.
#'
.emptyAction <- function (data) {
    return (data)
}

#' Function performs the request for data from the estimize.com website and
#' returns the result as a data frame.
#'
#' @param destination The Estimize API web method to invoke.
#'
#' @param params Any parameter pertaining to a specific web method invocation.
#'
#' @param andExecute A function that is be applied to the resultant data. This
#'  function takes a data param which is the result returned from the call to
#'  fromJSON. The function's return value is returned to the caller.
#'
.DoGet <- function (destination, params = NULL, andExecute = .emptyAction) {

    key <- .rEstimizeAPI.env$key

    message <- paste ("The key has not been set which must be done prior to ",
        "calling any of the other functions this package provides. Refer to ",
        "the SetKey function and contact sales@estimize.com to request a ",
        "key if you need one.",
        sep="")

    if (is.null (key)) {
        stop (message)
    }

    header <- list ("X-Estimize-Key" = key)

    data <- getURL (url=destination, httpheader=header)

    tryCatch (
        unconvertedResult <- RJSONIO::fromJSON(data, nullValue = "null"),
        error = function (exception) {
            stop (
                paste (
                    "An exception was thrown while attempting to transform ",
                    "the JSON returned from the destination URL: ",
                    destination,
                    exception,
                    ".",
                    sep="\n"))
        })

    if (is.null (andExecute)) {
        warning ("The action function is null so .emptyAction will be substituted.")
        andExecute <- .emptyAction
    }

    result <- andExecute (unconvertedResult)

    return (result)
}

#' Function sets the key which is required as a header in every web method call
#' to estimize.com -- contact \email{sales@@estimize.com} to obtain a key.
#'
#' @param key The API key as provided by Estimize.com.
#'
#' @examples
#'  \dontrun{
#'  SetKey ("ENTER YOUR PRIVATE API KEY HERE.")
#'  }
#'
#' @export
#'
SetKey <- function (key) {
    .rEstimizeAPI.env$key <- key
}

#' Function transforms the companies into a data frame.
#'
#' @param companies A list of companies.
#'
#' @return The companies as a data frame.
#'
.TransformCompanies <- function (companies) {

    names <- unlist (lapply (companies, '[', NAME))
    tickers <- unlist (lapply (companies, '[', TICKER))

    result <- data.frame (names, tickers)

    names (result) <- c(NAME, TICKER)

    return (result)
}

#' Returns the company information for all companies that Estimize covers.
#'
#' @return The companies information as a data frame.
#'
#' @examples
#'  \dontrun{
#'  result <- GetCompanies ()
#'  }
#'
#' @export
#'
GetCompanies <- function () {

    destination <- paste (.baseURL, "companies.json", sep="")

    result <- .DoGet (destination, andExecute = .TransformCompanies)

    return (result)
}

#' Function transforms the company into a data frame.
#'
#' @param company A single company.
#'
#' @return The company information as a data frame.
#'
.TransformCompany <- function (company) {

    result <- data.frame (company['name'], company['ticker'])

    names (result) <- c(COMPANY_NAME, TICKER)

    return (result)
}

#' Returns the company information for the company specified by the ticker.
#'
#' @param ticker The Estimize ticker -- for example "MSFT".
#'
#' @return The company information as a data frame.
#'
#' @examples
#'  \dontrun{
#'      result <- GetCompany ("msft")
#'  }
#'
#' @export
#'
GetCompany <- function (ticker) {

    if (is.null (ticker))
        stop ("The ticker symbol is required (it is currently NULL).")

    destination <- paste (.baseURL, "companies/", ticker, ".json", sep="")

    result <- .DoGet (destination, andExecute = .TransformCompany)

    return (result)
}

#' Function transforms the company releases into a data frame.
#'
#' @param companyReleases A releases list.
#'
#' @return The company releases as a data frame.
#'
.TransformCompanyReleases <- function (companyReleases) {

    fiscalYears <- unlist (
        lapply (companyReleases, '[', FISCAL_YEAR_KEY))
    fiscalQuarters <- unlist (
        lapply (companyReleases, '[', FISCAL_QUARTER_KEY))
    epss <- unlist (
        lapply (companyReleases, '[', EPS_KEY))
    revenues <- unlist (
        lapply (companyReleases, '[', REVENUE_KEY))
    consensusEPSEstimates <- unlist (
        lapply (companyReleases, '[', CONSENSUS_EPS_ESTIMATE_KEY))
    consensusRevenueEstimates <- unlist (
        lapply (companyReleases, '[', CONSENSUS_REVENUE_ESTIMATE_KEY))
    wallStreetRevenueEstimates <- unlist (
        lapply (companyReleases, '[', WALLSTREET_REVENUE_ESTIMATE_KEY))
    wallStreetEPSEstimates <- unlist (
        lapply (companyReleases, '[', WALLSTREET_EPS_ESTIMATE_KEY))
    releaseDates <- unlist (
        lapply (companyReleases, '[', RELEASE_DATE_KEY))
    ids <- unlist (
        lapply (companyReleases, '[', ID_KEY))

    result <- data.frame (
        fiscalYears,
        fiscalQuarters,
        epss,
        revenues,
        consensusEPSEstimates,
        consensusRevenueEstimates,
        wallStreetRevenueEstimates,
        wallStreetEPSEstimates,
        releaseDates,
        ids)

    names (result) <- c(
        FISCAL_YEAR,
        FISCAL_QUARTER,
        EPS,
        REVENUE,
        CONSENSUS_EPS_ESTIMATE,
        CONSENSUS_REVENUE_ESTIMATE,
        WALLSTREET_REVENUE_ESTIMATE,
        WALLSTREET_EPS_ESTIMATE,
        RELEASE_DATE,
        ID)

    return (result)
}

#' Returns the past financial releases for the specified company by the ticker,
#' for the specified fiscal year, and quarter.
#'
#' @param ticker The Estimize ticker -- for example "MSFT".
#' @param year A four-digit year -- for example 1999.
#' @param quarter A numeric value between 1 and 4 inclusive; note that if the
#'  quarter is used the year must also be set as well.
#'
#' @return The company releases as a data frame.
#'
#' @examples
#'  \dontrun{
#'  result <- GetCompanyReleases ("MSFT")
#'  result <- GetCompanyReleases ("MSFT", 2009)
#'  result <- GetCompanyReleases ("MSFT", 2009, 2)
#'  }
#'
#' @export
#'
GetCompanyReleases <- function (
    ticker,
    year=NULL,
    quarter=NULL
) {

    .CheckParams (ticker, year, quarter)

    if (is.null (year) && is.null (quarter)) {
        destination <- paste (.baseURL, "companies/", ticker, "/releases.json", sep="")
    } else if (!is.null (year) && is.null (quarter)) {
        destination <- paste (.baseURL, "companies/", ticker, "/releases/", year, ".json", sep="")
    } else if (!is.null (year) && !is.null (quarter)) {
        destination <- paste (.baseURL, "companies/", ticker, "/releases/", year, "/", quarter, ".json", sep="")
    } else {
        stop ("This function requires a year and a quarter (a quarter alone is not sufficient).")
    }

    result <- .DoGet (destination, andExecute = .TransformCompanyReleases)
    
    return (result)
}

#' Function transforms the company estimates into a data frame.
#'
#' @param companyEstimates An estimates list.
#'
#' @return The company estimates as a data frame.
#'
.TransformCompanyEstimates <- function (companyEstimates) {

    ticker <- unlist (
        lapply (companyEstimates, '[', TICKER_KEY))
    fiscalYears <- unlist (
        lapply (companyEstimates, '[', FISCAL_YEAR_KEY))
    fiscalQuarters <- unlist (
        lapply (companyEstimates, '[', FISCAL_QUARTER_KEY))
    epss <- unlist (
        lapply (companyEstimates, '[', EPS_KEY))
    revenues <- unlist (
        lapply (companyEstimates, '[', REVENUE_KEY))
    userNames <- unlist (
        lapply (companyEstimates, '[', USERNAME_KEY))
    createdAts <- unlist (
        lapply (companyEstimates, '[', CREATED_AT_KEY))
    ids <- unlist (
        lapply (companyEstimates, '[', ID_KEY))
    analystIds <- unlist (
        lapply (companyEstimates, '[', ANALYST_ID_KEY))
    
    result <- data.frame (
        ticker,
        fiscalYears,
        fiscalQuarters,
        epss,
        revenues,
        userNames,
        createdAts,
        ids,
        analystIds)
    
    names (result) <- c(
        TICKER,
        FISCAL_YEAR,
        FISCAL_QUARTER,
        EPS,
        REVENUE,
        USERNAME,
        CREATED_AT,
        ID,
        ANALYST_ID)
    
    return (result)
}

#' Returns all estimates for a company specified by the ticker.
#'
#' @param ticker The Estimize ticker -- for example "MSFT".
#' @param year A four-digit year -- for example 1999.
#' @param quarter A numeric value between 1 and 4 inclusive. Note that if the
#'  quarter is used the year must also be set as well.
#'
#' @return The company estimates as a data frame.
#'
#' @examples
#'  \dontrun{
#'  result <- GetCompanyEstimates ("MSFT")
#'  result <- GetCompanyEstimates ("MSFT", 2015)
#'  result <- GetCompanyEstimates ("MSFT", 2015, 2)
#'  }
#'
#' @export
#'
GetCompanyEstimates <- function (
    ticker,
    year=NULL,
    quarter=NULL
) {

    .CheckParams (ticker, year, quarter)

    if (is.null (year) && is.null (quarter)) {
        destination <- paste (.baseURL, "companies/", ticker, "/estimates.json", sep="")
    } else if (!is.null (year) && is.null (quarter)) {
        destination <- paste (.baseURL, "companies/", ticker, "/estimates/", year, ".json", sep="")
    } else if (!is.null (year) && !is.null (quarter)) {
        destination <- paste (.baseURL, "companies/", ticker, "/estimates/", year, "/", quarter, ".json", sep="")
    } else {
        stop ("This function requires a year and a quarter (a quarter alone is not sufficient).")
    }

    result <- .DoGet (destination, andExecute = .TransformCompanyEstimates)

    return (result)
}

#' Returns all estimates in the specified date-range for all companies.
#'
#' @param startDate The start date, which must be in the format YYYY-MM-DD, for
#'  example "2015-01-20"
#' @param endDate The end date, which must be in the format YYYY-MM-DD, for
#'  example "2015-02-15"
#'
#' @return The estimates as a data frame.
#'
#' @examples
#'  \dontrun{
#'      result <- GetEstimates ("2015-01-20", "2015-02-15")
#'  }
#'
#' @export
#'
GetEstimates <- function (startDate, endDate) {

    .assertNotNull ("startDate", startDate)
    .assertNotNull ("endDate", endDate)

    # TODO:
    #.CheckDateFormat (startDate)
    #.CheckDateFormat (endDate)

    destination <- paste (.baseURL, "estimates.json?start_date=", startDate,
        "&end_date=", endDate, sep="")

    result <- .DoGet (destination, andExecute = .TransformCompanyEstimates)

    return (result)
}

#' Transforms a single consensus entry into a data frame. Note that if the entry
#' is null or empty this function will display an error and then return.
#'
#' @param dataSource For example, Wall Stret.
#'
#' @param type For example EPS.
#'
#' @param entries The entries being transformed.
#'
#' @return The release consensus entry as a data frame.
#'
.TransformConsensusEntry <- function (dataSource, type, entry) {

    if (is.null (entry)) {
        warning (paste ("The entry parameter is null for the dataSource: ",
                        dataSource, " and type: ", type, ".", sep=""))
        
        return (NULL)
    }

    mean <- entry[MEAN_KEY]
    high <- entry[HIGH_KEY]
    low <- entry[LOW_KEY]
    standardDeviation <- entry[STANDARD_DEVIATION_KEY]
    count <- entry[COUNT_KEY]
    updatedAt <- entry[UPDATED_AT_KEY]

    result <- data.frame (
        mean,
        high,
        low,
        standardDeviation,
        count,
        updatedAt)
    
    names (result) <- c(
        MEAN,
        HIGH,
        LOW,
        STANDARD_DEVIATION,
        COUNT,
        UPDATED_AT)

    result <- transform (result, Source = dataSource)
    result <- transform (result, Type = type)

    return (result)
}

#' Transforms consensus entries into a data frame. Note that if the entries are
#' null or empty this function will display an error and then return.
#'
#' @param dataSource For example, Wall Stret.
#'
#' @param type For example EPS.
#'
#' @param entries The entries being transformed.
#'
#' @return The release consensus entries as a data frame.
#'
.TransformConsensusEntries <- function (dataSource, type, entries) {

    if (is.null (entries) | length (entries) == 0) {
        warning (paste ("The entries parameter is null or of size zero for ",
            "the dataSource: ", dataSource, " and type: ", type, ".", sep=""))

        return (NULL)
    }

    mean <- unlist (
        lapply (entries, '[', MEAN_KEY))
    high <- unlist (
        lapply (entries, '[', HIGH_KEY))
    low <- unlist (
        lapply (entries, '[', LOW_KEY))
    standardDeviation <- unlist (
        lapply (entries, '[', STANDARD_DEVIATION_KEY))
    count <- unlist (
        lapply (entries, '[', COUNT_KEY))
    updatedAt <- unlist (
        lapply (entries, '[', UPDATED_AT_KEY))

    result <- data.frame (
        mean,
        high,
        low,
        standardDeviation,
        count,
        updatedAt)

    names (result) <- c(
        MEAN,
        HIGH,
        LOW,
        STANDARD_DEVIATION,
        COUNT,
        UPDATED_AT)

    result <- transform (result, Source = dataSource)
    result <- transform (result, Type = type)

    return (result)
}

#' Function transforms a releaseConsensus result from the call to
#' api.estimize.com into a data frame.
#'
#' @param releaseConsensus The result returned from the call to api.estimize.com
#'  specifically from the release consensus web service.
#'
#' @return The release consensus as a data frame.
#'
.TransformReleaseConsensus <- function (releaseConsensus) {

    wallStreetEPS <- releaseConsensus$wallstreet$eps
    wallStreetRevenue <- releaseConsensus$wallstreet$revenue
    wallStreetRevenueRevisions <- releaseConsensus$wallstreet$revenue$revisions
    wallStreetEPSRevisions <- releaseConsensus$wallstreet$eps$revisions

    estimizeEPS <- releaseConsensus$estimize$eps
    estimizeRevenue <- releaseConsensus$estimize$revenue
    estimizeRevisions <- releaseConsensus$estimize$revisions

    wallStreetEPSDataFrame <- .TransformConsensusEntry (
        dataSource = WALL_STREET,
        type = EPS,
        entry = wallStreetEPS)

    wallStreetEPSRevenueDataFrame <- .TransformConsensusEntry (
        dataSource = WALL_STREET,
        type = EPS_REVENUE,
        entry = wallStreetRevenue)

    wallStreetEPSRevisionsDataFrame <- .TransformConsensusEntries (
        dataSource = WALL_STREET,
        type = EPS_REVISION,
        entries = wallStreetEPSRevisions)

    wallStreetRevenueRevisionsDataFrame <- .TransformConsensusEntries (
        dataSource = WALL_STREET,
        type = REVENUE_REVISION,
        entries = wallStreetRevenueRevisions)

    estimizeEPSDataFrame <- .TransformConsensusEntries (
        dataSource = ESTIMIZE,
        type = EPS,
        entries = estimizeEPS)

    estimizeRevenueDataFrame <- .TransformConsensusEntries (
        dataSource = ESTIMIZE,
        type = REVENUE,
        entries = estimizeRevenue)

    result <- rbind (
        wallStreetEPSDataFrame,
        wallStreetEPSRevenueDataFrame,
        wallStreetEPSRevisionsDataFrame,
        wallStreetRevenueRevisionsDataFrame,
        estimizeEPSDataFrame,
        estimizeRevenueDataFrame
    )

    return (result)
}

#' Returns the current consensus as well as the consensus history of the
#' specified release. The user can obtain the id from results returned
#' from the GetCompanyReleases function.
#'
#' @param id The company identifier.
#'
#' @return The release consensus as a data frame.
#'
#' @examples
#'  \dontrun{
#'  result <- GetReleaseConsensus ("535c963053c804e0d50002a1")
#'  }
#'
#' @export
#'
GetReleaseConsensus <- function (id) {

    .assertNotNull ("id", id)

    destination <- paste (.baseURL, "releases/", id, "/consensus.json", sep="")

    result <- .DoGet (destination, andExecute = .TransformReleaseConsensus)

    return (result)
}

#' Function prints some information about this package.
#'
#' @examples
#'  \dontrun{
#'      About()
#'  }
#'
#' @export
#'
About <- function () {
  cat (
    " ***********************************************************\n",
    "***                                                     ***\n",
    "***      Welcome to the R Estimize.com API Package      ***\n",
    "***                                                     ***\n",
    "***                   version 0.8.5.                    ***\n",
    "***                                                     ***\n",
    "***                Follow us on LinkedIn:               ***\n",
    "***                                                     ***\n",
    "***       https://www.linkedin.com/company/229316       ***\n",
    "***                                                     ***\n",
    "***                Follow us on Twitter:                ***\n",
    "***                                                     ***\n",
    "***        https://twitter.com/CoherentLogicCo          ***\n",
    "***                                                     ***\n",
    "***********************************************************\n")
}