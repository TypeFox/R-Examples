##' Access to World Bank WDI API
##'
##' The function downloads data from the World Bank API.
##'
##'
##' @param name The new name to be used in the column.
##' @param indicator The World Bank official indicator name.
##' @param startDate The start date for the data to begin
##' @param endDate The end date.
##' @param printURL Whether the url link for the data should be printed
##' @param getMetaData Whether the data definition and the meta data
##' should be downloaded as well.
##' @param printMetaData logical, print out the meta data information
##' @param saveMetaData logical, whether meta data should be saved as a
##' local csv file
##' @param outputFormat The format of the data, can be 'long' or 'wide'.
##' @return A list containing the following elements
##' \describe{
##'     \item{data}{The country level data}
##'     \item{aggregates}{The aggregates provided by the World Bank}
##'     \item{metaData}{The metaData associated with the data}
##'     \item{results}{The status of the download, whether success/failed}
##' }
##' @export
##'
##' @seealso \code{\link{getWDI}}, \code{\link{getFAO}},
##' \code{\link{getFAOtoSYB}}
##' @examples
##' ## pop.df = getWDItoSYB(name = "total_population",
##' ##                      indicator = "SP.POP.TOTL")
##'
##'

getWDItoSYB = function(indicator = "SP.POP.0014.TO.ZS", name = NULL,
    startDate = 1960, endDate = format(Sys.Date(), "%Y"), printURL = FALSE,
    getMetaData = TRUE, printMetaData = FALSE, saveMetaData = FALSE,
    outputFormat = c("wide", "long")){

    outputFormat = match.arg(outputFormat)

    n = length(indicator)
    if(is.null(name))
        name = indicator

    if(length(name) != length(indicator))
        stop("length of name need to be the same as indicator")
    downloadInfo = unique(data.frame(name = name, indicator = indicator,
        stringsAsFactors = FALSE))
    wbData = data.frame(ISO2_WB_CODE = character(), Country = character(),
                        Year = integer(), stringsAsFactors = FALSE)
  ## wbData = data.frame(ISO2_WB_CODE = "US", Country = "United States",
  ##   Year = startDate, stringsAsFactors = FALSE)
    results = data.frame(Name = name, indicator = indicator,
                         Success = logical(length(name)),
                         Reason = character(length(name)),
                         Time = as.POSIXct(rep(NA, length(name))),
                         stringsAsFactors = FALSE)
    if(!(all(colSums(table(na.omit(downloadInfo))) == 1) &
         all(rowSums(table(na.omit(downloadInfo))) == 1)))
        stop("The relationship between the names and indicators are not 1-to-1")
    printLab(paste("World Bank Data Download (", n, " in Total)", sep = ""))
    i = 1
    retry = 1
    success = vector("logical", NROW(downloadInfo))
    while(i <= n){
        if(retry == 1)
            cat(paste("(", i, "): Downloading variable ", name[i], " ... ",
                      sep = ""))
        if(is.na(indicator[i])){
            cat("FAIL\n\tError: WDI indicator name missing\n")
            results[i, "Success"] = FALSE
            results[i, "Reason"] = "WDI indicator name missing"
            i = i + 1
            retry = 1
        } else {
            tmp = try(getWDI(indicator = downloadInfo[i, "indicator"],
                             name = name[i], startDate = startDate,
                             endDate = endDate, printURL = printURL,
                             outputFormat = outputFormat))
            if(!inherits(tmp, "try-error")){
                ## colnames(tmp) = c("Country", "ISO2_WB_CODE", "Year", name[i])
                if(outputFormat == "wide"){
                    wbData = merge(wbData, tmp, all = TRUE)
                } else if(outputFormat == "long"){
                    wbData = rbind(wbData, tmp)
                }
                results[i, "Success"] = TRUE
                results[i, "Reason"] = "Download Successful"
                results[i, "Time"] = Sys.time()
                i = i + 1
                retry = 1
                cat("OK\n")
            } else if(retry <= 10){
                retry = retry + 1
            } else if(retry > 10){
                cat("Download fail after 10 tries\n")
                results[i, "Success"] = FALSE
                results[i, "Reason"] = attr(tmp, "condition")$message
                i = i + 1
                retry = 1
            }
        }
    }

    if(getMetaData){
        metaData = getWDImetaData(indicator = na.omit(indicator),
                                  printMetaData = printMetaData,
                                  saveMetaData = saveMetaData)
    } else {
        metaData = NULL
    }

    entity.df = arrange(with(wbData, wbData[ISO2_WB_CODE %in%
        FAOcountryProfile[, "ISO2_WB_CODE"], ]), ISO2_WB_CODE, Year)
    region.df = arrange(with(wbData, wbData[!(ISO2_WB_CODE %in%
        FAOcountryProfile[, "ISO2_WB_CODE"]), ]), ISO2_WB_CODE, Year)

    cat(paste("\nNumber of variables successfully downloaded: ",
              sum(results$Success), " out of ", NROW(results), "\n", sep = ""))
    list(entity = entity.df, aggregates = region.df,
         metaData = metaData, results = results)
}

utils::globalVariables(names = c("ISO2_WB_CODE", "Year"))
