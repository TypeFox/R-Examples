##' Access to FAO FAOSTAT API.
##'
##' A function to access FAOSTAT data through the FAOSTAT API.
##'
##' Need to account for multiple itemCode, currently only support one
##' single variable.
##'
##' @param name The name to be given to the variable.
##' @param domainCode The domain of the data.
##' @param elementCode The code of the element.
##' @param itemCode The code of the specific item.
##' @param yearRange A numeric vector containing the years to be downloaded.
##' @param countrySet The FAOSTAT codes of those countries to be downloaded.
##' @param query The object created if using the FAOsearch function.
##' @param printURL Whether the url link for the data should be printed.
##' @param useCHMT logical, whether the CHMT function should be
##' applied to avoid double counting of China.
##' @param outputFormat The format of the data, can be 'long' or 'wide'.
##' @param returnNames Logical, should the area, the element and the
##' item names be reported?.
##' @param returnFlags, Logical, whether the flags should be
##' returned. Only work with outputFormat long.
##' @return Outputs a data frame containing the specified data.
##' @export
##' 
##'
##' @seealso \code{\link{getWDI}}, \code{\link{getWDItoSYB}},
##' \code{\link{getFAOtoSYB}}, \code{\link{FAOsearch}}
##'

getFAO = function(name = NULL, domainCode = "RL", elementCode = 5110,
    itemCode = 6621, query, printURL = FALSE,
    useCHMT = TRUE, outputFormat = "wide", returnNames = FALSE,
    returnFlags = FALSE, yearRange = NULL, countrySet = NULL){
    
    ## Year range
    if (!is.null(yearRange)) {
      if (!is.numeric(yearRange)) {
        stop("Please, provide a numeric vector for the year range.")
      } else {
        yearRange = paste(yearRange, collapse = ":")
      }
    }
    ## Country set
    if (!is.null(countrySet)) {
      if (!is.numeric(countrySet)) {
        stop("Please, provide a numeric vector for the year range.")
      } else {
        countrySet = paste(countrySet, collapse = ":")
      }
    }
    ## Query
    if(!missing(query)){
        if(NROW(query) > 1)
            stop("Use 'getFAOtoSYB' for batch download")
        domainCode = query$domainCode
        itemCode = query$itemCode
        elementCode = query$elementCode
        if(is.null(query$name)){
            name = with(query, paste(domainCode, itemCode, elementCode, sep = "_"))
        } else {
            name = query$name
        }
    }
    ## Name
    if(is.null(name))
        name = paste(domainCode, itemCode, elementCode, sep = "_")

    base = c("http://fenix.fao.org/wds/api?", "http://faostat3.fao.org/wds/api?", "http://fenixapps.fao.org/wds/api?")
    database = "db=faostat2&"
    selection = "select=A.AreaCode[FAOST_CODE],D.year[Year],D.value[Value]"
    from = "&from=data[D],element[E],item[I],area[A]&"
    condition = paste0("where=D.elementcode(", elementCode, "),D.itemcode(",
        itemCode, "),D.domaincode('", domainCode, "')")
    if (!is.null(yearRange)) {
      condition = paste0(condition, ",D.year(", yearRange, ")")
    }
    if (!is.null(countrySet)) {
      condition = paste0(condition, ",A.AreaCode(", countrySet, ")")
    }
    join = ",JOIN(D.elementcode:E.elementcode),JOIN(D.itemcode:I.itemcode),JOIN(D.areacode:A.areacode)&orderby=E.elementnamee,D.year"

    ## Flags
    if(returnFlags){
        outputFormat = "long"
        selection = paste0(selection, ",D.Flag[Flags]")
    }
    ## Names
    if(returnNames)
        selection = paste0(selection, "A.AreaNameE[AreaName],E.elementnamee[ElementName],I.itemnamee[ItemName]")
    ## Call
    out = "out=csv&"
    url = paste0(base, out, database, selection, from, condition, join)
    if(printURL)
        print(url)

    ## Allowing multiple server if any failed.
    for(i in 1:length(url)){
        faoData = suppressWarnings(try(read.csv(file = url[i],
            stringsAsFactors = FALSE), silent = TRUE))
        if(!inherits(faoData, "try-error"))
            break
    }
    faoData$FAOST_CODE = as.integer(faoData$FAOST_CODE)
    faoData$Year = as.integer(faoData$Year)
    ## CHMT
    if(useCHMT)
        faoData = CHMT(var = "Value", data = faoData, year = "Year")
    ## Output format
    if(outputFormat == "long" & NROW(faoData) != 0){
        faoData$domainCode = domainCode
        faoData$itemCode = itemCode
        faoData$elementCode = elementCode
        faoData$name = name
        faoData$Value <- as.numeric(gsub("n.a.", "", faoData$Value))
    } else if(outputFormat == "wide"){
        colnames(faoData)[colnames(faoData) == "Value"] = name
        faoData[, name] <- as.numeric(gsub("n.a.", "", faoData[, name]))
    }
    faoData
}
