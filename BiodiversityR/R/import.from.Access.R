if(.Platform$OS.type == "windows") {

`import.from.Access` <-
function(file=file.choose(), data.type="community", table=NULL, sitenames="sites", column="species", value="abundance", factor="", level="", cepnames=FALSE) {
#    if (!require(RODBC)) {stop("Requires package RODBC")}
    dataplace <- RODBC::odbcConnectAccess(file)
    if (is.null(data.type) == TRUE) {data.type <- table}
    TYPES <- c("community", "environmental", "stacked")
    data.type <- match.arg(data.type, TYPES)
    if (is.null(table) == TRUE) {table <- data.type}
    if (data.type == "stacked") {
        stackeddata <- RODBC::sqlFetch(dataplace, table)  
        result <- makecommunitydataset(stackeddata, row=sitenames, column=column, value=value, factor=factor, level=level)
    }else{
         result <- RODBC::sqlFetch(dataplace, table, rownames=sitenames)
    }
    close(dataplace)
    rownames(result) <- make.names(rownames(result), unique=T)
    if (cepnames == TRUE) {
        colnames(result) <- make.cepnames(colnames(result))
    }else{
        colnames(result) <- make.names(colnames(result), unique=T)
    }
    return(result)
}

`import.from.Access2007` <-
function(file=file.choose(), data.type="community", table=NULL, sitenames="sites", column="species", value="abundance", factor="", level="", cepnames=FALSE) {
#    if (!require(RODBC)) {stop("Requires package RODBC")}
    dataplace <- RODBC::odbcConnectAccess2007(file)
    if (is.null(data.type) == TRUE) {data.type <- table}
    TYPES <- c("community", "environmental", "stacked")
    data.type <- match.arg(data.type, TYPES)
    if (is.null(table) == TRUE) {table <- data.type}
    if (data.type == "stacked") {
        stackeddata <- RODBC::sqlFetch(dataplace, table)  
        result <- makecommunitydataset(stackeddata, row=sitenames, column=column, value=value, factor=factor, level=level)
    }else{
         result <- RODBC::sqlFetch(dataplace, table, rownames=sitenames)
    }
    close(dataplace)
    rownames(result) <- make.names(rownames(result), unique=T)
    if (cepnames == TRUE) {
        colnames(result) <- make.cepnames(colnames(result))
    }else{
        colnames(result) <- make.names(colnames(result), unique=T)
    }
    return(result)
}

}