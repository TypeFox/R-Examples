if(.Platform$OS.type == "windows") {

`import.from.Excel` <-
function(file=file.choose(), data.type="community", sheet=NULL, sitenames="sites", 
    column="species", value="abundance", factor="", level="", cepnames=FALSE, 
    write.csv=FALSE, csv.file=paste(data.type, ".csv", sep="") )
{
#    if (!require(RODBC)) {stop("Requires package RODBC")}
    dataplace <- RODBC::odbcConnectExcel(file)
    if (is.null(data.type) == TRUE) {data.type <- sheet}
    TYPES <- c("community", "environmental", "stacked")
    data.type <- match.arg(data.type, TYPES)
    if (is.null(sheet) == TRUE) {sheet <- data.type}
    if (data.type == "stacked") {
        stackeddata <- RODBC::sqlFetch(dataplace,sheet)  
        result <- makecommunitydataset(stackeddata,row=sitenames,column=column,value=value,factor=factor,level=level)
        data.type <- "community"
    }else{
        result <- RODBC::sqlFetch(dataplace,sheet,rownames=sitenames)
    }
    close(dataplace)
    rownames(result) <- make.names(rownames(result),unique=T)
    if (cepnames == TRUE && data.type == "community") {
        colnames(result) <- make.cepnames(colnames(result))
    }else{
        colnames(result) <- make.names(colnames(result),unique=T)
    }
    if (write.csv == TRUE) {utils::write.table(x=result, file=csv.file, row.names=T, col.names=T, sep=',')}
    return(result)
}

`import.from.Excel2007` <-
function(file=file.choose(), data.type="community", sheet=NULL, sitenames="sites", 
    column="species", value="abundance", factor="", level="", cepnames=FALSE, 
    write.csv=FALSE, csv.file=paste(data.type, ".csv", sep="") )
{
#    if (!require(RODBC)) {stop("Requires package RODBC")}
    dataplace <- RODBC::odbcConnectExcel2007(file)
    if (is.null(data.type) == TRUE) {data.type <- sheet}
    TYPES <- c("community", "environmental", "stacked")
    data.type <- match.arg(data.type, TYPES)
    if (is.null(sheet) == TRUE) {sheet <- data.type}
    if (data.type == "stacked") {
        stackeddata <- RODBC::sqlFetch(dataplace,sheet)  
        result <- makecommunitydataset(stackeddata,row=sitenames,column=column,value=value,factor=factor,level=level)
    }else{
         result <- RODBC::sqlFetch(dataplace,sheet,rownames=sitenames)
    }
    close(dataplace)
    rownames(result) <- make.names(rownames(result),unique=T)
    if (cepnames == TRUE && data.type == "community") {
        colnames(result) <- make.cepnames(colnames(result))
    }else{
        colnames(result) <- make.names(colnames(result),unique=T)
    }
    if (write.csv == TRUE) {utils::write.table(x=result, file=csv.file, row.names=T, col.names=T, sep=',')}
    return(result)
}

}