setOldClass(c("quitte","data.frame"))


as.quitte <- function(x) {
  UseMethod("as.quitte",x)
}

as.quitte.quitte <- function(x) {
  if(is.quitte(x, warn=FALSE)) {
    return(x)
  } else {
    class(x) <- class(x)[class(x)!="quitte"]
    return(as.quitte.data.frame(x))
  }
}

setMethod("as.quitte", "quitte", as.quitte.quitte)

as.quitte.data.frame <- function (x) {
  store_attributes <- copy.attributes(x,0)
  mandatory_columns <- c("model","scenario","region","variable","unit","period","value")
  factor_columns <- c("model","scenario","region","variable","unit")
  colnames(x)[colnames(x)=="Cell"] <- "cell"
  colnames(x)[colnames(x)=="Region"] <- "region"
  colnames(x)[colnames(x)=="Year"] <- "period"
  colnames(x)[colnames(x)=="Value"] <- "value"
  colnames(x)[colnames(x)==paste0("Data",1)] <- "scenario"
  colnames(x)[colnames(x)==paste0("Data",2)] <- "model"
  colnames(x)[colnames(x)==paste0("Data",3)] <- "variable"
  
  if(!all(mandatory_columns %in% colnames(x))) {
    if(!("model"    %in% colnames(x))) x <- cbind(x,model=as.factor(NA))
    if(!("scenario" %in% colnames(x))) x <- cbind(x,scenario=as.factor(NA))
    if(!("region"   %in% colnames(x))) x <- cbind(x,region=as.factor("GLO"))
    if(!("variable" %in% colnames(x))) x <- cbind(x,variable=as.factor(NA))
    if(!("unit"     %in% colnames(x))) x <- cbind(x,unit=as.factor(NA))
    if(!("period"   %in% colnames(x))) x <- cbind(x,period=as.POSIXct(NA))
    if(!("value"    %in% colnames(x))) stop("Data frame cannot be converted. A column \"value\" has to be provided!")      
  }
  factor_check <- sapply(x[,factor_columns],is.factor)
  if(!all(factor_check)) {
    for(i in names(factor_check)[!factor_check]) x[[i]] <- as.factor(x[[i]])
  }
  ISOyear <- make.ISOyear()
  if(!("POSIXct" %in% attr(x$period,"class"))) x$period <- ISOyear(x$period)
  if(!is.numeric(x$value)) stop("Value column must contain numeric data!")
  
  #rearrange data for better readability
  reorder <- c(mandatory_columns[mandatory_columns!="value"], names(x)[!(names(x) %in% mandatory_columns)], "value")
  x <- x[reorder]
  
  class(x) <- c("quitte","data.frame")
  return(copy.attributes(store_attributes,x))
}

setMethod("as.quitte", "data.frame", as.quitte.data.frame)

as.quitte.magpie <- function (x) {
  x <- clean_magpie(x,what="sets")
  store_attributes <- copy.attributes(x,0)
  d <- dimnames(x)
  if(!is.null(names(d)[[3]])) {
    datanames <- strsplit(names(d)[[3]],"\\.")[[1]]
    datanames <- make.unique(c("cell","region","year","value",datanames),sep="")[-(1:4)]
  } else {
    datanames <- NULL
  }
  x <- as.data.frame(x)
  if(all(is.na(x$Cell))) x$Cell <- NULL
  colnames(x)[colnames(x)=="Cell"] <- "cell"
  colnames(x)[colnames(x)=="Region"] <- "region"
  colnames(x)[colnames(x)=="Year"] <- "period"
  colnames(x)[colnames(x)=="Value"] <- "value"
  if(length(datanames)>0) {
    for(i in 1:length(datanames)) colnames(x)[colnames(x)==paste0("Data",i)] <- datanames[i]
  } else {
    if("Data1" %in% colnames(x)) if(all(levels(x$Data1)=="NA")) x$Data1 <- NULL
  }
  if(all(x$period==0)) {
    levels(x$period) <- NA  
  } else {
    ISOyear <- make.ISOyear()
    x$period <- ISOyear(x$period)
  }
  return(copy.attributes(store_attributes,as.quitte.data.frame(x)))
}

setMethod("as.quitte", "magpie", as.quitte.magpie)
setMethod("as.quitte", "quitte", as.quitte.quitte)
setMethod("as.quitte", "data.frame", as.quitte.data.frame)