SeqDataFrames <- function(...) new("SeqDataFrames", data=list(...))

## Access methods

setMethod("samplesize", "SeqDataFrames", function(object) as.numeric(lapply(object@data,function(y) nrow(y))))
setMethod("obsDim", "SeqDataFrames", function(object) ncol(object@data[[1]]))
setMethod("runs", "SeqDataFrames", function(object) length(object@data))
setMethod("obsdimnames", "SeqDataFrames", function(object) colnames(object@data[[1]]))
setMethod("names", "SeqDataFrames", function(x) names(x@data))
setMethod("runnames", "SeqDataFrames", function(object) names(object@data))

## Replacement methods

setReplaceMethod("obsdimnames", "SeqDataFrames", function(object,value) {colnames(object@data[[1]])<- value; object})
setReplaceMethod("names", "SeqDataFrames", function(x,value) {names(x@data)<- value; x})
setReplaceMethod("runnames", "SeqDataFrames", function(object,value) {names(object@data)<- value; object})


