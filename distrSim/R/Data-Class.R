#########
## about Dataclass: slot Data
##
### changed from version 1.8 on:
## ith observation in ith line of datamatrix/array
## jth item/dimension of each observation in jth column of datamatrix/array
## kth run/time of each observation in kth slide of datamatrix/array

## ++old
## +ith run in ith line of datamatrix
## +jth samples of each run in jth column of datamatrix


################################
##
## Class: Dataclass
##
################################

## Access methods

setMethod("filename", "Dataclass", function(object) object@filename)
setMethod("Data", "Dataclass", function(object) object@Data)
setMethod("runs", "Dataclass", function(object) object@runs)
setMethod("samplesize", "Dataclass", function(object) object@samplesize)
setMethod("obsDim", "Dataclass", function(object) object@obsDim)### new v.1.8
setMethod("name", "Dataclass", function(object) object@name)### new v.1.8
setMethod("getVersion", "Dataclass", 
           function(object) object@version)### new v.1.8

setMethod("Dataclass","DataframeorSeqDataFrames", function(Data, filename = NULL, name = "Data-Set"){ 
 if(!is(Data,"SeqDataFrames"))
    Data <- SeqDataFrames(Data)
 runs0 <-  runs(Data) 
 obsDim0 <- obsDim(Data)
 samplesize0 <- samplesize(Data)
 new("Dataclass", filename = filename, Data = Data, runs = runs0, 
     obsDim = obsDim0, samplesize = samplesize0, name = name)
})

setMethod("Dataclass","ArrayorNULLorVector", function(Data, filename = NULL, name = "Data-Set") 
{if(is.null(Data))
    stop("generating an object of class \"Dataclasss\" requires data of type \"array\" or \"vector\"")
 
 runs0 <- 1
 obsDim0 <- 1
 samplesize0 <- length(Data)
 dimnames0 <- NULL
 rnames <- NULL
 dnames <- NULL
 snames <- NULL
 if(!is.null(names(Data)))  
    {snames <- names(Data)
     dimnames0 <- list(snames, dnames, rnames)
     }
 Data0 <- array(data = Data, dim = c(samplesize0,obsDim0,runs0),
                dimnames = dimnames0)       
 new("Dataclass", filename = filename, Data = Data0, runs = runs0, 
     obsDim = obsDim0, samplesize = samplesize0, name = name)
})

setMethod("Dataclass","matrix", function(Data, filename = NULL, name = "Data-Set") 
{if(is.null(Data))
    stop("generating an object of class \"Dataclasss\" requires data of type \"array\" or \"vector\"")
 
 runs0 <- 1
 rnames <- NULL
 dimnames0 <- NULL
 obsDim0 <- ncol(Data)   
 samplesize0 <- nrow(Data)
 if(!is.null(dimnames(Data))) 
    {dnames <- colnames(Data)
     snames <- rownames(Data)
     dimnames0 <- list(snames, dnames, rnames)
 }    
 Data0 <- array(data = Data, dim = c(samplesize0,obsDim0,runs0),
                dimnames = dimnames0)       
 new("Dataclass", filename = filename, Data = Data0, runs = runs0, 
     obsDim = obsDim0, samplesize = samplesize0, name = name)      
})

setMethod("Dataclass","array", function(Data, filename = NULL, name = "Data-Set") 
{if(is.null(Data))
    stop("generating an object of class \"Dataclasss\" requires data of type \"array\" or \"vector\"")
 if(length(dim(Data))>3) 
    stop("generating an object of class \"Dataclasss\" requires data of type \"array\" in atmost 3 dimensions")

 samplesize0 <- dim(Data)[1] 
 snames  <- dimnames(Data)[1]

 if(length(dim(Data))==1) {
    obsDim0 <- 1 
    dnames  <- NULL
 }else{
    obsDim0 <- dim(Data)[2] 
    dnames  <- dimnames(Data)[2]
 }
 if(length(dim(Data))<=2) {
    runs0  <- 1 
    rnames <- NULL
 }else{
    runs0  <- dim(Data)[3] 
    rnames <- dimnames(Data)[3]
 }
 dimnames0 <- list(snames, dnames, rnames)
 Data0 <- array(data = Data, dim = c(samplesize0,obsDim0,runs0),
                dimnames = dimnames0)       
 new("Dataclass", filename = filename, Data = Data0, runs = runs0, 
     obsDim = obsDim0, samplesize = samplesize0, name = name)      
})

 


## Replacement methods

setReplaceMethod("name", "Dataclass", 
                  function(object, value)
                          { object@name <- value; object}) ### new 1.8
setReplaceMethod("filename", "Dataclass", 
                  function(object, value){ object@filename <- value; object})

setReplaceMethod("Data", "Dataclass",
                 function(object, value){
                   datdim <- dim(value)
                   object <- new("Dataclass",
                                 filename = filename(object),
                                 Data = value, samplesize=datdim[1],
                                 obsDim = datdim[2], 
                                 runs = datdim[3])
                   object
                 })




