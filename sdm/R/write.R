# Author: Babak Naimi, naimi.b@gmail.com
# Date :  April 2016
# Version 1.1
# Licence GPL v3
#--------

if (!isGeneric("read.sdm")) {
  setGeneric("read.sdm", function(filename,...)
    standardGeneric("read.sdm"))
}


setMethod('read.sdm', signature(filename='character'),
          function(filename,...) {
            o <- readRDS(file=filename,...)
            if (!class(o) %in% c('sdmdata','sdmModels','.sdmCorSetting')) stop('the file is not a sdm object!')
            o
          }
)

if (!isGeneric("write.sdm")) {
  setGeneric("write.sdm", function(x,filename,overwrite,...)
    standardGeneric("write.sdm"))
}


setMethod('write.sdm', signature(x='sdmModels',filename='character'),
          function(x,filename,overwrite,...) {
            if (missing(overwrite)) overwrite <- FALSE
            filename <- .trim(filename)
            if (extension(filename) == '') filename <- paste(filename,'.sdm',sep='')
            else if (extension(filename) != '.sdm') filename <- paste(filename,'.sdm',sep='')
            
            if (!overwrite && file.exists(filename)) stop('a file with the same name exists; choose different name or use overwrite=TRUE')
            saveRDS(x,file=filename,...)
          }
)
setMethod('write.sdm', signature(x='sdmdata',filename='character'),
          function(x,filename,overwrite,...) {
            if (missing(overwrite)) overwrite <- FALSE
            filename <- .trim(filename)
            if (extension(filename) == '') filename <- paste(filename,'.sdd',sep='')
            else if (extension(filename) != '.sds') filename <- paste(filename,'.sds',sep='')
            if (!overwrite && file.exists(filename)) stop('a file with the same name exists; choose different name or use overwrite=TRUE')
            saveRDS(x,file=filename,...)
          }
)

setMethod('write.sdm', signature(x='.sdmCorSetting',filename='character'),
          function(x,filename,overwrite,...) {
            if (missing(overwrite)) overwrite <- FALSE
            filename <- .trim(filename)
            if (extension(filename) == '') filename <- paste(filename,'.sds',sep='')
            else if (extension(filename) != '.sds') filename <- paste(filename,'.sds',sep='')
            if (!overwrite && file.exists(filename)) stop('a file with the same name exists; choose different name or use overwrite=TRUE')
            saveRDS(x,file=filename,...)
          }
)