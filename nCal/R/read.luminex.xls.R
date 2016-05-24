## Author: Krisztian Sebestyen

## download and install perl interpreter from http://www.activestate.com/activeperl/downloads
## more options here http://www.perl.org/get.html#win32

## assumes a specific format for xls file
## luminex header block - stored in attr( . ,'info')
## luminex data block - returned as data frame, assumed to have a column-name-header of its own
## luminex tail block - stored in attr( . ,'summary')

## read all xls sheets
read.luminex.xls<-function(file,verbose = FALSE,sheets = NULL, assay_id=NULL, na.strings = c("NA", "#DIV/0!"), ..., perl = "perl") 
{
    #require(gdata)
    if(is.null(sheets)) sheets <- try(sheetNames(file,perl = perl))
    if(inherits(sheets,'try-error')) sheets <- 1

    if (verbose) print("read.luminex.xls 200")
    
    docs <- NULL
    for(i in 1:length(sheets))
        docs[[i]] <- read.luminex.sheet.xls(file, sheet = sheets[i], verbose = verbose, na.strings = na.strings, ..., perl = perl) 
    names(docs) <- sheets
    
    if(length(docs) == 1)
        docs <- docs[[1]]
    else 
        docs <- do.call(rbind, docs)

    colnames(docs) <- tolower(colnames(docs))
    colnames(docs)[colnames(docs)=="description"]="sample_id"
    colnames(docs)[colnames(docs)=="exp.conc"]="expected_conc"
    colnames(docs)[colnames(docs)=="type"]="well_role"  
    
    docs$well_role[startsWith(docs$well_role, "X")]="Unknown"
    docs$well_role[startsWith(docs$well_role, "S")]="Standard"
    docs$well_role[startsWith(docs$well_role, "B")]="Background"
    docs$well_role[startsWith(docs$well_role, "C")]="QC"
    
    if (is.null(assay_id)) docs$assay_id=round(runif(1)*10000)
    
    docs
}


read.luminex.fixed.sheet.xls<-function (xls, sheet = 1, lengths = list(header = 7,data = 1),verbose = FALSE, na.strings = c("NA", "#DIV/0!"), ..., perl = "perl") 
{
    con <- tfn <- NULL
    on.exit({
        err <- FALSE
        if (inherits(con, "connection")) {
            tryCatch(op <- isOpen(con), error = function(x) err <<- TRUE)
            if (!err && op) close(con)
        }
        if (file.exists(tfn)) file.remove(tfn)
    })
    xls <- path.expand(xls)
    perl <- if (missing(perl)) 
        findPerl(verbose = verbose)
    else findPerl(perl, verbose = verbose)
    
    # if wanted to automatically find the data block and tail block then 
    # should not remove blank lines and should find blank lines - not as easy to do it efficiently in R
    #   con <- xls2sep(xls, sheet, verbose = verbose, ..., method = method,perl = perl,blank.lines.skip = FALSE)
    
    con <- xls2sep(xls, sheet, verbose = verbose, as.is=TRUE, ..., method = 'csv',perl = perl)
    tfn <- summary(con)$description
    close(con)
    
    luminex.header <- readLines(tfn,n = lengths$header)
    luminex.data <- read.csv(tfn, na.strings = na.strings,skip = lengths$header,nrows = lengths$data, as.is=TRUE, ...)
    luminex.tail <- scan(tfn,skip = sum(unlist(lengths)),what = 'character')
    attr(luminex.data,'info') <- luminex.header
    attr(luminex.data,'summary') <- luminex.tail
    luminex.data
}

## assumes a specific format for xls file
## luminex header block - stored in attr( . ,'info')
## luminex data block - returned as data frame, assumed to have a column-name-header of its own
## luminex tail block - stored in attr( . ,'summary')
read.luminex.sheet.xls<-function (xls, sheet = 1, lengths = NULL, verbose = FALSE, na.strings = c("NA", "#DIV/0!"), ..., perl = "perl") 
{
    con <- tfn <- NULL
    on.exit({
        err <- FALSE
        if (inherits(con, "connection")) {
            tryCatch(op <- isOpen(con), error = function(x) err <<- TRUE)
            if (!err && op) close(con)
        }
        if (file.exists(tfn)) file.remove(tfn)
    })
    xls <- path.expand(xls)
    perl <- if (missing(perl)) 
        findPerl(verbose = verbose)
    else findPerl(perl, verbose = verbose)
    
    # if wanted to automatically find the data block and tail block then 
    # should not remove blank lines and should find blank lines - not as easy to do it efficiently in R
    
    if (verbose) print("read.luminex.sheet.xls 200")
    
    con <- xls2sep(xls, sheet, verbose = verbose, blank.lines.skip = FALSE, method = 'csv',perl = perl)
    tfn <- summary(con)$description
    close(con)
    
    # find blank lines that separate luminex header, luminex data and luminex tail blocks
    # a blank line in a .csv file consists of commas only
    if(!length(lengths)){
        x <- readLines(tfn)
        lengths <- grep("^,+$",x)
    }
    luminex.header <- readLines(tfn,n = lengths[1])
#   luminex.header <- x[1:(lengths[1]-1)]
    ndata <- lengths[2] - lengths[1] - 1 - 1
    luminex.data <- read.csv(tfn, na.strings = na.strings,skip = lengths[1],nrows = ndata,as.is=TRUE,...)
    luminex.tail <- scan(tfn,skip = sum(unlist(lengths)),what = 'character', quiet=TRUE)
    #attr(luminex.data,'info') <- luminex.header
    #attr(luminex.data,'summary') <- luminex.tail
    names(luminex.data)=tolower(names(luminex.data))
    luminex.data
}

#lumi2rumi <- function(data){
#        
#    cn.data[which(cn.data=='analyte')] <- 'bead_type'
#    cn.data[which(cn.data=='well')] <- 'replicate'
#    cn.data[which(cn.data=='type')] <- 'well_role'
#    cn.data[which(cn.data=='obs.conc')] <- 'starting_conc'
#    cn.data[which(cn.data=='exp.conc')] <- 'expected_conc'
#    cn.data[which(cn.data=='description')] <- 'id'
#    
#    all.valid.rumi.colums <- c('bead_type','fi','well_role','dilution',
#    'starting_conc','expected_conc','assay_id','sample_id','replicate')
#    
#    
#    cx <- match(all.valid.rumi.colums,cn.data,nomatch = NA) 
#    
#    # if(any(is.na(cx)){ 
#        # tx <- "Some ruminex columns could not be matched in data :"
#        # txt <- paste(txt,setdiff(all.valid.rumi.colums,cn.data))
#        # warning(txt)
#    # }
#    data <- data[,cx]
#    
#    data
#}


## findPerl is copied (since not exported) from gdata
findPerl <- function (perl, verbose = "FALSE") 
{
    errorMsg <- "perl executable not found. Use perl= argument to specify the correct path."
    if (missing(perl)) {
        perl = "perl"
    }
    perl = Sys.which(perl)
    if (perl == "" || perl == "perl") 
        stop(errorMsg)
    if (.Platform$OS == "windows") {
        if (length(grep("rtools", tolower(perl))) > 0) {
            perl.ftype <- shell("ftype perl", intern = TRUE)
            if (length(grep("^perl=", perl.ftype)) > 0) {
                perl <- sub("^perl=\"([^\"]*)\".*", "\\1", perl.ftype)
            }
        }
    }
    if (verbose) 
        cat("Using perl at", perl, "\n")
    perl
}
