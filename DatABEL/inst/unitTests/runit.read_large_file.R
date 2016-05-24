### --- Test setup ---

if(FALSE) {
    ## Not really needed, but can be handy when writing tests
    library("RUnit")
    library("DatABEL")
}
#test.empty <- function(){}
### do not run
#stop("SKIP THIS TEST")
###

# common functions and data

quiet <- TRUE

data_prefix = "big_data/ERF.merlin.22.collected.ped.out.mldose"

databel_large_file_reads <- function(fv_data)
{
    # collect check data
    varsum_filename <- paste(data_prefix,"_varsum",sep="")
    varsum <- scan(varsum_filename)
    obssum_filename <- paste(data_prefix,"_obssum",sep="")
    obssum <- scan(obssum_filename)
    varnames_filename <- paste(data_prefix,"_varnames",sep="")
    varnames <- scan(varnames_filename,what="character")
    obsnames_filename <- paste(data_prefix,"_obsnames",sep="")
    obsnames <- scan(obsnames_filename,what="character")
    
    
    # check properties are there
    realdim <- c(length(obsnames),length(varnames))
    checkEquals(realdim,dim(fv_data))
    checkEquals(realdim[1]*realdim[2],length(fv_data))
    checkTrue(is.list(dimnames(fv_data)))
    checkEquals(data_prefix,backingfilename(fv_data))
    checkTrue(is.list(get_dimnames(fv_data)))
    checkEquals(2,length(get_dimnames(fv_data)))
    checkTrue(is.integer(cachesizeMb(fv_data)))
    
    # test that names are ok    
    checkEquals(varnames,get_dimnames(fv_data)[[2]])
    checkEquals(obsnames,get_dimnames(fv_data)[[1]])
    checkEquals(varnames,dimnames(fv_data)[[2]])
    checkEquals(obsnames,dimnames(fv_data)[[1]])
    
    # test read all    
    accumulated_obssum <- rep(0,dim(fv_data)[1])
    
    gsize <- 1000
    for (i in seq(1,dim(fv_data)[2],gsize)) {
        eow <- min(dim(fv_data)[2],i+gsize-1)
        dta <- as(fv_data[,i:eow],"matrix")
        dimnames(dta) <- NULL
        checkEquals(varsum[i:eow],apply(dta,FUN=sum,MAR=2),tolerance=5*sqrt(.Machine$double.eps))
        cat("checked",eow,"out of",dim(fv_data)[2],"variables\n")
        accumulated_obssum <- accumulated_obssum + apply(dta,FUN=sum,MAR=1)
    }
    checkEquals(obssum,as(accumulated_obssum,"vector"),tolerance=5*sqrt(.Machine$double.eps))
    
    # as ...
    nCol <- round(runif(1,min=10,max=30))
    nRow <- round(runif(1,min=10,max=30))
    rancols <- sample(1:dim(fv_data)[2],nCol)
    ranrows <- sample(1:dim(fv_data)[1],nRow)
    testmatr <- as(fv_data[ranrows,rancols],"matrix")
    checkTrue(is.matrix(testmatr))
    checkTrue(is.numeric(testmatr))
    checkTrue(is.double(testmatr))
    checkEquals(c(nRow,nCol),dim(testmatr))
    
    
}

### --- Test functions ---

test.databel_large_file_reads <- function()
{
    
    #library("RUnit")
    #library("DatABEL")
    
    ## try to read from non-existing file
    checkException(databel("no_such_file"))
    checkException(databel("no_such_file1",128))
    checkException(databel(c("no_such_file2","aaa")))
    checkException(databel(127))
    checkException(databel(TRUE))
    
    # cann I read from base file name, *.fvi, *fvd?
    tmpname1 <- paste(data_prefix,".fvi",sep="")
	if (!quiet) print(tmpname1)
#    if (file.exists(tmpname1)) {
#        fv_data <- databel(tmpname1)
#    }
    tmpname2 <- paste(data_prefix,".fvd",sep="")
	if (!quiet) print(tmpname2)
#    if (file.exists(tmpname2)) {
#        fv_data <- databel(tmpname2)
#    }
    
    if (!(file.exists(tmpname1) && file.exists(tmpname1))) {
        warning(paste("file",tmpname1,"is missing; stopping tests"))
        return(NULL)
    }
    fv_data <- databel(data_prefix)
    
    # check if nothing works with disconnected
    disconnect(fv_data)
    checkTrue(is.null(disconnect(fv_data)))
    #checkException(fv_data[,1])
    #checkException(fv_data[1,])
    #checkException(dim(fv_data))
    #checkException(get_dimnames(fv_data))
    #checkException(cachesizeMb(fv_data) <- 12)
    
    # check that some properties are still there
    checkTrue(is.integer(cachesizeMb(fv_data)) && cachesizeMb(fv_data) > 0)
    checkEquals(data_prefix,backingfilename(fv_data))
    
    connect(fv_data)
    checkTrue(is.null(connect(fv_data)))
    
    # test reading with different chache size
    databel_large_file_reads(fv_data)
	if (!quiet) print(fv_data)
    checkException(cachesizeMb(fv_data) <- "aaa")
    checkException(cachesizeMb(fv_data) <- -1)
    cachesizeMb(fv_data) <- 0
    checkTrue(cachesizeMb(fv_data)>0)
    databel_large_file_reads(fv_data)
    cachesizeMb(fv_data) <- 4
	if (!quiet) print(cachesizeMb(fv_data))
    checkEquals(4,cachesizeMb(fv_data))
    cachesizeMb(fv_data) <- 180
    checkEquals(180,cachesizeMb(fv_data))
    cachesizeMb(fv_data) <- 96
    checkEquals(96,cachesizeMb(fv_data))
    databel_large_file_reads(fv_data)
    
}
