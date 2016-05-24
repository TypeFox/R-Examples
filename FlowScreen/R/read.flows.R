#' Read .csv or .Rdata file of streamflows 
#' 
#' Reads .csv or .Rdata files of daily streamflow time series.  Recognizes 
#' several formats, including those used by Water Survey Canada and the United 
#' States Geological Survey. Uses read.csv() or load() functions from base package, 
#' returns data frame with ID, Date, and Flow, and, if available,
#' associated quality codes and source agency. Replaces negative values that are 
#' sometimes used to denote missing data with NAs.
#' @param filename name of .csv file to be read from.
#' @author Jennifer Dierauer
#' @export

read.flows <- function(filename) {

    my.split <- unlist(strsplit(filename, ".", fixed=T))
    suff <- tolower(my.split[length(my.split)])
    
    if (length(grep(suff, c("rdata", "csv"))) == 0) {
        stop("File must be a .Rdata or .csv file")
    }
    
    if (suff == "csv") {
        ddata <- utils::read.csv(filename)
    }
    
    if (suff == "rdata") {
        the.object <- load(filename)
        ddata <- get(the.object)
    }
    
    
    ## change ID from factor to character
    mID <- grep("id", tolower(colnames(ddata)))
    ddata[,mID] <- as.character(ddata[,mID])
    
    mdatec <- grep("date", tolower(colnames(ddata)))  ## look for "date" column
    mdates <- as.character(ddata[,mdatec])
    
    
    ## identify file as WSC based on PARAM column, remove extra end rows, and format dates
    parcol <- grep("^PARAM$", colnames(ddata))
    
    if (length(parcol) != 0) {
        
        cut <- length(ddata[,1])-2
        ddata <- ddata[1:cut,]
        mdates <- mdates[1:cut]
        
        if (substr(mdates[1],3,3)=="/") {
            
            d.split <- unlist(strsplit(mdates, "/", fixed=T))
            mseq <- seq(from=1, to=length(d.split), by=3)
            
            if (max(as.numeric(d.split[mseq])) > 12) {
                
                mdates <- as.Date(mdates, format="%d/%m/%Y")
                mdates <- as.Date(mdates, format(mdates, "%Y/%m/%d"))
                
            } else {
            
                mdates <- as.Date(mdates, format="%m/%d/%Y")
                mdates <- as.Date(mdates, format(mdates, "%Y/%m/%d"))
            
            }
            
        } else {
            
            if (substr(mdates[1],5,5)=="/") {
                mdates <- as.Date(mdates, format="%Y/%m/%d")
            }
        }
        
        ## format data frame for return
        mFlow <- grep("flow", tolower(colnames(ddata)))
        msym <- grep("sym", tolower(colnames(ddata)))
        
        ddata <- data.frame(ID=ddata[,mID], Date=mdates, Flow=ddata[,mFlow], 
                            SYM=ddata[,msym], Agency="WSC")
        
    }
    
    ## identify files with date column, but no PARAM column
    ##  %m/%d/%Y date format 
    ## or %Y-%m-%d
    
    if (length(mdatec) != 0 && length(parcol) == 0) {
        
        
        ## check formatting of date
        if (substr(mdates[1],3,3)=="/") {
            
            mdates <- as.Date(mdates, format="%m/%d/%Y")
            mdates <- as.Date(mdates, format(mdates, "%Y/%m/%d"))
        }
        
        if (substr(mdates[1],5,5)=="-") {
            mdates <- as.Date(mdates, format="%Y-%m-%d")
            mdates <- as.Date(mdates, format(mdates, "%Y/%m/%d"))
        }
                    
        ## get agency column if it exists
        magc <- grep("agency", tolower(colnames(ddata)))
        if (length(magc) == 0) {
            magency <- rep(NA, length(ddata[,1]))
        } else {
            magency <- as.character(ddata[,magc])
        }
        
        mFlow <- grep("flow", tolower(colnames(ddata)))  ## get flow column
        
        if (length(mFlow) == 0) {
            ## if no flow column, look for val column
            mFlow <- grep("val", tolower(colnames(ddata))) 
        }
        
        ## find quality codes in a sym, flag, or code column, if exist
        ## fill column with NA if no matching column is found
        msymc <- grep("sym", tolower(colnames(ddata))) ##look for sym column
        
        if (length(msymc) == 0) {
            msymc <- grep("code", tolower(colnames(ddata))) ## if none, look for code column
        }
        
        if (length(msymc) == 0) {
            msymc <- grep("flag", tolower(colnames(ddata))) ## if still none, look for flag column
        }
        
        if (length(msymc) == 0) {
            msym <- rep(NA, length(ddata$Flow))  ## if not match, fill column with NA
        } else {
            msym <- ddata[,msymc]  ## if match, get values from correct column
        }
        
        ## format data frame for return
        ddata <- data.frame(ID=ddata[,mID], Date=mdates, Flow=ddata[,mFlow], SYM=msym, 
                            Agency=magency)
            
    }
    
    ## identify files with Day, Month, Year columns and format to match Water Survey Canada files
    if (length(grep("^day$", tolower(colnames(ddata)))) != 0) {
        mday <- grep("^day$", tolower(colnames(ddata)))
        mmonth <- grep("^month$", tolower(colnames(ddata)))
        myr <- grep("^year$", tolower(colnames(ddata)))
        
        mdates <- as.Date(paste(ddata[,myr], ddata[,mmonth], ddata[,mday], sep="/"),
                          format="%Y/%m/%d")
        mFlow <- grep("flow", tolower(colnames(ddata)))
        
        ## get agency column if it exists
        magc <- grep("agency", tolower(colnames(ddata)))
        if (length(magc) == 0) {
            magency <- rep(NA, length(ddata$Flow))
        } else {
            magency <- as.character(ddata[,magc])
        }
        
        ## find quality codes in a sym, flag, or code column, if exist
        ## fill column with NA if no matching column is found
        msymc <- grep("sym", tolower(colnames(ddata))) ##look for sym column
        
        if (length(msymc) == 0) {
            msymc <- grep("code", tolower(colnames(ddata))) ## if none, look for code column
        }
        
        if (length(msymc) == 0) {
            msymc <- grep("flag", tolower(colnames(ddata))) ## if still none, look for flag column
        }
        
        if (length(msymc) == 0) {
            msym <- rep(NA, length(ddata$Flow))  ## if not match, fill column with NA
        } else {
            msym <- ddata[,msymc]  ## if match, get values from correct column
        }
        
        ## format data frame for return
        ddata <- data.frame(ID=ddata[,mID], Date=mdates, Flow=ddata[,mFlow], SYM=msym, 
                            Agency=magency)
    }
    
    flow <- as.character(ddata$Flow)
    msym <- as.character(ddata$SYM)
    
    if (length(grep("_", flow)) != 0) {
        for (i in 1:length(flow)) {
            d.split <- unlist(strsplit(flow[i], "_", fixed=T))
            flow[i] <- d.split[1]
            msym[i] <- d.split[2]
        }
        ddata$Flow <- flow
        ddata$SYM <- msym
    }
    
    flow <- as.numeric(flow)
    
    Station <- ddata[1,1]
    if (is.na(ddata$Agency[1])) {
        if (Station %in% USGS.site.info$STAID) {ddata$Agency <- "USGS"}
    }
    
    ##replace no data values (e.g. -99) with NA
    ddata$Flow[ddata$Flow < 0] <- NA

    return(ddata)

}
