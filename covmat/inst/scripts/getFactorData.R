rm(list = ls())

library(quantmod)

#' Function to download factor data. This is a utility function and can be used
#' to update factor data in the data folder for the package. After executing the
#' function store the factor.data object in the data folder.
#' 
#' @details
#' This method downloads monthly data for factors that can be used for FMMC. FMMC needs
#' a risk model. We will use the 5 factors form Fama-French website. Momentum is 
#' also available from Ken French's website. We will add liquidity factor as pointed
#' out in Pastor and Stambaugh (JPE 2003). We will also include a volatility factor
#' with difference in VIX index. We will use data past 1996-04-01.
#' @author Rohit Arora
#' 
get.factor.data <- function() {
    
    # setup temp folder and extract
    temp.folder = paste(getwd(), 'temp', sep='/')
    dir.create(temp.folder, F)
    temp.folder = paste(getwd(), '/', 'temp', sep='')
    shell('del /F /S /Q temp\\*.*', wait = TRUE)
    
    # 1 Fama-French Data
    filenames <- c("F-F_Research_Data_5_Factors_2x3", "F-F_Momentum_Factor")
    data <- sapply(filenames, function(filename) {
 
        filename.txt <- paste(filename,".txt", sep="")
        filename.zip <- paste(filename,".zip", sep="")
        
        url = paste("http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/", 
                    filename.zip, sep="")
        
        download.file(url, filename)
        files = unzip(filename, exdir=temp.folder)    
        
        # read data
        file = paste(temp.folder, '/', filename.txt, sep='')
        
        out <- readLines(file)
        index <- tail(head(which(nchar(out) == 0),3),2)
        if(is.na(index[2])) index[2] <- length(out) + 1
        
        linesAtStart <- index[1]; linesForHeader <- 1
        linesToRead <- (index[2] - 1) - linesAtStart - linesForHeader    
        
        data.ff <- read.table(file, skip=index[1], header=T, nrows=linesToRead)
        xts(data.ff, as.Date(paste(rownames(data.ff),01,sep=""), format="%Y%m%d"))        
    })
    
    ff.factors <- na.omit(do.call(cbind,data))/100
    
    # 2 Liquidity Factor data based on Stambaugh 2003
    out <- read.table(paste("http://faculty.chicagobooth.edu/lubos.pastor/research/",
                      "liq_data_1962_2014.txt",sep=""), skip = 14) 
    liq <- xts(out[,4],as.Date(paste(out[,1],01,sep=""), format="%Y%m%d"))
    
    factors <- cbind(ff.factors, LIQ=liq)
    factors <- na.omit(factors)
    
    # 3 Volatility Factor as difference in vix
    getSymbols("^VIX", adjust=TRUE, from ="1996-04-01",to=tail(index(factors),1))
    vix    <- to.monthly(get("VIX"), indexAt='endof', drop.time=FALSE)
    assign(x="VIX", value=vix, envir =  .GlobalEnv)
    
    vix   <- diff(Cl(get("VIX", envir =  .GlobalEnv))); colnames(vix) <- "VOL"
    
    dates.vix.monthly <- format(index(vix), "%Y%m")
    dates.factors.monthly <- format(index(factors), "%Y%m")
    index(vix) <- index(factors)[which(dates.factors.monthly %in% dates.vix.monthly)]

    factors <- cbind(factors, vix)
    factors <- na.omit(factors)
    
    return( factors )
}

factor.data <- get.factor.data() 
remove(VIX); remove(get.factor.data)