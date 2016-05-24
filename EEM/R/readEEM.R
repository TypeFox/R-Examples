#' Read raw files and return a list
#' 
#' Read raw files from fluorescence spectrometer
#' 
#' @param path path to the files or folders which contains raw files (accept a vector). 
#' 
#' @return \code{readEEM} returns a list containing each raw files
#' 
#' @details The supported format is *.txt, *.csv and *.dat files from FP-8500 (JASCO), F-7000 (Hitachi Hi-tech), RF-6000 (Shimadzu) and
#' Aqualog (Horiba) fluorescence spectrometer. 
#' It is likely that outputs from different machines of the same companies are supported by this function.
#' Please send a word or pull request to add support for other formats. 
#' 
#' @export
#' 
#' @importFrom utils read.delim
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom R.utils isDirectory isFile
readEEM <-
    function(path = NULL){
        
        if (nargs() == 0) stop("Folder or file has not been specified.")
        
        # if file names are provided, use them
        if (all(isFile(path))) fileList <- path
        
        # if folder paths are provided, use them
        if (all(isDirectory(path))) {
            acceptableFileExtension <- "\\.csv$|\\.txt$|\\.dat$"
            fileList <- list.files(path = path, pattern = acceptableFileExtension, 
                                   ignore.case = TRUE, full.names = TRUE)
        }
        
        N = length(fileList)
        
        # stop if there is no file in the directory
        if (N == 0) stop("There is no files in the directory. ")
        
        ## initialize output
        EEM <- list() # store each sample data in list format. call each one using EEM[[i]]    
        
        ## get data from each file in loop
        for (i in 1:N){
            file <- fileList[i]
            EEM[[i]] <- readSingleEEM(file)
        }
        
        # delete those that were not meant to be read
        if (sum(grepl("DO NOT READ", EEM)) > 0) {
            idx <- grep("DO NOT READ", EEM)
            EEM <- EEM[-idx]
            fileList <- fileList[-idx]
        }
        
        # add name and class
        names(EEM) <- sapply(fileList, 
                             function(x) basename(file_path_sans_ext(x)), 
                             USE.NAMES = FALSE)
        class(EEM) <- "EEM"
        return(EEM)
    }

## create function to read each file 
readSingleEEM <- function(file){
    
    # check format of file  
    fileExtension <- tolower(file_ext(file))
    
    # read in the lines 
    tmpData = readLines(file, warn = FALSE)
    switch(fileExtension, csv = {SEP <- ","}, txt = {SEP <- ""}, dat = {SEP <- "\t"})

    # check if tmpdata contains non ASCII
    if (sum(grepl("UTF-8|unknown", sapply(tmpData, Encoding))) > 0) {
        tmpData <- sapply(tmpData, iconv, from = "UTF-8", to = "ASCII", sub = "byte")
        }
    
    # check for "EX/EM" (Shimadzu), if present transpose
    if (length(grep("<97><e3><8b>N<94>g<92><b7>/<8c>u<8c><f5><94>g<92><b7>", tmpData)) > 0) {
        toTranspose <- TRUE
    } else toTranspose <- FALSE
    
    # find the line index that contains either of the following word
    pat_FP8500 <- "XYDATA" # FP-8500 file: "XYData"
    # F-7000 file: "Data Points" or "Data list (in Japanese)"
    pat_F7000 <- "Data Points"
    pat_F7000_J_xls <- "<ef><be><83><ef><be><9e><ef><bd><b0><ef><be><80><ef><be><98><ef><bd><bd><ef><be><84>"
    pat_F7000_J_txt <- "<c3><de><b0><c0><d8><bd><c4>"
    pat_RF6000 <- "RawData|CorrectionData"
    pattern <- paste(pat_FP8500, pat_F7000, pat_F7000_J_xls, pat_F7000_J_txt, pat_RF6000, sep = "|")
    index <- grep(pattern, tmpData, ignore.case = TRUE)
    if (length(index) == 0) {
        if (fileExtension %in% "dat") { 
            # add exception for Aqualog 
            index <- 0
        } else {
            warning(paste0("'", basename(file), "' does not have the right format. So it will not be read."))
            data <- "DO NOT READ"
            return(data)   
        }
    }
    # read data

        # for txt or csv files
        data_noRowNames <- read.delim(file, sep = SEP, skip = index, 
                                      check.names = FALSE, row.names = NULL)
        
        # check if comma-delimited file was camouflaged as txt file
        if (fileExtension == "txt"){
            test <- read.delim(file, sep = ",", skip = index, 
                               check.names = FALSE, row.names = NULL)
            # replace data_noRowNames with test if test contains more columns
            if (dim(data_noRowNames)[2] < dim(test)[2]) {
                data_noRowNames <- test
            }
        }
    
    # output
    data <- as.matrix(data.frame(c(data_noRowNames[,-1]), 
                                 row.names = as.numeric(as.character(data_noRowNames[,1])), 
                                 check.names = FALSE)) 
    
    # delete NA or blank rows if present
    NA_rows <- which(is.na(data[,1])|data[,1] %in% "")
    if (length(NA_rows) > 0) data <- data[-NA_rows,]
    
    # delete NA or blank columns if present
    NA_col <- which(is.na(colnames(data))|colnames(data) %in% ""|grepl("NA", colnames(data)))
    if (length(NA_col) > 0) data <- data[,-NA_col]
    
    # make column names into numeric
    colnames(data) <- as.numeric(colnames(data))
    
    # transpose data if required
    if (toTranspose) data <- t(data)
    
    # sort columns so that the wavelength will continue to increase
    col_order <- order(as.numeric(colnames(data)))
    if (col_order[1] != 1) data <- data[, col_order]
    
    return(data)
}