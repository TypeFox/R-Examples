
WigRepeatAnalyzer <- function(filename, inputFolder, outputFolder="./", repeatMaskerFilePath, isRegex=FALSE) 
{
    
    ## Check that files to read exist
    
    # repeat_masker_file_path
    if(!file.exists(repeatMaskerFilePath))
    {
        stop("The argument repeatMaskerFilePath does not refer to a valid file, please check that the path is correct and that the permissions policy allows to read it")
    }
    
    # Can we find files with regex in input_folder ? 
    # (tested in C code)
    
    ## Create output folder safely
    .safeCreateFolder(outputFolder)
    
    ## Search for the file in case of regular expression is given
    if( isRegex){
        cat( paste( "\nSearching input file in ", inputFolder), "\n")
        file_list <- list.files( path=inputFolder, pattern=filename)
        cat( "|--File(s) found : \n")
        cat( "|-- ", file_list, "\n")
        if( length(file_list)!=1){
            warning( "|--Could not locate unique file with regex :",filename)
            stop()
        }
        true_filename  <-  file_list[[1]]
    }else{
        true_filename  <-  filename
    }
    
    ## Launch C function
    
    # Arguments : ( char** r_file_name_regex, char** r_input_folder, char** r_output_folder, char** r_repeat_masker_file_path, int* is_regex)
    returnValues <- .C("C_AnalyzeRepeat", true_filename, inputFolder, outputFolder, repeatMaskerFilePath)
    
    Sys.sleep(10)
    
    cat("\nLooking for base filename  in :", true_filename)
    index_regex  <-  regexpr( ".wig", true_filename)[1]
    
    if( index_regex > 0){
        wigFileBaseName  <-  substring(true_filename, 1, index_regex-1)
    }
    else{
        wigFileBaseName  <-  true_filename
    }
    
    cat("\nBase filename is :", wigFileBaseName, "\n")
    
    # ----------------------------------------------------------------------
    # DATA ON THE DISTRIBUTION OF REPEAT CLASSES
    # ----------------------------------------------------------------------
    
    
    print( "----------------------------------------------------------------------")
    print( "1. Computing distribution of repeat class signal coverage and weight")
    print( "----------------------------------------------------------------------")
    print(outputFolder)
    outputFolder  <-  as.character(outputFolder)
    print(outputFolder)
    
    # Define the file about repeat class signal coverage
    regex_coverage  <-  paste(wigFileBaseName,"_repeatClassCoverage\\.txt", sep="")
    print( paste( "Searching input files in ", outputFolder))
    file_list_class_coverage <- list.files( path=outputFolder, pattern=regex_coverage)
    print( "File(s) found : ")
    print( file_list_class_coverage)
    if( length(file_list_class_coverage)!=1){
        warning( paste( "Could not locate unique file with regex :", regex_coverage))
        stop()
    }
    file_class_coverage  <-  paste( outputFolder, file_list_class_coverage[[1]], sep="/")
    
    # Define the file about repeat class signal weight
    regex_weight  <-  paste(wigFileBaseName,"_repeatClassWeight\\.txt", sep="")
    file_list_class_weight <- list.files( path=outputFolder, pattern=regex_weight)
    print( "File(s) found : ")
    print( file_list_class_weight)
    if( length(file_list_class_weight)!=1){
        warning( paste( "Could not locate unique file with regex :", regex_weight))
        stop()
    }
    file_class_weight  <-  paste( outputFolder, file_list_class_weight[[1]], sep="/")
    
    # Read the file containing the information on repeat class signal coverage and weight
    cat( "\n|-- reading source files (class): '", file_class_coverage,"'",sep="")
    classCoverageList  <-  tryCatch({
        read.table( file=file_class_coverage, stringsAsFactors=FALSE)
    }, error = function(e) {
        cat("\n|-- |-- Error reading source file : ")
        print(e)
        return (NA)
    })
    
    cat( "\n|-- reading source files (weight): '", file_class_weight,"'",sep="")
    classWeightList  <-  tryCatch({
        read.table( file_class_weight)
    }, error = function(e) {
        cat("\n|-- |-- Error reading source file : ")
        print(e)
        return( NA)
    })
    
    # If data is available, compute the coverage and weigth statistics
    if( !is.na(classCoverageList) && !is.na(classWeightList)){
        
        # Define the output image
        base_output_name  <-  substring( file_class_coverage, first=0, last=nchar(file_class_coverage)-24)
        print( base_output_name)
        svg( paste( base_output_name, "ClassCoverageAndWeight.svg", sep="_"), width=14, height=9)
        par(mfrow = c(1,2))
        
        # Draw pie chart for signal coverage
        print( "|-- drawing chart for coverage")
        class_coverage <- as.numeric( classCoverageList$V3)
        signal_length <- as.numeric( classCoverageList$V4)
        total_signal_length <- as.numeric( classCoverageList$V5)
        slices_class <- 100 * class_coverage/signal_length
        class_total_coverage <- 100 * sum( class_coverage) / total_signal_length
        names_class <- paste( classCoverageList$V1, classCoverageList$V2, sep=" ")
        colors <- rainbow(length(names_class)/2)
        colors <- c(rbind(colors,colors))
        bp <- barplot( slices_class, xlim=c(0,25), names.arg=names_class, col=colors,las=1, main= "Repeat classes signal coverage", horiz=TRUE)
        text(slices_class, bp, paste( signif(slices_class,2),"%",sep=""), pos=4)
        message <-  paste( "Total coverage = ", paste( signif(class_total_coverage,3),"%",sep=""), sep="")
        mtext( message, side=3, cex=0.8)
        bp
        
        # Draw bar chart for signal weight
        print( "|-- drawing chart for weight")
        class_weight  <-  as.numeric( classWeightList$V3)
        signal_weight  <-  as.numeric( classWeightList$V4)
        total_signal_weight  <-  as.numeric( classWeightList$V5)
        slices_class  <-  100 * class_weight/signal_weight
        class_total_weight  <-  100 * sum( class_weight)/total_signal_weight
        names_class  <-  paste( classWeightList$V1, classWeightList$V2, sep=" ")
        colors  <-  rainbow(length(names_class)/2)
        colors<-c(rbind(colors,colors))
        bp <- barplot( slices_class, xlim=c(0,25), names.arg=names_class, col=colors,las=1, main= "Repeat classes signal weight", horiz=TRUE)
        text(slices_class, bp, paste( signif(slices_class,2),"%",sep=""), pos=4)
        message <-  paste( "Total weight = ", paste( signif(class_total_weight,3),"%",sep=""), sep="")
        mtext( message, side=3, cex=0.8)
        bp
        
        dev.off()
        
        # Release memory
        print( "|-- release memory")
        rm( classCoverageList)
        rm( classWeightList)
        gc()
    }
    
    
    # ----------------------------------------------------------------------
    # DATA ON THE DISTRIBUTION OF REPEAT FAMILY
    # ----------------------------------------------------------------------
    print( "")
    print( "----------------------------------------------------------------------")
    print( "2. Computing distribution of repeat family signal coverage and weight")
    print( "----------------------------------------------------------------------")
    
    # Define the file about repeat family signal coverage
    regex_coverage  <-  paste(wigFileBaseName,"_repeatFamilyCoverage\\.txt", sep="")
    print( paste( "Searching input files in ", outputFolder))
    file_list_family_coverage <- list.files( path=outputFolder, pattern=regex_coverage)
    print( "File(s) found : ")
    print( file_list_family_coverage)
    if( length(file_list_family_coverage)!=1){
        warning( paste( "Could not locate unique file with regex :", regex_coverage))
        stop()
    }
    file_family_coverage  <-  paste( outputFolder, file_list_family_coverage[[1]], sep="/")
    
    # Define the file about repeat family signal weight
    regex_weight  <-  paste(wigFileBaseName,"_repeatFamilyWeight\\.txt", sep="")
    file_list_family_weight <- list.files( path=outputFolder, pattern=regex_weight)
    print( "File(s) found : ")
    print( file_list_family_weight)
    if( length(file_list_family_weight)!=1){
        warning( paste( "Could not locate unique file with regex :", regex_weight))
        stop()
    }
    file_family_weight  <-  paste( outputFolder, file_list_family_weight[[1]], sep="/")
    
    # Read the file containing the information on repeat family signal coverage and weight
    cat( "\n|-- reading source files (class): '", file_family_coverage,"'",sep="")
    familyCoverageList  <-  tryCatch({
                read.table( file=file_family_coverage, stringsAsFactors=FALSE)
            }, error = function(e) {
                cat("\n|-- |-- Error reading source file : ")
                print(e)
                return( NA)
            })
    
    familyWeightList  <-  tryCatch({
                read.table( file=file_family_weight, stringsAsFactors=FALSE)
            }, error = function(e) {
                cat("\n|-- |-- Error reading source file : ")
                print(e)
                return(NA)
            })
    
    # If data is available, compute the coverage and weigth statistics
    if( !is.na( familyCoverageList) && !is.na( familyWeightList)){
    
        # Define the output image
        base_output_name  <-  substring( file_family_coverage, first=0, last=nchar(file_family_coverage)-25)
        print( base_output_name)
        svg( paste( base_output_name, "FamilyCoverageAndWeight.svg", sep="_"), width=17, height=17)
        par(mfrow = c(1,2))
        
        # Draw pie chart for signal coverage
        print( "|-- drawing chart for coverage")
        family_coverage  <-  as.numeric( familyCoverageList$V3)
        signal_length  <-  as.numeric( familyCoverageList$V4)
        total_signal_length  <-  as.numeric( familyCoverageList$V5)
        slices_family  <-  100 * family_coverage/signal_length
        family_total_coverage  <-  100 * sum( family_coverage) / total_signal_length
        names_family  <-  paste( familyCoverageList$V1, familyCoverageList$V2, sep=" ")
        colors  <-  rainbow(length(names_family)/2)
        colors<-c(rbind(colors,colors))
        bp <- barplot( slices_family, xlim=c(0,25), names.arg=names_family, col=colors,las=1, main= "Repeat families signal coverage", horiz=TRUE)
        text(slices_family, bp, paste( signif(slices_family,2),"%",sep=""), pos=4)
        message <-  paste( "Total coverage = ", paste( signif(family_total_coverage,3),"%",sep=""), sep="")
        mtext( message, side=3, cex=0.8)
        bp
        
        # Draw bar chart for signal weight
        print( "|-- drawing chart for weight")
        family_weight  <-  as.numeric( familyWeightList$V3)
        signal_weight  <-  as.numeric( familyWeightList$V4)
        total_signal_weight  <-  as.numeric( familyWeightList$V5)
        slices_family  <-  100 * family_weight/signal_weight
        family_total_weight  <-  100 * sum( family_weight)/total_signal_weight
        names_family  <-  paste( familyWeightList$V1, familyWeightList$V2, sep=" ")
        colors  <-  rainbow(length(names_family)/2)
        colors<-c(rbind(colors,colors))
        bp <- barplot( slices_family, xlim=c(0,25), names.arg=names_family, col=colors,las=1, main= "Repeat families signal weight", horiz=TRUE)
        text(slices_family, bp, paste( signif(slices_family,2),"%",sep=""), pos=4)
        message <-  paste( "Total weight = ", paste( signif(family_total_weight,3),"%",sep=""), sep="")
        mtext( message, side=3, cex=0.8)
        bp
        
        dev.off()
        
        # Release memory
        print( "|-- release memory")
        rm( familyCoverageList)
        rm( familyWeightList)
        gc()
    }
    
    cat("\n\nFinished.")
    
    # Return arguments values after execution, even if no changes were made on them
    return(returnValues)
}
