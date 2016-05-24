#############################
#RNA-Seq results reading
#Reads wig files
#Romain Fenouil - 21/03/2009
#Nicolas Descostes - 16/02/2010
#############################

# #########################################################################################
# Function that permits to read a WIG file and build a R structure from it
# fileName = the path to the WIG file to read
# #########################################################################################

readWIG  <-  function(fileName){
    
    #Return the wig File as a dataFrame object
    df <- read.table(fileName,sep="^",header=FALSE,quote="")
    
    #Converts the dataFrame into a vector of characters/strings
    charData <- as.character(df$V1)
    
    #Retrieves the index of each  line containing the informations about a chromosome
    indexTrack <- grep("fixedStep chrom=chr",charData)
    
    #Retrieves each line containing the name of a chromosome
    lines <- charData[indexTrack]
    
    #Add the last index of the last chromosome data. +2 is added regarding the line 23 of the script
    indexTrack <- c(indexTrack,length(charData)+2)
    
    #Retrieves the names of the chromosomes, "\\1" is retrieving what is between the parenthesis of the regular expression
    chromNames <- gsub("fixedStep chrom=chr(\\w+) start=\\d+ step=\\d+","\\1",lines,perl=TRUE)      
    
    #For each chromosome, retrieves the probe values the levels of the list are the names of the chromosomes.
    wigList <- list()
    num <- 0
    
    for(i in 1:length(chromNames)){
        num <- num + 1
        startOfChromi <- (indexTrack[i]+1) #The beginning of a chromosome list of data is the first line after the chromosome declaration
        endOfChromi <- (indexTrack[i+1]-2) #The end of a chromosome is the beginning of the next one minus the track line declaration (therefore, the last data is two lines above the next chromosome declaration)
        wigList[[chromNames[i]]] <- as.numeric(charData[startOfChromi:endOfChromi])               #Puts the data of a chromosome in a list which level will be the name of the chromosome
        #print(paste(num,") chromosome",chromNames[i]," :",startOfChromi," -> ",endOfChromi,sep=""))
    }
    return(wigList)
}


# #########################################################################################
# Function that permits to write a WIG structure to file
# wigData = the WIG structure to write
# fileName = the base name of the WIG file to write (name without extension)
# folder = the path to the output folder
# fixedStep = the size of the bins in the WIG
# #########################################################################################

writeWIG  <-  function(wigData, fileName, folder=NULL, fixedStep=50, addExtension=TRUE)
{
    # Check if the extension has to be added
    if(addExtension && (!grepl(".*\\.wig$", fileName[1], ignore.case=TRUE, perl=TRUE))) fileName <- paste(fileName,".wig",sep="")
    
    # creating output file
    output_file_name  <-  ifelse(is.null(folder), fileName, file.path(folder, fileName))
    cat( "Writing file : ", output_file_name, "\n")
    file.create( output_file_name)
    fileCon <- file( output_file_name,open="a")
    
    for(chr in names(wigData))
    {
        trackLine <- paste("track type=wiggle_0 name=\"",fileName,"\"",sep="")
        writeLines(trackLine, con=fileCon)
        
        descLine <- paste("fixedStep chrom=chr", chr, " start=1 step=", fixedStep, sep="")
        writeLines(descLine, con=fileCon)
        writeLines(format(wigData[[chr]], scientific=FALSE, digits=6, trim=TRUE, drop0trailing = TRUE), con=fileCon)
    }
    
    close(fileCon)
}


