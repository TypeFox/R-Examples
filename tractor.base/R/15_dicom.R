#' The DicomMetadata class
#' 
#' This class represents DICOM metadata, which typically contains detailed
#' information about the scan parameters and subject.
#' 
#' @field source String naming the source file
#' @field tags Data frame of tag information
#' @field tagOffset Starting offset for tags in the file
#' @field dataOffset Starting offset for pixel data in the file
#' @field dataLength Pixel data length
#' @field explicitTypes Logical value indicating whether explicit types are
#'   used in the file
#' @field endian String naming the endianness of the file
#' @field asciiFields Character vector containing the contents of the ASCII
#'   header, if requested and present in the file.
#' 
#' @export
DicomMetadata <- setRefClass("DicomMetadata", contains="SerialisableObject", fields=list(source="character",tags="data.frame",tagOffset="integer",dataOffset="integer",dataLength="integer",explicitTypes="logical",endian="character",asciiFields="character"), methods=list(
    getAsciiFields = function (regex = NULL)
    {
        "Retrieve the value of one or more fields in the ASCII header. Returns NA if no fields match"
        if (is.null(regex))
            return (asciiFields)
        
        regex <- ore("^\\s*(\\S*", regex, "\\S*)\\s*=\\s*(.+)\\s*$")
        groups <- groups(ore.search(regex, asciiFields, simplify=FALSE), simplify=TRUE)
        if (all(is.na(groups)))
            return (NA)
        
        values <- groups[,ncol(groups),]
        if (all(values %~% ore(number)))
            values <- as.numeric(values)
        names(values) <- groups[,ncol(groups)-1,]
        return (values)
    },
    
    getDataLength = function () { return (dataLength) },
    
    getDataOffset = function () { return (dataOffset) },
    
    getEndianness = function () { return (endian) },
    
    getSource = function () { return (source) },
    
    getTags = function () { return (tags) },
    
    getTagOffset = function () { return (tagOffset) },
    
    getTagValue = function (group, element)
    {
        "Retrieve the value of a given tag, using an appropriate R type. Returns NA if the tag is missing"
        valueRow <- subset(tags, (tags$groups == group & tags$elements == element))
        if (dim(valueRow)[1] == 0 || valueRow$values == "")
            return (NA)
        else
        {
            value <- unlist(strsplit(as.vector(valueRow$values), "\\", fixed=TRUE, useBytes=TRUE))
            if (capabilities("iconv") == TRUE)
                value <- iconv(value, "", "LATIN1", sub="byte")
            value <- gsub("^\\s*(.+?)\\s*$", "\\1", value, perl=TRUE)
            
            if (as.vector(valueRow$types) %in% .Dicom$convertibleTypes)
                return (as.numeric(value))
            else
                return (value)
        }
    },

    nTags = function () { return (nrow(tags)) }
))

#' @export
print.DicomMetadata <- function (x, descriptions = FALSE, ...)
{
    tags <- x$getTags()
    nTags <- nrow(tags)
    
    if (nTags > 0)
    {
        if (descriptions)
        {
            cat("DESCRIPTION", rep(" ",19), "VALUE\n", sep="")
            for (i in seq_len(nTags))
            {
                description <- getDescriptionForDicomTag(tags$groups[i], tags$elements[i])
                cat(" ", substr(description, 1, 27), sep="")
                nSpaces <- max(3, 30-nchar(description))
                cat(rep(" ",nSpaces), sep="")
                cat(implode(x$getTagValue(tags$groups[i],tags$elements[i]), sep=", "))
                cat("\n")
            }
        }
        else
        {
            cat("GROUP    ELEMENT  VALUE\n")
            for (i in seq_len(nTags))
            {
                cat(sprintf(" 0x%04x   0x%04x   ", tags$groups[i], tags$elements[i]))
                cat(implode(x$getTagValue(tags$groups[i],tags$elements[i]), sep=", "))
                cat("\n")
            }
        }
    }
}

getDescriptionForDicomTag <- function (groupRequired, elementRequired)
{
    dictionaryRow <- subset(dictionary, (dictionary$group==groupRequired & dictionary$element==elementRequired))
    if (nrow(dictionaryRow) == 0)
        description <- sprintf("Unknown (0x%04x, 0x%04x)", groupRequired, elementRequired)
    else
        description <- as.character(dictionaryRow$description)
    
    return (description)
}

#' Read a DICOM file into a DicomMetadata object
#' 
#' This function reads a DICOM file into a \code{\link{DicomMetadata}} object.
#' Only DICOM files from magnetic resonance scanners are supported.
#' 
#' @param fileName The name of a DICOM file.
#' @param checkFormat If \code{TRUE}, the function will check for the magic
#'   string \code{"DICM"} at byte offset 128. This string should be present,
#'   but in reality not all files contain it.
#' @param stopTag An integer vector giving the group and element numbers (in
#'   that order) of a DICOM tag, or \code{NULL}. If not \code{NULL}, the
#'   function will stop parsing the DICOM file if the specified tag is
#'   encountered. This can be used to speed up the process if a specific tag is
#'   required.
#' @param ignoreTransferSyntax If \code{TRUE}, any transfer syntax stored in
#'   the file will be ignored, and the code will try to deduce the transfer
#'   syntax using heuristics. This may occasionally be necessary for awkward
#'   DICOM files, but is not generally recommended.
#' @param ascii If \code{TRUE}, the function will attempt to read an embedded
#'   Siemens ASCII header, if one exists.
#' @return \code{readDicomFile} returns a \code{\linkS4class{DicomMetadata}}
#'   object, or \code{NULL} on failure.
#' 
#' @author Jon Clayden
#' @seealso The DICOM standard, found online at \url{http://dicom.nema.org/}.
#'   (Warning: may produce headaches!) Also \code{\link{readDicomDirectory}}
#'   for information on how to create \code{\linkS4class{MriImage}} objects
#'   from DICOM files.
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @export
readDicomFile <- function (fileName, checkFormat = TRUE, stopTag = NULL, ignoreTransferSyntax = FALSE, ascii = TRUE)
{
    fileName <- expandFileName(fileName)
    
    if (!file.exists(fileName))
        report(OL$Error, "DICOM file ", fileName, " not found")
    
    # DICOM is sufficiently complicated that this can really only be interpreted to mean "probably" or "probably not"
    isDicomFile <- !checkFormat
    endian <- "little"
    explicitTypes <- TRUE
    deferredTransferSyntax <- NULL
    tagOffset <- 0
    dataOffset <- dataLength <- NA
    asciiFields <- NULL
    
    connection <- file(fileName, "rb")
    on.exit(close(connection))
    
    if (checkFormat)
    {
        # A well-formed DICOM file should contain a magic number at offset 128
        seek(connection, where=128)
        str <- rawToChar(stripNul(readBin(connection, "raw", n=4)))
        if (str == "DICM")
        {
            isDicomFile <- TRUE
            tagOffset <- 132
        }
        else
            seek(connection, where=0)
    }
    
    # Read the first group number
    group <- readBin(connection, "integer", n=1, size=2, signed=FALSE, endian=endian)
    
    # Some GE files have no magic number and no group 2, starting with group 8 instead
    if (checkFormat && !isDicomFile)
    {
        if (isTRUE(group == 0x0008))
            isDicomFile <- TRUE
        else if (isTRUE(group == 0x0800))
            isDicomFile <- TRUE
    }
    
    # Assume that the first group shouldn't be huge, to guess the endianness
    if (isTRUE(group > 0x00ff))
        endian <- setdiff(c("big","little"), endian)
    
    # Pretty fragile, this, since it assumes we know the type of the first tag
    # Seems to work in practice, but in theory we can assume implicit little-endian
    # after group 2 if there is no transfer syntax specified
    seek(connection, where=2, origin="current")
    type <- rawToChar(stripNul(readBin(connection, "raw", n=2)))
    if (isTRUE(type == "UL"))
        explicitTypes <- TRUE
    else
        explicitTypes <- FALSE
    
    if (isDicomFile)
    {
        groups <- numeric(0)
        elements <- numeric(0)
        types <- character(0)
        values <- character(0)
        
        sequenceLevel <- 0
        duplicateTags <- FALSE
        
        seek(connection, where=tagOffset)
        repeat
        {
            currentGroup <- readBin(connection, "integer", n=1, size=2, signed=FALSE, endian=endian)
            if (length(currentGroup) == 0)
                break
            else if (currentGroup == 0x0800)
            {
                currentGroup <- 0x0008
                endian <- setdiff(c("big","little"), endian)
            }
            currentElement <- readBin(connection, "integer", n=1, size=2, signed=FALSE, endian=endian)
            
            # Group 2 should always be explicit little-endian; the transfer syntax applies after that
            if (!ignoreTransferSyntax && !is.null(deferredTransferSyntax) && currentGroup > 2)
            {
                endian <- .Dicom$transferSyntaxes[[deferredTransferSyntax]]$endian
                explicitTypes <- .Dicom$transferSyntaxes[[deferredTransferSyntax]]$explicitTypes
                deferredTransferSyntax <- NULL
            }
            
            # Sequence related tags are always untyped
            if (currentGroup == 0xfffe)
            {
                # End of sequence delimiter
                if (sequenceLevel > 0 && currentElement == 0xe0dd)
                    sequenceLevel <- sequenceLevel - 1
                
                lengthSize <- 4
                type <- "UN"
            }
            else if (explicitTypes)
            {
                type <- rawToChar(stripNul(readBin(connection, "raw", n=2)))
                if (any(.Dicom$longTypes == type))
                {
                    seek(connection, where=2, origin="current")
                    lengthSize <- 4
                }
                else
                    lengthSize <- 2
            }
            else
            {
                lengthSize <- 4
                type <- dictionary$type[dictionary$group==currentGroup & dictionary$element==currentElement]
                if (length(type) == 0)
                    type <- "UN"
            }
            
            length <- readBin(connection, "integer", n=1, size=lengthSize, signed=ifelse(lengthSize > 2, TRUE, FALSE), endian=endian)
            
            report(OL$Debug, "Group ", sprintf("0x%04x",currentGroup), ", element ", sprintf("0x%04x",currentElement), ", type ", type, ", length ", length, ifelse(sequenceLevel>0," (in sequence)",""))
            
            if (any(c("OX","OW","OB","UN") == type) || (type == "SQ" && sequenceLevel > 0))
            {
                if ((currentGroup == 0x7fe0) && (currentElement == 0x0010))
                {
                    dataOffset <- seek(connection, where=NA)
                    dataLength <- length
                }
                
                if (type == "SQ" && length == -1)
                    sequenceLevel <- sequenceLevel + 1
                else if (length > 0)
                    seek(connection, where=length, origin="current")
                
                next
            }
            else if (any(groups==currentGroup & elements==currentElement))
            {
                duplicateTags <- TRUE
                if (length > 0)
                    seek(connection, where=length, origin="current")
                
                next
            }
            
            groups <- c(groups, currentGroup)
            elements <- c(elements, currentElement)
            types <- c(types, type)
            
            # Handle sequences of indeterminate length (to date only seen in Philips data)
            if (type == "SQ" && length == -1)
            {
                if (sequenceLevel == 0)
                    values <- c(values, "(Sequence)")
                sequenceLevel <- sequenceLevel + 1
                next
            }
            
            if (any(.Dicom$nonCharTypes$codes == type))
            {
                loc <- which(.Dicom$nonCharTypes$codes == type)
                size <- .Dicom$nonCharTypes$sizes[loc]
                nValues <- length/size
                
                if (.Dicom$nonCharTypes$rTypes[loc] == "integer")
                    value <- readBin(connection, "integer", n=nValues, size=size, signed=ifelse(size > 2, TRUE, .Dicom$nonCharTypes$isSigned[loc]), endian=endian)
                else
                    value <- readBin(connection, "double", n=nValues, size=size, endian=endian)
                
                if (nValues > 1)
                    values <- c(values, implode(value,sep="\\"))
                else if (nValues == 0)
                    values <- c(values, "")
                else
                    values <- c(values, as.character(value))
            }
            else
                values <- c(values, rawToChar(stripNul(readBin(connection, "raw", n=length))))
            
            if (!ignoreTransferSyntax && currentGroup == 0x0002 && currentElement == 0x0010)
            {
                deferredTransferSyntax <- values[length(values)]
                if (!(deferredTransferSyntax %in% names(.Dicom$transferSyntaxes)))
                    report(OL$Error, "Transfer syntax ", deferredTransferSyntax, " is not supported")
            }
            
            if (!is.null(stopTag) && currentGroup == stopTag[1] && currentElement == stopTag[2])
                break
        }
        
        if (duplicateTags)
            flag(OL$Warning, "Duplicated DICOM tags detected - only the first value will be kept")
        
        if (ascii)
        {
            asciiHeader <- ore.search(ore("### ASCCONV BEGIN[^#]+###(.+)### ASCCONV END ###",options="m"), ore.file(fileName,binary=TRUE))[,1]
            asciiFields <- ore.split(ore("\n",syntax="fixed"), asciiHeader)
        }
        
        tags <- data.frame(groups=groups, elements=elements, types=types, values=values, stringsAsFactors=FALSE)
        invisible (DicomMetadata$new(source=fileName, tags=tags, tagOffset=as.integer(tagOffset), dataOffset=as.integer(dataOffset), dataLength=as.integer(dataLength), explicitTypes=explicitTypes, endian=endian, asciiFields=as.character(asciiFields)))
    }
    else
        invisible (NULL)
}
