#'
#' createColumnIndices.R
#'


#'
#' grep the column names of the header to find the right columns
#'
createColumnIndices <- function(col.names, foreground=NULL, background=NULL) {

    # check if the needed columns are given
    # otherwise search for the  pattern "F.*Mean" and "B.*Mean"
    # so we get only the first occurence of it
    if(is.null(foreground)) {
        foreground <- "F.*Mean"
        cat("Warning: No foreground signal column specified. Taking the first column matching 'F.*Mean'.\n")
    }

    if(is.null(background)) {
        background <- "B.*Mean"
        cat("Warning: No foreground signal column specified. Taking the first column matching 'B.*Mean'.\n")
    }


    colIndices <- vector("list", length=7)
    names(colIndices) <- c("F", "B", "Block", "Row", "Column", "ID", "Flags")

    # grep for the column in the header
    # foreground and background can be specified, in case it is not, take the first match for the above patterns
    colIndices$F <- grep(paste("^\\s*", foreground, "\\s*$",sep=""), col.names, perl=TRUE)[1]
    colIndices$B <- grep(paste("^\\s*", background, "\\s*$",sep=""), col.names, perl=TRUE)[1]

    colIndices$Block <- grep("Block", col.names, perl=TRUE)
    colIndices$Row <- grep("Row", col.names, perl=TRUE)
    colIndices$Column <- grep("Column", col.names, perl=TRUE)
    colIndices$ID <- grep("ID", col.names, perl=TRUE)
    colIndices$Flags <- grep("Flags", col.names, perl=TRUE)

    # test if we found extactly one column for each colum name
    # since we take the first element for F and B it can be NA if no F and B column was found by the match
    testLength <- as.logical(sapply(colIndices, function(xx) length(xx)!=1 || is.na(xx)))

    if(any(testLength)) {
        
        colNames <- names(colIndices)[testLength]
        stop("For column(s) ", paste(colNames, collapse=","), " there is no hit or more than one in the data file.")
    
    }



    return(colIndices)

}
