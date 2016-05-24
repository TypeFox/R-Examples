# Miscellaneous common functions and data

#### Data ####

#' Example Change-O database
#'
#' A small example database subset from Laserson and Vigneault et al, 2014.
#'
#' @format   A data.frame containing Change-O style data.
#' 
#' @references
#' \enumerate{
#'   \item  Laserson U and Vigneault F, et al. High-resolution antibody dynamics of 
#'            vaccine-induced immune responses. 
#'            Proc Natl Acad Sci USA. 2014 111:4928-33.
#' }
#'
#' @seealso  See \code{\link{HS1FDistance}} for the human 1-mer distance matrix.
"InfluenzaDb"


#### Transformation functions ####

# Converts a matrix to a vector
#
# \code{clearConsole} clears the console.
# 
# @examples
# # Generate a sampel mutations_array 
# sample_matrix <- matrix(sample(20,4),nrow=2, dimnames=list( c("CDR","FWR"), c("R","S") ))
# collapseMatrixToVector(sample_matrix)
#
collapseMatrixToVector <- function(mat, byrow = FALSE){    
    # Get the row and column names
    rnames <- rownames(mat)
    cnames <- colnames(mat)
    if (is.null(rnames)) { rnames <- paste0("Row", 1:nrow(mat)) }
    if (is.null(cnames)) { cnames <- paste0("Column", 1:ncol(mat)) }
    # Combine the row and columns names
    cobminedNames <- outer(rnames, cnames, paste, sep = "_")
    
    # Collapse the matrix to a vector
    if (byrow) {
        collapsed_mat <- c(t(mat))
        names(collapsed_mat) <- c(t(cobminedNames))
    }
    else{
        collapsed_mat <- c(mat)
        names(collapsed_mat) <- c(cobminedNames)
    }
    return(collapsed_mat)
}

#### Logging and error checking functions ####

# Make a progress bar when using %dopar% with a parallel backend
# Adapated from http://stackoverflow.com/questions/5423760/how-do-you-create-a-progress-bar-when-using-the-foreach-function-in-r
#
# @return  NULL
doparProgressBar <- function(n){
    pb <- txtProgressBar(min=1, max=n-1,style=3)
    count <- 0
    function(...) {
        count <<- count + length(list(...)) - 1
        setTxtProgressBar(pb,count)
        #Sys.sleep(0.01)
        flush.console()
        rbind(...)
        
    }
}


# Check data.frame for valid columns and issue message if invalid
#
# @param   data    data.frame to check
# @param   columns  vector of column names to check
# @param   logic   one of "all" or "any" controlling whether all or at least one of
#                  the columns must be valid
# @return  TRUE is columns are valid and a string message if not.
checkColumns <- function(data, columns, logic=c("all", "any")) {
    # Check arguments
    logic <- match.arg(logic)
    
    data_names <- names(data)
    if (logic == "all") {
        # Check that all columns exist
        for (f in columns) {
            if (!(f %in% data_names)) { 
                msg <- paste("The column", f, "was not found") 
                return(msg)
            }
        }        
        # Check that all values are not NA
        for (f in columns) {
            if (all(is.na(data[[f]]))) { 
                msg <- paste("The column", f, "contains no data") 
                return(msg)
            }
        }
    } else if (logic == "any") {
        # Check that columns exist
        if (!any(columns %in% data_names)) {
            msg <- paste("Input must contain at least one of the columns:", paste(columns, collapse=", "))
            return(msg)
        }
        # Check that all values are not NA
        invalid <- sapply(columns, function(f) all(is.na(data_names[[f]])))
        if (all(invalid)) { 
            msg <- paste("None of the columns", paste(columns, collapse=", "), "contain data") 
            return(msg)
        }
    }
    
    # Return TRUE if all checks pass
    return(TRUE)
}


#### OS interaction functions ####

#' Determines the OS platform being used
#'
#' @return  The OS platform.
#' 
#' @examples
#' getPlatform()
#' 
#' @export
getPlatform <- function() {
    return(.Platform$OS.type)
}


#' Determines the numbers of CPU cores available
#'
#' @return  The number of cores available. Returns 1 if undeterminable.
#'
#' @examples
#' getnproc()
#' 
#' @export
getnproc <-function(){
    platform <- getPlatform()
    if (platform == "windows") {
        nproc <- as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS'))
    } else if (platform == "unix") { 
        nproc <- parallel::detectCores()
    } else {
        nproc <- 1    
    }

    return(nproc)
}


#' Clears the console
#'
#' \code{clearConsole} clears the console.
#' 
#' @examples
#' clearConsole()
#'
#' @export
clearConsole <- function(){
    platform <- getPlatform()
    if (platform == "windows") {
        cat("\014")
    } else if (platform == "unix") {
        system('clear')
    }
}
