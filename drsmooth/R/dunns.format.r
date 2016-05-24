#' @name dunns.format
#' @title Format Dunn's Test Output
#' @param model  Dunn's test model object
#' @param test   Name of the test executed
#' @keywords internal

dunns.format <- function (model, test) {

    outputlist <- list()
    
    # Format output into matrices, 1 for the title, 1 for the results
    # Title
    title <- matrix(nrow=1, ncol=1)
    title[1,1] <- test
    
    # Results
    outputlist[[1]] <- title
    outputlist[[2]] <- model[[3]]
    return(outputlist)
}
