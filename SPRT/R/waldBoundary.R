waldBoundary <-
function(type1 = 0.05, type2 = 0.2, boundary = NULL, log = TRUE) {
    
    ###############
    # Verify inputs
    
    type1 <- as.numeric(type1)
    type2 <- as.numeric(type2)    
    
    if (!is.null(boundary)) {
        boundary <- as.character(boundary)
        
        if (length(boundary) != 1 || !(boundary %in% c("A","B"))) {
            stop("The boundary parameter is incorrect. Please enter either \"A\" or \"B\", or leave it empty to return both \"A\" and \"B\".")
        }
    }
    
    #################
    # Wald boundaries
    
    if (!is.null(boundary)) {
        
        if (boundary == "A") {
            output <- (1-type2)/type1
            
        } else if (boundary == "B") {
            output <- type2/(1-type1)
        }
        
    } else {
        output <- c((1-type2)/type1, type2/(1-type1))
        attr(output, "names") <- c("A", "B")
    }
    
    # Log
    if (log == TRUE) output <- log(output)
    
    # Return output
    output
}
