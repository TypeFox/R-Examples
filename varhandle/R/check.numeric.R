## Function Description:
##     A function to assess if a vector can be interpreted as numbers
# * it can get a parameter to trim white space
# * it can get a parameter to return the numeric vector or just one simple TRUE or FALSE

check.numeric <- function(v=NULL, rm.na=TRUE){
    #----[ checking the input ]----#
    if(is.null(v)){
        stop("The parameter \"v\" is not defined. It can be character vector or factor vector.")
    }else if(!class(v) %in% c("character", "factor")){
        if(class(v) %in% c("numeric", "integer")){
            warning(paste("The input is already of class", class(v)))
        }else{
            stop("The parameter \"v\" can only be a character vector or factor vector.")
        }
        
    }
    
    if(rm.na!=TRUE & rm.na!=FALSE){
        stop("The parameter \"rm.na\" should be either TRUE or FALSE.")
    }
    
    
    #----[ pre-processing ]----#
    # convert to character if it is vector
    if(class(v)=="factor"){
        v <- as.character(v)
    }
    
    
    #----[ processing ]----#
    if(rm.na){
        # if it has some NAs
        if(any(is.na(v))){
            # remove NAs
            v <- v[-pin.na(v)]
        }
        
    }
    if(length(grep(pattern="^(-|\\+)?\\d+(\\.?\\d+)?$", x=v, invert=TRUE))==0){
        # everything is numbers, so change it to numeric
        output <- TRUE
    }else{
        output <- FALSE
    }
    
    # return the result
    return(output)
}