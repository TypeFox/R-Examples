## Function Description:
##     This function finds NAs in a vector, data.frame or matrix and returns a
##     matrix which contains two columns that includes the row number and column
##     number of each NA.

pin.na <- function(x=NULL){
    #----[ checking the input ]----#
    if(is.null(x)){
        stop("Please provide parameter \"x\". It can be a matrix, data.frame or vector.")
    }
    
    
    #----[ pre-processing ]----#
    # convert to character if it is vector
    if(class(x)=="factor"){
        x <- as.character(x)
    }
    
    
    #----[ processing ]----#
    # if x is either data.frame or matrix
    if(class(x) %in% c("data.frame", "matrix")){
        if(length(tmp <- which(is.na(x)))!=0){
            clmn <- (tmp %/% dim(x)[1])+1
            rw <- tmp %% dim(x)[1]
            clmn[rw==0] <- clmn[rw==0]-1
            rw[rw==0] <- dim(x)[1]

            output <- cbind(rw,clmn)
        }else{
            ## no NA was found
            output <- FALSE
        }
    }else if(class(x) %in% c("numeric", "integer", "character")){
        if(length(tmp <- which(is.na(x)))!=0){
            output <- tmp
        }else{
            ## no NA was found
            output <- FALSE
        }
    }else{
        # if the input is something other than classes above, we can not handle it!!
        stop("Please provide the d which can be a matrix, data.frame or vector.")
    }
    
    return(output)
}
