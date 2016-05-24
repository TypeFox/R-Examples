## create good info table for variables
# http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session
# the total memory the variables are using.
# sorting, alphabetically and by memory consumption (acending or deceding)

var.info <- function(list="ALL", envir=.GlobalEnv, human.readable=TRUE, sortby="size", decreasing=TRUE, n=Inf){
    #----[ checking the input ]----#
    if(is.null(list)){
        stop("The parameter \"list\" should be defined. It should contain a list of variable names with the length of 1 or more.")
    }else if(list=="ALL"){
        # if user has selected nothing for list or "ALL", the function will consider it as ls() for provided environment
        list=ls(envir=envir)
    }
    
    # If the keep has a name that is not in ls(), show error with the location of bad variable name
    if(any(!is.element(list, ls(envir = as.environment(envir))))){
        bad_var_location <- which(!is.element(list, ls(envir = as.environment(envir))))
        stop(paste("All the items in the \"list\" should be a real existing variable.\nThe item number", bad_var_location, "is not among variables of the selected environment!\n", sep=" "))
    }
    
    # check if the sortby is among valid column names
    sortby.index <- match(sortby, c("name", "class", "size", "detail"))
    if(is.na(sortby.index) | length(sortby)!=1){
        stop("The column specified by the parameter \"sortby\" is not valid. Valid columns are:\n  name, class, size, detail\nOne column should be selected.")
    }
    
    # check human.readable input
    if(!human.readable %in% c(TRUE, FALSE)){
        stop("The parameter \"human.readable\" should ither be TRUE or FALSE.")
    }
    
    # check decreasing input
    if(!decreasing %in% c(TRUE, FALSE)){
        stop("The parameter \"decreasing\" should ither be TRUE or FALSE.")
    }
    
    # check the n input
    if(class(n)!="numeric" | n<1 | length(n)!=1){
            stop("The parameter \"n\" should be a single positive integer number.")
    }
    
    
    #----[ processing ]----#
    # define variables to fill in for loop
    var.size.readable <- vector()
    var.size.byte <- vector()
    var.class <- vector()
    var.detail <- vector()
    for(i in list){
        ## get the variable
        the.var <- get(i, envir = as.environment(envir))
        
        
        ## get the size
        # size in byte format
        var.size.byte <- c(var.size.byte, object.size(x=the.var))
        if(human.readable){
            # this will be in human readable format
            var.size.readable <- c(var.size.readable, format(object.size(x=the.var), units="auto"))
        }
        
        
        ## get the class
        var.class <- c(var.class, class(get(i, envir = as.environment(envir))))
        
        
        ## get dim for matrix/dataframe
        if(class(the.var) %in% c("data.frame", "matrix")){
            var.detail <- c(var.detail, paste("dimension:", paste0(dim(the.var), collapse=", ")))
        }else if(class(the.var) %in% c("integer", "numeric", "character", "factor", "logical")){
            var.detail <- c(var.detail, paste("length:", length(the.var)))
        }else{
            var.detail <- c(var.detail, NA)
        }
    }
    
    ## construct the output
    if(human.readable){
        output <- data.frame(name=list, class=var.class, size=var.size.readable, detail=var.detail)
    }else{
        output <- data.frame(name=list, class=var.class, size=var.size.byte, detail=var.detail)
    }
    
    ## sort based on user's preference
    if(sortby=="size" & human.readable){
        output <- output[order(var.size.byte, decreasing=decreasing), ]
    }else{
        output <- output[order(output[, sortby.index], decreasing=decreasing), ]
    }
    
    ## correct the rownames
    row.names(output) <- 1:nrow(output)
    
    ## select the number of desired rows
    if(is.finite(n)){
        #  convert to integer
        n <- as.integer(n)
        # if the n is larger than number of rows we have
        if(nrow(output)<n){
            n=nrow(output)
        }
        # select the columns
        output <- output[1:n, ]
    }
    
    ## return the output    
    return(output)
}
