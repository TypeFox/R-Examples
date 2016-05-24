KFparest <- function(   data, # data ie Y 
                        ...){ # delete and paste in OTHER ARGUMENTS TO BE PASSED TO KFfit
    start <- Sys.time()
    
    inits <- PUT INITIAL VALUES FOR PARAMETER VECTOR HERE
    
    # use optim to find optimal parameters
    oppars <- optim(inits,
                    KFfit,
                    data=data,
                    OTHER ARGUMENTS TO BE PASSED TO KFfit,  # delete and paste in OTHER ARGUMENTS 
                                                            # TO BE PASSED TO KFfit
                    optim=TRUE,
                    control=list(trace=100))    
    
    end <- Sys.time()
    cat("\n")
    cat("Time Taken",difftime(end,start,units="mins"),"\n")
    cat("\n")
    
    return(oppars)
}

