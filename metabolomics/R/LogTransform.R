LogTransform <- function(inputdata, base=exp(1), saveoutput=FALSE,
    outputname="log.results")
{
    Group <- inputdata[, 1]
    
    # Remove groups for data processing
    prelog_data <- inputdata[, -1]
    
    #    Log transform the data
    log_data <- log(prelog_data, base)
    
    # Reattach groups information
    outdata <- cbind(Group, log_data)
    if (saveoutput) {
        write.csv(output, paste(c(outputname, ".csv"), collapse=""))
    }
    
    output <- list()
    output$output <- outdata
    output$groups <- Group
    output$samples <- row.names(inputdata)

    return(structure(output, class="metabdata"))
}
