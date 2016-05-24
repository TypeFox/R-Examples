`annotate.samples` <-
function (x,y){

    if(!all(names(x) %in% c("expression", "background", "arraydescription", "Flags", "localization"))) {
        stop("The first argument is not a valid RPPA data list.")
    }

    ## order the sample information according the rownames of the expression matrix
    line<-match(rownames(x$expression),y[,1])
    sample.data<-y[line,]
   
    ## add the localization information
    sample.data <- cbind(sample.data, x$localization)
   
    ## set the column names
    colnames(sample.data) <- c(colnames(y), colnames(x$localization))
    
    ## assemble new RPPA data list containing all information, including the sample description
    rppa <- list(expression=x$expression,background=x$background,
            arraydescription=x$arraydescription,sampledescription=sample.data,Flags=x$Flags)
    
    
    ## print the head of the sample description
    print(head(sample.data))
    
    ## return data
    return(rppa)
}

