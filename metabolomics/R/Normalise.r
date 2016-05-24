# require(crmn)
Normalise <-function(inputdata, 
    method=c("median", "mean", "sum", "ref", "is", "nomis", "ccmn", "ruv2"),
    refvec=NULL, ncomp=NULL, k=NULL, nc=NULL, saveoutput=FALSE, 
    outputname=NULL)
{
    # Match the method
    method <- match.arg(method)

    # If method is not one of the above listed, then stop
    if (!is.element(method, 
        c("median", "mean", "sum", "ref", "is", "nomis", "ccmn", "ruv2"))
    ) {
        stop("Invalid normalization method")
    }
    
    # Reference vector should be a vector
    if (!is.null(refvec)) {
        if (class(refvec) %in% c("data.frame", "list", "matrix")) {
            stop("Reference should be a vector")
        }
    }
                                     
    # If there is no refvec is given, get them to enter the reference vector
    if ((method == "ref"| method == "is") & is.null(refvec)) {
        stop("Please enter the reference vector")
    }
    # If no k, get them to enter it 
    if (method == "ruv" & is.null(k)) {
        stop("Please enter the number of unwanted variation factors")
    }
    # If there is no nc, get them to enter it
    if (method == "ruv2" & is.null(nc)) {
        stop(
            paste("Please enter a logical vector indicating",
                "the non-changing metabolites"
            )
        )
    }

    
    # If there is no nc, get them to enter it for ccmn
    if ((method == "ccmn"|method == "nomis") & is.null(nc)) {
        stop(
            paste("Please enter a logical vector indicating",
                "the internal standards"
            )
        )
    }
    
    if (method == "ccmn"){
        warning(paste("The ccmn method uses the grouping structure in the normalisation method, 
                      therefore, should not be used for those unsupervised 
                      methods where the groups must be treated as unknown."))      
    }
    
    if (method == "ruv2"){
      warning("The ruv2 method generates a matrix of unwanted 
                variation using LinearModelFit(). For identifying 
                differentially abundant metabolites, use LinearModelFit() 
                directly with ruv2=TRUE")     
    }
    
    # Get groups information
    Group <- inputdata[, 1]

    # Get normalisation vector according to the method and inputs, and
    # remove groups and internal standard for data processing

    # Remove groups for processing
    pre_norm <- inputdata[, -1]

    # Median vector
    if (method == "median") {
        norm_vector <- apply(pre_norm, 1, median, na.rm=TRUE)
    # Mean vector
    } else if (method == "mean") {
        norm_vector <- rowMeans(pre_norm, na.rm=TRUE)
    # Sum vector
    } else if (method == "sum") {
        norm_vector <- rowSums(pre_norm, na.rm=TRUE)
    } else if (method == "ref" | method == "is") {
        norm_vector <- refvec
    }
    
    # Prepare an empty matrix
    norm_data <- matrix(NA, nrow=nrow(pre_norm), ncol=length(pre_norm))
    rownames(norm_data) <- rownames(pre_norm)
    colnames(norm_data) <- colnames(pre_norm)

    if (!is.null(nc)){
        ncvec<-logical(ncol(pre_norm))
        ncvec[nc]<-TRUE      
    }
      
    if (method == "ccmn") {
        norm_data <- t(normalize(t(pre_norm), "crmn", 
            factor=model.matrix(~-1 + Group), standards=ncvec, ncomp=ncomp, 
            lg=FALSE)
        )
    } else if (method == "nomis") {
        norm_data <- t(normalize(t(pre_norm), "nomis", standards=ncvec, 
            lg=FALSE)
        )
    } else if (method == "ruv2") {
        norm_data<-LinearModelFit(datamat=data.matrix(pre_norm),
            factormat=model.matrix(~-1 + Group), 
            ruv2=TRUE,
            k=k, nc=nc)$uvmat
    } else {
        norm_data <- sweep(pre_norm, 1, norm_vector, "-")
    }
    
    # Reattach groups information
    outdata <- data.frame(Group, norm_data)
    #Edit column names
    outdata <- editcolnames(outdata)
    
    # Generate the output matrix in .csv format
    if (saveoutput) {
        write.csv(outdata,
            if (!is.null(outputname)) {
                paste(c(outputname, ".csv"), collapse="")
            } else {
                paste(c("normalized_", method, ".csv"), collapse="")
            }
        )
    }
    
    output <- list()
    output$output <- outdata
    output$groups <- Group
    output$samples <- row.names(inputdata)
    
    return(structure(output, class="metabdata"))
}
