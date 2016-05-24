# A method which uses the k-Means algorithm to binarize a real-valued vector.
# <nstart> controls how many random sets should be chosen by the 'kmeans' method and <iter.max> is the
# maximum number of iterations that are allowed (see also the help for kmeans)
binarize.kMeans <- function(vect, nstart=1, iter.max=10, dip.test=TRUE, na.rm=FALSE){
    #some checks of the arguments
    
    if(!is.numeric(vect))
        stop("The input vector must consist of numerical values!")

    if (!na.rm && any(is.na(vect)))
      stop("Cannot binarize in the presence of NA values!")
    else
    if (na.rm)
    {
      vect <- vect[!is.na(vect)]
    }
    
    if (any(!is.finite(vect)))
     stop("Cannot binarize Inf values!")
        
    if(!is.numeric(nstart))
        stop("'nstart' must be numeric!")
    if(nstart < 0)
        stop("'nstart' must be >= 0!")
        
    if(!is.numeric(iter.max))
        stop("'iter.max' must be numeric!")
    if(iter.max < 0)
        stop("'iter.max' must be >= 0!")
    
    if(length(vect) < 3)
        stop("The input vector must have at least 3 entries!")
    
    if(length(unique(vect))==1)
        stop("The input vector is constant!")
  
    if (dip.test)
    {      
      p.value <- dip.test(vect)$p.value
    }
    else
      p.value <- as.numeric(NA)
    
    #start the standard kmeans method to do all the calulations
    km_res <- kmeans(vect, 2, nstart = nstart, iter.max = iter.max)
    
    #the center with greater value should get the binarized value 1, the other 0.
    if(km_res$centers[1] > km_res$centers[2]){
        binarizeddata <- abs(km_res$cluster - 2)
    }
    else{
        binarizeddata <- km_res$cluster - 1
    }
    
    #calculate the threshold as mean of the calculated centers
    #threshold <- min(km_res$centers) + dist(km_res$centers)[1] / 2
    #threshold <- mean(km_res$centers)
    threshold <- (max(vect[!as.logical(binarizeddata)]) + min(vect[as.logical(binarizeddata)])) / 2
    
    #put all computed results into a 'BinarizationResult' object and return it
    return(new("BinarizationResult",
            originalMeasurements = vect,
            binarizedMeasurements = as.integer(binarizeddata),
            threshold = threshold,
            p.value = p.value,
            method = "k-Means"))
}

