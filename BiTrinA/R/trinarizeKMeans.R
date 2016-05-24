# A method which uses the k-Means algorithm to trinarize a real-valued vector.
# <nstart> controls how many random sets should be chosen by the 'kmeans' method and <iter.max> is the
# maximum number of iterations that are allowed (see also the help for kmeans)
trinarize.kMeans <- function(vect, nstart=1, iter.max=10, dip.test=TRUE, na.rm=FALSE){
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
    km_res <- kmeans(vect, 3, nstart = nstart, iter.max = iter.max)
    
    #the center with greater value should get the binarized value 2, then 1, the other 0.
    cent <- order(km_res$centers)
    trinarizeddata <- km_res$cluster

    trinarizeddata[km_res$cluster == 1] <- which(cent ==1)
    trinarizeddata[km_res$cluster == 2] <- which(cent ==2)
    trinarizeddata[km_res$cluster == 3] <- which(cent ==3)

    trinarizeddata <- trinarizeddata-1

    #calculate the threshold as mean of the calculated centers
    threshold1 <- (max(vect[trinarizeddata == 0]) + min(vect[trinarizeddata == 1])) / 2
    threshold2 <- (max(vect[trinarizeddata == 1]) + min(vect[trinarizeddata == 2])) / 2
    
    #put all computed results into a 'TrinarizationResult' object and return it
    return(new("TrinarizationResult",
            originalMeasurements = vect,
            trinarizedMeasurements = as.integer(trinarizeddata),
            threshold1 = threshold1,
            threshold2 = threshold2,
            p.value = p.value,
            method = "k-Means"))
}

