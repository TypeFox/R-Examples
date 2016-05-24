posPredictionInterval = function(krige_object, p = 95, 
	value =  median(krige_object$krige_output$var1.pred))
# This function gives the position of the "p"% prediction interval
# relative to "value". 
{
    if(!inherits(krige_object,"autoKrige")) stop("Invalid input object, object should be of class autoKrige")
   
    # It is assumed that the kriging error is distributed normally around the kriging prediction.
    p_in_tail = (100 - p) / 2             # The probability left in each tail (upper and lower)
                         
    # An assumption is made that the kriging prediction and variance are located in the
    # first and second column of krige_object.                         
                                      # Prediction                                                # Standard deviation
    krige_object$krige_output$lower = krige_object$krige_output[[1]] + qnorm(p_in_tail / 100)*sqrt(krige_object$krige_output[[2]])
    krige_object$krige_output$upper = krige_object$krige_output[[1]] + qnorm((100 - p_in_tail) / 100)*sqrt(krige_object$krige_output[[2]])
    
    position = function(upper,lower,value)
    # This function returns the status for each gridcell in nl_grid
    # "above", "below", "unclear" for the p% prediction interval
    {
        if(lower > value) return("higher")
        if(upper < value) return("lower")
        else return("not distinguishable")
    }
    krige_object$krige_output$position = factor(mapply(FUN = position, krige_object$krige_output$upper, krige_object$krige_output$lower,value))
    
    result = list(pos_prediction = krige_object$krige_output[c("upper","lower","position")], p = p, value = value)
    class(result) = c("posPredictionInterval","list")

    return(result)
}
