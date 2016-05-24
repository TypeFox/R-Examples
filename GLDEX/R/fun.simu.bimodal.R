"fun.simu.bimodal" <-
function(result1, result2, prop1, prop2, len=1000,no.test = 1000, param1, param2) 
{
    if (missing(prop2)) {
        prop2 <- 1 - prop1
    }
   
    no.1 <- round(len * no.test * prop1)
    sample.fitted1 <- rgl(no.1, result1[1], result1[2], result1[3], 
        result1[4], param1)
    sample.fitted2 <- rgl(len * no.test - no.1, result2[1], result2[2], 
        result2[3], result2[4], param2)
    sample.fitted <- split(c(sample.fitted1, sample.fitted2), 
        1:no.test)
return(sample.fitted)

}

