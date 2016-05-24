"fun.moments.bimodal" <-
function (result1, result2, prop1, prop2, len=1000,no.test = 1000, param1, param2) 
{
    if (missing(prop2)) {
        prop2 <- 1 - prop1
    }
   
sample.fitted <- fun.simu.bimodal(result1, result2, prop1, prop2, len,no.test, param1, param2)

  result.o <- lapply(1:no.test, function(i, sample.fitted) fun.moments( 
        sample.fitted[[i]]), sample.fitted)

result.o<-do.call("rbind",result.o)
dimnames(result.o)<-list(NULL,c("mean", "variance", "skewness", "kurtosis"))
return(result.o)

}

