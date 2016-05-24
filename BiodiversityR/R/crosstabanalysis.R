`crosstabanalysis` <-
function(x,variable,factor){
    cross <- table(x[,variable]>0,x[,factor])
    result <- chisq.test(cross) 
    return(result)
}

