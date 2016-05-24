compair <- function(x, predlog10mfi.var = "predicted.log10_mfi", 
                    residual.var = "residuals"){
    if( all(c(predlog10mfi.var,residual.var)%in%names(x)) ) {
        mat <- apply( x[,c(predlog10mfi.var,residual.var)], 1, 
                    function(y) is.na(y))
        ans <- apply(mat,2,function(y) sum(y))
        ans <- all(ans>=1)
    } else {
        ans <- TRUE
    }
    return(ans)
} 
