## Wrap for test goodnes of fit as it can produce an error
get_modelfit <- function(object, x, method, grouping) {
    ans <- try(get_neilltest(object, x, method, grouping),silent=TRUE)
    if(inherits(ans,"try-error")){
        ans <- NA
    } else{
        if(is.character(ans)){
           if(ans=="error1") ans <- "Too many groups in 'grouping'"
           if(ans=="error2") ans <- "Too few groups in 'grouping'"
        } else{
          ans <- as.numeric(round(ans[,2],3))  
       }        
    }
    return(ans)
}

