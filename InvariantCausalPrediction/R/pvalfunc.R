pvalfunc <- function( x, y, test="normal"){

    if(is.function( test)){
        pval <- test(x,y)
    }else{
        if(test %in% c("ranks","ks")){
            if( test=="ranks"){
                pval <- 2*min( t.test(x,y)$p.value, var.test(x,y)$p.value)
            }else{
                pval <- ks.test(x,y)$p.value
            }
        }else{
            pval <- 2*min( t.test(x,y)$p.value, var.test(x,y)$p.value)
        }
    }
    return(pval)
}
