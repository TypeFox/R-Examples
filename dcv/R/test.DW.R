test.DW <-
function(x, y){
    if(!length(x) == length(y)){
        stop("\"x\" and \"y\" have different length.")
    }
    summary.dw <- dwtest(formula = x ~ y)
    return(summary.dw)
}
