belongsTo <- function(a, B){
    for(i in 1:length(B)) if(setequal(a, B[[i]])) return(i)
    return(FALSE)
}
