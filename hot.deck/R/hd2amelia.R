hd2amelia <-
function(object){
    class(object) <- "amelia"
    names(object)[which(names(object) == "data")] <- "imputations"
    return(object)
}
