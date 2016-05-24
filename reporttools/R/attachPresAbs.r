attachPresAbs <- function(v){
    v <- as.numeric(v)
    res <- factor(v, levels = 0:1, labels = c("absent", "present"))
    return(res)
}
