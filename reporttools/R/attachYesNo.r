attachYesNo <- function(v){
    res <- factor(v, levels = 0:1, labels = c("no", "yes"))
    return(res)
}
