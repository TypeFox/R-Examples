correlated <- function(result, limit=0.85, plot.scatter=FALSE, keep=NA){

    answer <- list()
    answer$measures <- correl(result$measures, plot.scatter=plot.scatter, limit=limit, keep=keep)
    answer$measures.uniform <- correl(result$measures.uniform, plot.scatter=plot.scatter, limit=limit, keep=keep)
    return(answer)
}

