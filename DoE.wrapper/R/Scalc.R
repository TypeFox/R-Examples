Scalc <- function(design){
    if (!(is.matrix(design) | is.data.frame(design)))
        stop("design must be a matrix or a data frame")
    if (is.data.frame(design))
        if ("design" %in% class(design)) {
           rn <- response.names(design)
           design <- desnum(design)
           if (!is.null(rn)) design <- design[,-which(colnames(design) %in% rn)]
           }
        else design <- as.matrix(design)
    if (!is.numeric(design))
        stop("design must be a numeric matrix, or a data frame with numeric columns only.")
    choose(nrow(design),2)/sum(1/dist(design))
}