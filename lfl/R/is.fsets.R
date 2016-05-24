is.fsets <- function(x) {
    return(inherits(x, 'fsets') && 
           is.matrix(x) &&
           is.vector(vars(x)) &&
           ncol(x) == length(vars(x)) &&
           all(colnames(x) == names(vars(x))) &&
           is.matrix(specs(x)) &&
           ncol(x) == ncol(specs(x)) &&
           ncol(x) == nrow(specs(x)) &&
           all(colnames(x) == colnames(specs(x))) &&
           all(colnames(x) == rownames(specs(x))))
}
