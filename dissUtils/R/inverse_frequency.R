inverse.frequency <- function(X){

    X <- as.factor(X);

    return(1/table(X)[X]);
}
