what.else <- function(nn) {
    if(is.pnn.trained && !is.pnn.smoothed) {
        message(paste(
            "PNN trained with", nn$n, "observations describing",
            length(nn$categories), "categories with", nn$k, "variables."
        ), sep=" ")
        message("Ready to be smoothed with function 'smooth'")
        return(TRUE)
    } else {
        message(
            "This object is not a Probabilist neural network from the package PNN."
        )
        message("Please, use another appropriate function.")
        return(FALSE)
    }
}

is.data.frame.with.real.variables <- function(df) {
    for( i in 1:length(df[1,]) ) { if( !is.numeric(df[,i]) ) return(FALSE) }
    return(TRUE)
}

is.trainingset.correct <- function(nn) {
    if(
        is.data.frame(nn)
        && is.factor(nn$set[,nn$category_column])
        && is.data.frame.with.real.variables(nn$set[,-nn$category_column])
    ) return(TRUE) else return(FALSE)
}

is.pnn.trained <- function(nn) {
    if(
        is.data.frame(nn)
        && !is.null(nn$model)
        && !is.null(nn$set)
        && !is.null(nn$category.column)
        && !is.null(nn$categories)
        && !is.null(nn$k)
        && !is.null(nn$n)
    ) return(TRUE) else return(FALSE)
}

is.pnn.smoothed <- function(nn) {
    if(
        is.pnn.trained(nn)
        && !is.null(nn$sigma)
    ) return(TRUE) else return(FALSE)
}

is.pnn.performed <- function(nn) {
    if(
        is.pnn.smoothed(nn)
        && !is.null(nn$observed)
        && !is.null(nn$guessed)
        && !is.null(nn$success)
        && !is.null(nn$fails)
        && !is.null(nn$success_rate)
        && !is.null(nn$bic)
    ) return(TRUE) else return(FALSE)
}