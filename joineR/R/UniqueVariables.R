UniqueVariables<-function (data, var.col, id.col = "ID") 
{
    patid <- unique(data[[id.col]])
    npat <- length(patid)
    n.cov <- length(var.col)
    n.col <- n.cov + 1
    new <- matrix(ncol = n.col, nrow = npat)
    new <- as.data.frame(new)
    if (is.numeric(id.col)) {
        names(new)[1] <- names(data)[id.col]
    }
    else {
        names(new)[1] <- id.col
    }
    if (is.numeric(var.col)) {
        names(new)[2:n.col] <- names(data)[var.col]
    }
    else {
        names(new)[2:n.col] <- var.col
    }
    id.col <- names(new)[1]
    var.col <- names(new)[2:n.col]
    new <- new[, c(1, order(names(new)[2:n.col]) + 1)]
    for (i in 1:n.col) {
        tt <- names(new)[i]
        if (class(data[[tt]]) == "factor") {
            new[[tt]] <- as.factor(new[[tt]])
            levels(new[[tt]]) <- levels(data[[tt]])
        }
        else {
            class(new[[tt]]) <- class(data[[tt]])
        }
    }
    var.col <- names(new)[2:n.col]
    new[[id.col]] <- patid
    for (i in 1:(npat)) {
        data.i <- data[data[[id.col]] == patid[i], names(data) %in% 
            var.col]
        if (length(var.col) == 1) {
           names(data.i) <- var.col
           data.i <- data.i[order(names(data.i))]
           tt.i <- unique(data.i)
        }
        else {
           data.i <- data.i[, order(names(data.i))]
           tt.i <- apply(data.i, 2, function(x) {
           unique(x)
           })
        }
        if (sum(unlist(lapply(tt.i, function(x) {
            length(x) > 1
        }))) > 0) {
            stop("No consistency on the variables information")
        }
        if (length(var.col) == 1){
            new[new[[id.col]] == patid[i], -which(names(new) == id.col)] <- data.i[1]
        }
        else {
        new[new[[id.col]] == patid[i], -which(names(new) == id.col)] <- data.i[1, ]
        }
    }
    return(new)
    cat("\n")
}
