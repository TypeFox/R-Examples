gwplot <-
function (x, rep = NULL, max = NULL, min = NULL, file = NULL, 
    selected.covariate = 1, selected.response = 1, highlight = FALSE, 
    type = "p", col = "black", ...) 
{
    net <- x
    if (is.null(net$generalized.weights)) 
        stop("generalized weights were not calculated", call. = F)
    if (!is.null(file)) {
        if (!is.character(file)) 
            stop("'file' must be a string")
        if (file.exists(file)) 
            stop(sprintf("%s already exists", sQuote(file)))
    }
    if (!is.numeric(selected.covariate)) 
        for (i in 1:length(net$model.list$variables)) if (net$model.list$variables[i] == 
            selected.covariate) 
            selected.covariate = i
    if (!is.numeric(selected.covariate) || selected.covariate < 
        1 || selected.covariate > ncol(net$covariate)) 
        stop("'selected.covariate' does not exist")
    if (!is.numeric(selected.response)) 
        for (i in 1:length(net$model.list$response)) if (net$model.list$response[i] == 
            selected.response) 
            selected.response = i
    if (!is.numeric(selected.response) || selected.response < 
        1 || selected.response > ncol(net$response)) 
        stop("'selected.response' does not exist")
    if (!is.null(rep)) {
        if (rep == "best") 
            rep <- as.integer(which.min(net$result.matrix["error", 
                ]))
        if (length(net$generalized.weights) < rep) 
            stop("'rep' does not exist")
    }
    covariate <- as.vector(net$covariate[, selected.covariate])
    variablename <- net$model.list$variables[selected.covariate]
    column <- (selected.response - 1) * ncol(net$covariate) + 
        selected.covariate
    if (is.null(rep)) {
        matrix <- as.matrix(sapply(net$generalized.weights, function(x) rbind(x[, 
            column])))
        item.to.print <- min(which.min(net$result.matrix["error", 
            ]))
    }
    else {
        highlight = F
        matrix <- as.matrix(net$generalized.weights[[rep]][, 
            column])
        item.to.print <- 1
    }
    if (is.null(max)) 
        max <- max(matrix)
    if (is.null(min)) 
        min <- min(matrix)
    ylim <- c(min, max)
    if (!highlight || item.to.print != 1 || ncol(matrix) == 1) 
        plot(x = covariate, y = matrix[, 1], ylim = ylim, xlab = variablename, 
            ylab = "GW", type = type, col = col, ...)
    else plot(x = covariate, y = matrix[, 2], ylim = ylim, xlab = variablename, 
        ylab = "GW", type = type, col = col, ...)
    if (ncol(matrix) >= 2) {
        for (i in 2:ncol(matrix)) if (!highlight || (i != item.to.print)) 
            lines(x = covariate, y = matrix[, i], type = type, 
                col = col, ...)
    }
    if (highlight) {
        lines(x = covariate, y = matrix[, item.to.print], type = type, 
            col = "red", ...)
        legend("topright", paste("Minimal Error: ", round(net$result.matrix["error", 
            item.to.print], 3), sep = ""), col = "red", ...)
    }
    title(paste("Response: ", net$model.list$response[selected.response], 
        sep = ""))
    if (!is.null(file)) {
        weight.plot <- recordPlot()
        save(weight.plot, file = file)
    }
}
