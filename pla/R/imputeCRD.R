.imputeCRD <- function(data,
                       response = "Response",
                       sample = "SampleStep", ...) {
    responses <- split(data[, c(response)], data[, sample])
    means <- lapply(responses, mean, na.rm = TRUE)
    names <- names(means)
    for (i in 1:length(names)) {
        name <- names[i]
        value <- means[i]
        index <- is.na(data[, response]) & data[, sample] == name
        if (any(index))
            data[index, response] <- value
    }
    return(data)
}
