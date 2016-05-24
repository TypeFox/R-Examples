trainTest <-
function(data, trainpart = 0.6, stratify = FALSE) {
    assertClass(data, "llama.data")
    assertNumeric(trainpart)

    if(stratify) {
        stratifier = sapply(data$best, paste, collapse="-")
    } else {
        stratifier = rep.int(TRUE, nrow(data$data))
    }

    tmp = do.call(c, by(1:nrow(data$data), stratifier, function(x) {
        n = length(x)
        c(rep.int(1, round(n*trainpart)), rep.int(2, n-round(n*trainpart)))[sample(n, n)]
    }))
    parts = split(1:nrow(data$data), tmp)

    newdata = data
    newdata$train = list(parts[[1]])
    newdata$test = list(parts[[2]])
    attr(newdata, "hasSplits") = TRUE
    return(newdata)
}
