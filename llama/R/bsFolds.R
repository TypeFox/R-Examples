bsFolds <-
function(data, nfolds = 10L, stratify = FALSE) {
    assertClass(data, "llama.data")
    assertInteger(nfolds)

    if(stratify) {
        stratifier = sapply(data$best, paste, collapse="-")
    } else {
        stratifier = rep.int(TRUE, nrow(data$data))
    }

    trainIdxs = lapply(1:nfolds, function(i) {
        tmp = by(1:nrow(data$data), stratifier, function(x) {
            unique(sample(x, length(x), replace = TRUE))
        })
        as.integer(unlist(tmp))
    })

    newdata = data
    newdata$train = trainIdxs
    newdata$test = lapply(trainIdxs, function(x) { setdiff(1:nrow(data$data), x) })
    attr(newdata, "hasSplits") = TRUE
    return(newdata)
}
