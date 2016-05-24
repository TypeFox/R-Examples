misclassificationPenalties <-
function(data=NULL, model=NULL, addCosts=NULL) {
    if(is.null(data) || is.null(model)) {
        stop("Need both data and model to calculate misclassification penalties!")
    }

    if(attr(model, "hasPredictions")) {
        predictions = model$predictions
    } else {
        if(length(data$test) > 0) {
            predictions = rbind.fill(lapply(data$test, function(x) {
                data$data = data$data[x,]
                data$best = data$best[x]
                model(data)
            }))
        } else {
            predictions = model(data)
        }
    }

    optfun = if(data$minimize) { min } else { max }
    perfs = data$data[data$performance]
    opts = apply(perfs, 1, optfun)
    predictions$iid = match(do.call(paste, predictions[data$ids]), do.call(paste, data$data[data$ids]))
    predictions$pid = match(predictions$algorithm, data$performance)
    predictions$score = apply(predictions, 1, function(x) {
        pid = as.numeric(x[["pid"]])
        if(is.na(pid)) {
            0
        } else {
            iid = as.numeric(x[["iid"]])
            as.numeric(abs(as.numeric(perfs[iid,pid]) - opts[iid]))
        }
    })
    agg = aggregate(as.formula(paste("score~", paste(c(data$ids, "iteration"), sep="+", collapse="+"))), predictions, function(ss) { ss[1] })
    agg$score
}
class(misclassificationPenalties) = "llama.metric"
attr(misclassificationPenalties, "minimize") = TRUE
