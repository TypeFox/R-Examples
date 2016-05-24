successes <-
function(data=NULL, model=NULL, timeout=NULL, addCosts=NULL) {
    if(is.null(data) || is.null(model)) {
        stop("Need both data and model to calculate successes!")
    }
    if(is.null(data$success)) {
        stop("Need successes to calculate successes!")
    }

    if(is.null(addCosts)) {
        ac = attr(model, "addCosts")
        if(is.null(ac) || ac == TRUE) {
            addCosts = TRUE
        } else {
            addCosts = FALSE
        }
    }

    hp = attr(model, "hasPredictions")
    if(is.null(hp) || hp != TRUE) {
        if(length(data$test) > 0) {
            predictions = rbind.fill(lapply(data$test, function(x) {
                data$data = data$data[x,]
                data$best = data$best[x]
                model(data)
            }))
        } else {
            predictions = model(data)
        }
    } else {
        predictions = model$predictions
    }

    if(is.null(timeout)) {
        # if timeout value wasn't given, assume maximum from data set
        timeout = max(data$data[data$performance])
    }

    if(addCosts || !is.null(data$costs)) {
        if(is.null(data$costGroups)) {
            usedFeatures = intersect(data$cost, sapply(data$features, function(x) { paste(x, "cost", sep="_") }))
        } else {
            usedFeatures = subset(data$cost, sapply(data$cost, function(x) { length(intersect(data$costGroups[[x]], data$features)) > 0 }))
        }
    }

    perfs = data$data[data$performance]
    successes = data$data[data$success]
    if(!addCosts || is.null(data$cost)) {
        costs = rep.int(0, nrow(perfs))
    } else {
        costs = apply(data$data[usedFeatures], 1, sum)
    }
    predictions$iid = match(do.call(paste, predictions[data$ids]), do.call(paste, data$data[data$ids]))
    predictions$pid = match(predictions$algorithm, data$performance)
    predictions$score = apply(predictions, 1, function(x) {
        pid = as.numeric(x[["pid"]])
        if(is.na(pid)) {
            FALSE
        } else {
            iid = as.numeric(x[["iid"]])
            score = as.numeric(perfs[iid,pid]) + costs[iid]
            if(!as.logical(successes[iid,pid]) || score > timeout) {
                FALSE
            } else {
                TRUE
            }
        }
    })
    agg = aggregate(as.formula(paste("score~", paste(c(data$ids, "iteration"), sep="+", collapse="+"))), predictions, function(ss) { ss[1] })
    agg$score
}
class(successes) = "llama.metric"
attr(successes, "minimize") = FALSE
