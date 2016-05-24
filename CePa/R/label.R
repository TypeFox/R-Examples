
sampleLabel = function (label, treatment, control) {
    if (sum(label == treatment) == 0 | sum(label == control) == 0) {
        stop("Can not find treatment label or control label.")
    }
    res = list(label = label, treatment = treatment, control = control)
    class(res) = "sampleLabel"
    return(res)
}

.treatment = function(sl) {
    return(which(sl$label == sl$treatment))
}

.control = function(sl) {
    return(which(sl$label == sl$control))
}

.permutate = function(sl) {
    sl$label = sample(sl$label, length(sl$label), replace = FALSE)
    return(sl)
}

.factor = function(sl) {
    return(factor(sl$label))
}