# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

createtemplate = function (features, classes){
    if (nrow(features) != length(classes)) 
        stop("Formant and vowel dimensions do not match.")
    vs = levels(as.factor(classes))
    nvs = length(vs)
    means = matrix(0, nvs, ncol(features))
    for (i in 1:ncol(features)) means[, i] = tapply(features[, 
        i], classes, mean)
    rownames(means) = vs
    colnames(means) = paste("f", 1:ncol(features), sep = "")
    tmp = features
    for (i in 1:nvs) tmp[classes == vs[i], ] = features[classes == 
        vs[i], ] - matrix(means[i, ], nrow(tmp[classes == vs[i], 
        ]), ncol(means), byrow = TRUE)
    covariance = var(tmp)
    ranges = matrix (0,ncol(features),2)
    for (i in 1:ncol(features)) ranges[i,] = range (features[,i])
    output = list(classes = vs, means = means, covariance = covariance, ranges = ranges, territory = NULL)
    class(output) = "template"
    return(output)
}

