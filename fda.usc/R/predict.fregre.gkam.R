################################################################################
predict.fregre.gkam<-function (object, newx = NULL, type = "response", ...)
{
    namesx <- names(object$result)
    nvars = length(namesx)
    nr = nrow(newx[[namesx[1]]])
    pr = matrix(NA, nrow = nr, ncol = nvars + 3)
    colnames(pr) = c(colnames(object$effects), "eta", "mu")
    pr[, "Intercept"] = rep(object$effects[1, "Intercept"], nr)
    for (i in 1:nvars) {
      pr[, namesx[i]] = predict(object$result[[namesx[i]]],newx[[namesx[i]]])
    }
    if (nr == 1)
        pr[, "eta"] <- sum(pr[, 1:(nvars + 1)])
    else pr[, "eta"] = rowSums(pr[, 1:(nvars + 1)]) 
#pr[, "eta"] = apply(pr[, 1:(nvars + 1)], 1, sum)
    pr <- switch(type, response = object$family$linkinv(pr[,
        "eta"]), link = pr[, "eta"])
    return(pr)
}





