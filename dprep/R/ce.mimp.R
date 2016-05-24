ce.mimp <-
function (w.cl, method = c("mean", "median"), atr, nomatr = 0) 
{
    p = dim(w.cl)
    index.na = which(is.na(w.cl), arr.ind = TRUE)
    o = order(index.na[, 1], index.na[, 2])
    index.na = index.na[o, ]
    dimnames(index.na) = NULL
    var.na = sort(as.numeric(names(table(index.na[, 2]))))
    var.na = var.na[var.na %in% atr]
    if (length(var.na) == 0) 
        stop("Error: No missing values occur in relevant variables!")
    index.atr = matrix(index.na[index.na[, 2] %in% atr[]], , 
        2)
    class.na = as.matrix(w.cl[index.atr[, 1], p[2]])
    dimnames(class.na) = NULL
    class.na = cbind(index.atr, class.na)
#print(table(w.cl[index.na[,1],p[2]]))
#print(dim(index.na))
#    classes = sort(as.numeric(names(table(w.cl[index.na[, 1],p[2]]))))
#print(classes)
classes=table(w.cl[index.na[,1],p[2]])
#    num.class = length(classes)
    replace.na = rep(0, 0)
    for (i in 1:dim(class.na)[1]) {
        sub = w.cl[w.cl[, p[2]] == class.na[i, 3], ]
        tempo=as.numeric(class.na[i,2])
        if (tempo %in% nomatr) 
            imput.col = moda(sub[, tempo])[1]
        else if (method == "mean") 
            imput.col = mean(sub[, class.na[i, 2]], na.rm = TRUE)
        else if (method == "median") 
            imput.col = median(sub[, class.na[i, 2]], na.rm = TRUE)
        replace.na = rbind(replace.na, imput.col)
    }
    dimnames(replace.na) = NULL
    class.na = cbind(class.na, replace.na)
#print(dim(class.na))
    for (i in 1:dim(class.na)[1])
{ i1=as.numeric(class.na[i,1])
i2=as.numeric(class.na[i,2])
w.cl[i1, i2] = class.na[i, 4]}
    cat("\nSummary of imputations using substitution of ", method, 
        "(mode for nominal features):\n")
    colnames(class.na) = c("Row", "Column", "Class", "Imput.value")
    print(class.na)
    cat("\nTotal number of imputations per class: \n")
print(classes)
#    for (i in classes) {
#        amount = sum(class.na[, 3] == i)
#        cat("Class ", i, ": ", amount, "\n")
#    }
    cat("\nTotal number of imputations: ", dim(class.na)[1], 
        "\n")
    return(w.cl)
}
