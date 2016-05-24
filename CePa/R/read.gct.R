read.gct <-
function (file) {
    expr = read.table(file, skip = 2, header = TRUE, sep = "\t", quote = "")
    rownames(expr) = expr[,1]

    checkName = table(expr[,1])
    if(max(checkName) > 1) {
        stop(paste("Genes in gct file should be unique: ", names(which.max(checkName)), sep = " "))
    }
    expr = expr[,-c(1,2)]
    expr = as.matrix(expr)
    
    return(expr)
}
