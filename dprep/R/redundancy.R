redundancy <-
function (data) 
{
    data = as.matrix(data)
    p = dim(data)[2]
    data = data[, -p]
    datau = unique(data)
    rowsu = dim(datau)[1]
    cat("\nnumber of distinct instances", rowsu)
    cat("\n")
    frequ = rep(0, rowsu)
    rowsrep = 1:rowsu
    for (j in 1:rowsu) {
        frequ[j] = length(row.matches(datau[j, ], data))
    }
    rowsrep = rowsrep[frequ > 1]
    if (length(rowsrep) == 1) {
        datarep = c(datau[rowsrep, ], frequ[rowsrep])
    }
    else {
        datarep = cbind(datau[rowsrep, ], frequ[rowsrep])
        colnames(datarep) = c(colnames(datau), "freq")
    }
    return(datarep)
}
