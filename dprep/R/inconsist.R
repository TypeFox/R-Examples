inconsist <-
function (data) 
{
    data = as.matrix(data)
    dimnames(data) = NULL
    n = dim(data)[1]
    p = dim(data)[2]
    if (p == 2) {
        cardi = table(data[, 1])
        unicos = unique(data[, 1])
        mfrecu = max(table(data[data[, 1] == unicos[1], 2]))
        for (j in unicos[-1]) {
            tempo = max(table(data[data[, 1] == j, 2]))
            mfrecu = c(mfrecu, tempo)
        }
        mfrecu = mfrecu[mfrecu > 0]
    }
    else {
        datau = unique(data[, 1:(p - 1)])
        rowsu = dim(datau)[1]
        tempo = row.matches(datau[1, 1:(p - 1)], data[, 1:(p - 
            1)])
        ind = tempo
        mfrecu = max(table(data[ind, p]))
        cardi = length(ind)
        datatempo = as.matrix(data[-ind, ])
        for (j in 2:(rowsu - 1)) {
            tempov = datau[j, 1:p - 1]
            tempo = row.matches(tempov, datatempo[, 1:p - 1])
            ind = tempo
            mfrecu = c(mfrecu, max(table(datatempo[ind, p])))
            datatempo = datatempo[-ind, ]
            cardi = c(cardi, length(ind))
        }
        if (length(dim(datatempo)) == 0) {
            cardi = c(cardi, dim(t(as.matrix(datatempo)))[1])
            mfrecu = c(mfrecu, 1)
        }
        else {
            cardi = c(cardi, dim(datatempo)[1])
            mfrecu = c(mfrecu, max(table(datatempo[, p])))
        }
    }
    incon = sum(cardi - mfrecu)/n
    incon
}
