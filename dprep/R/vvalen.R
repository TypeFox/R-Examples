vvalen <-
function (data) 
{
    p = dim(data)[2]
    data[,p]=as.numeric(factor(data[,p]))
    clases = length(table(data[, p]))
    testlist = vvalen1(data, 1)
    testlist = list(testlist)
    for (i in 2:clases) {
        a = vvalen1(data, i)
        testlist = c(testlist, list(a))
    }
    z = kruskal.test(testlist)
    return(z)
}
