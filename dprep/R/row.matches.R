row.matches <-
structure(function (y, X) 
{
    i <- seq(nrow(X))
    j <- 0
    while (length(i) && (j <- j + 1) <= ncol(X)) i <- i[X[i, 
        j] == y[j]]
    i
}, source = c("function(y, X) {", "#********************************************", 
"#This functions finds out the rows in matrix X", "#that are equal to the vector y", 
"#y: a vector", "#X: a matrix", "#**************************************************", 
"        i <- seq(nrow(X)) ", "        j <- 0 ", "        while(length(i) && (j <- j + 1) <= ncol(X)) ", 
"                i <- i[X[i, j] == y[j]] ", "        i ", "}"
))
