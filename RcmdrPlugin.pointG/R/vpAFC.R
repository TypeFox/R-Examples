vpAFC<-
function (X) 
{

    result <- dudi.coa(data.frame(matrix(X, nrow = nrow(X))), 
        scannf = FALSE, nf = 1)$eig[1]
    result
}
