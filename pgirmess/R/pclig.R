"pclig" <-
function(matr){
    mat<-na.omit(matr)
    mat/rowSums(mat)
}

