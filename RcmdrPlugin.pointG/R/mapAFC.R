mapAFC<-
function(X){
Y<-data.frame(matrix(X, nrow = nrow(X)))
rownames(Y)<-rownames(X)
colnames(Y)<-colnames(X)
scatter(dudi.coa(Y,scannf = FALSE, nf = 2))
}
