univariate0 <-function (X) 
{
    p <- ncol(X)
    noms <- colnames(X)
    dev.new()
    if (p<10) {
par(mfrow = n2mfrow(p))}
else{
par(mfrow=c(3,3))}

    for (j in 1:p) {
                    if (is.ordered(X[, j])) {
                barplot(table(X[, j]), main = noms[j], ylab="Effectif",las = 1,col=rev(brewer.pal(3,name="PuRd"))[1] )
            }
            else {
                if (is.factor(X[, j])) {
if(length(table(X[,j]))<6){pie(rev(sort(table(X[,j]))),main = noms[j],col=rev(brewer.pal(3,name="PuRd"))[c(1,2,3,2,3)])}
else{
                  www <- as.numeric(table(X[, j]))
                  names(www) <- names(table(X[, j]))
                  dotchart(sort(www), xlab="Effectif",main = noms[j],pch=16,color=rev(brewer.pal(3,name="PuRd"))[1])
}
                }
                else {
                  if (is.numeric(X[, j])) {
			if(length(table(X[, j]))>11){
                    hist(X[, j], main = noms[j],ylab="Effectif", xlab = "", 
                      las = 1,col=rev(brewer.pal(3,name="PuRd"))[1])}
else{
plot(as.numeric(names(table(X[, j]))),table(X[, j]),ylim=c(0,max(table(X[, j]))),main = noms[j],ylab="Effectif",xlab="",col=rev(brewer.pal(3,name="PuRd"))[1],type="h",lwd=3)
axis(2)

}
                  }
                }
            
        }
   
if(((j%%9)==0)&(p!=9)){
dev.new()
par(mfrow=c(3,3))
}
 }
}