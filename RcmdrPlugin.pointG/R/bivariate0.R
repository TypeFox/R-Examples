bivariate0<-function (X,y,...) 
{
noms <- colnames(X)
p <- ncol(X)
    
# Dimensionnement des fenêtres
dev.new()
if(p<3){
par(mfrow = n2mfrow(p))}
else{
par(mfrow = c(2,2))}

   
for (j in 1:p) {
        if (is.factor(X[, j]) & is.factor(y)) {
            mosaicplot(table(y, X[, j]), main = noms[j],
                las = 1,col=rev(brewer.pal(3,name="PuRd")),...)
        }
        if (is.numeric(X[, j]) & is.numeric(y)) {
            plot(y, X[, j], main = noms[j], ylab = "", 
                las = 1,col=rev(brewer.pal(3,name="PuRd"))[1],...)
        }
        if (is.numeric(X[, j]) & is.factor(y)) {
            boxplot(X[, j] ~ y, main = noms[j], ylab = "",col=rev(brewer.pal(3,name="PuRd"))[1],...)
        }
        if (is.numeric(y) & is.factor(X[, j])) {
            boxplot(y ~ X[, j], main = noms[j], ylab = "", 
                horizontal = TRUE, las = 1,col=rev(brewer.pal(3,name="PuRd"))[1],...)
        }

if(((j%%4)==0)&(p!=4)){
dev.new()
par(mfrow=c(2,2))
}



    }
}
