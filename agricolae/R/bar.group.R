`bar.group` <-
function(x,horiz=FALSE, ...) {
y<-x[,2]
names(y)<-x[,1]
nivel<-x[,3]
n<-length(y)
indice<-barplot(y,horiz=horiz, ...)
tope<-max(y)/10
for ( i in 1:n) {
if (horiz) text(y[i]+tope,indice[i],nivel[i])
else text(indice[i],y[i]+tope,nivel[i])
}
invisible(list(index=indice,means=y))
}

