`bar.err` <-
function(x,variation=c("SE","SD","range"),horiz=FALSE, bar=TRUE,...) {
variation<-match.arg(variation)
y<-x[,1]
names(y)<-rownames(x)
if( variation=="SE" ) {
names(x)-> nombre
if(sum(nombre=="SE")==0) std.err<-x$"std"/sqrt(x$"r")
else std.err<-x$"SE"
nivel0<-y-std.err
nivel1<-y+std.err
}
if( variation=="SD" ) {
nivel0<-y-x$"std"
nivel1<-y+x$"std"
}
if( variation=="range" ) {
nivel0<-x$"Min"
nivel1<-x$"Max"
}
n<-length(y)
if (bar) {
indice<-barplot(y,horiz=horiz, ...)
tope<-max(nivel1)/20
}
else {
indice<-barplot(y,horiz=horiz, border=0, ...)
if(horiz)lines(y,indice,type="b")
else lines(indice,y,type="b")
}
for ( i in 1:n) {
if (horiz)  {
lines(rbind(c(nivel0[i],indice[i]),c(nivel1[i],indice[i])))
text( cex=1,nivel0[i],indice[i],"[")
text( cex=1,nivel1[i],indice[i],"]")
}
else {
lines(rbind(c(indice[i],nivel0[i]),c(indice[i],nivel1[i])))
text( cex=1,indice[i],nivel0[i],"---")
text( cex=1,indice[i],nivel1[i],"---")
}
}
invisible(list(index=indice,means=y))
}

