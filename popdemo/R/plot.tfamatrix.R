plot.tfamatrix <-
function(x,xvar=NULL,yvar=NULL,mar=c(1.1,1.1,0.1,0.1),...){
order<-dim(x$p)[1]
lay<-rbind(rep(max(x$layout)+1,order),x$layout)
layout(lay)
par(mar=mar)
if(is.null(xvar)){
    xv<-x$p
    xvar<-"p"
}
else{
    if(xvar=="p") xv<-x$p
    if(xvar=="lambda") xv<-x$lambda
    if(xvar=="inertia") xv<-x$inertia
}
if(is.null(yvar)){
    yv<-x[[length(x)-1]]
    yvar<-names(x)[length(x)-1]
}
else{
    if(yvar=="p") yv<-x$p
    if(yvar=="lambda") yv<-x$lambda
    if(yvar=="inertia") yv<-x$inertia
}
elementrow<-x$rows
elementcol<-x$cols
for(i in 1:order){
    for(j in 1:order){
        if(x$layout[i,j]!=0){
            plot(xv[i,j,],yv[i,j,],type="l",xaxt="n",yaxt="n",...)
            axis(side=1,tck=0.05,padj=-1.7,cex.axis=0.9)
            axis(side=2,tck=0.05,padj=1.4,cex.axis=0.9)
        }
    }
}
par(mar=c(0,0,0,0))
plot(0,type="n",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
text(1,0,paste(yvar,"~",xvar),cex=2)
}
