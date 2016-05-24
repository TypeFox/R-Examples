plot.tfa <-
function(x,xvar=NULL,yvar=NULL,...){
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
    yv<-x[[length(x)]]
    yvar<-names(x)[length(x)]
}
else{
    if(yvar=="p") yv<-x$p
    if(yvar=="lambda") yv<-x$lambda
    if(yvar=="inertia") yv<-x$inertia
}
plot(xv,yv,type="l",xlab=xvar,ylab=yvar,...)
}
