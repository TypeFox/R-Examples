plot.AssetPricing <- function(x,witch=c("price","expVal","vdot"),xlim=NULL,
                              ylim=NULL,lty=NULL,cols=NULL,xlab=NULL,ylab=NULL,
                              main=NULL,main.panel=NULL,groups=NULL,add=FALSE,
                              gloss=FALSE,glind=NULL,extend=0.3,col.gloss=1,
                              cex.gloss=0.8,mfrow=NULL,...) {
#
witch <- match.arg(witch)
xxx   <- switch(EXPR=witch,price=x[["x"]],
                           expVal=x[["v"]],
                           vdot=x[["vdot"]])
plot(xxx,xlim=xlim,ylim=ylim,lty=lty,cols=cols,xlab=xlab,ylab=ylab,
     main=main,main.panel=main.panel,groups=groups,add=add,
     gloss=gloss,glind=glind,extend=extend,col.gloss=col.gloss,
     cex.gloss=cex.gloss,mfrow=mfrow,...)
}
