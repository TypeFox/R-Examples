`plot.modTempEff` <-function(x, which=c("cold","heat"), add=FALSE, new=TRUE, var.bayes=FALSE, delta.rr=TRUE, level=0.95,
    updown=TRUE, col.shade=NULL, leg=NULL, ...){
    new.args<-list(...)
    if(length(x$ToTheat)<=0){
        if(is.null(new.args$xlab)) new.args$xlab<-"Time"
        if(is.null(new.args$ylab)) new.args$ylab<-"Fitted Values"
        if(length(x$fit.seas)<=0) {stop("the model does not include csdl() or seas() terms")
          } else {
          if(add) lines(exp(x$fit.seas),...) else plot(exp(x$fit.seas), 
            xlab=new.args$xlab, ylab=new.args$ylab, type="l", ...)
          return(invisible(x))
        }
      }

    which<-match.arg(which, several.ok = TRUE)
    cold<-match("cold",which,nomatch = 0)>0
    heat<-match("heat",which,nomatch = 0)>0
    if(add) {
      new<-FALSE
      if((cold+heat)>=2) stop("'add=TRUE' works only with a *single* DL curve")
      }
    if(new) dev.new() 
    if(updown) par(mfrow=c(cold+heat, 1)) else par(mfrow=c(1,cold+heat))
    z<-qnorm((1+level)/2)
    xf<- -x$betaCold
    xc<-x$betaHeat
    if(var.bayes) {
        SE.f<-x$SE.c.bayes
        SE.c<-x$SE.h.bayes
        } else {
          SE.f<-x$SE.c
          SE.c<-x$SE.h
               }
    if(cold){
      etich<-"logRR for Cold"
      rrr<- cbind(xf,xf-z*SE.f,xf+z*SE.f)
      if(!add){
          if(delta.rr) {
              etich<-"% change in mortality for Cold"
              rrr<-100*(exp(rrr)-1)
              }
        if(is.null(new.args$ylab)) new.args$ylab<-etich
        if(is.null(new.args$xlab)) new.args$xlab<-"Lag (day)"
        if(is.null(new.args$lwd)) new.args$lwd<-1
          matplot(0:(length(xf)-1),rrr,
            type="l",lty=c(1,2,2),col=1,ylab=new.args$ylab, xlab=new.args$xlab, lwd=new.args$lwd)
          if(!is.null(col.shade)){
            x.shade<-c(0:(length(xf)-1), rev(0:(length(xf)-1))) #xx<-c(x, rev(x))
            y.shade<-c(rrr[,2],rev(rrr[,3]))
            #matplot(x, cbind(y1,y2))
            polygon(x.shade,y.shade,col=col.shade,border=NA)
            lines(0:(length(xf)-1),rrr[,1],...)
          }
          abline(h=0,lty=3,col=1)
          if(!is.null(leg)) legend("topright",leg[1],bty="n")
          } else {
        matlines(0:(length(xf)-1),rrr,col=1,...)
        }
      }
    if(heat){
      etich<-"logRR for Heat"
      rrr<-cbind(xc,xc-z*SE.c,xc+z*SE.c)
      if(!add){
        if(delta.rr) {
            rrr<-100*(exp(rrr)-1)
            etich<-"% change in mortality for Heat"
            }
        if(is.null(new.args$ylab)) new.args$ylab<-etich
        if(is.null(new.args$xlab)) new.args$xlab<-"Lag (day)"
        if(is.null(new.args$lwd)) new.args$lwd<-1
        matplot(0:(length(xc)-1),rrr,type="l",
          lty=c(1,2,2),col=1,ylab=new.args$ylab, xlab= new.args$xlab, lwd=new.args$lwd )
        if(!is.null(col.shade)){
            x.shade<-c(0:(length(xc)-1), rev(0:(length(xc)-1))) #xx<-c(x, rev(x))
            y.shade<-c(rrr[,2],rev(rrr[,3]))
            #matplot(x, cbind(y1,y2))
            polygon(x.shade,y.shade,col=col.shade,border=NA)
            lines(0:(length(xc)-1),rrr[,1],...)
          }
        abline(h=0,lty=3,col=1)
        if(!is.null(leg)) legend("topright",leg[2],bty="n")
        } else {
          matlines(0:(length(xc)-1),rrr,col=1,...)}
      }
    } #end_funct

