plot.isat <-
function(x, col=c("red","blue"),
  lty=c("solid","solid"), lwd=c(1,1), coef.path=TRUE, ...)
{

#  ##check if mean quation:
#  if( is.null(x$mean.results) ){
#    cat("No mean equation to plot\n")
#  }

  ##if fitted mean:
  if(!is.null(x$mean.fit)){

    ##check line width:
    if(length(lwd)==1){
      print("lwd needs two arguments, but only one provided. Single argument applied to all lines plotted.")
      lwd=rep(lwd,2)
    }else if (length(lwd)>2){
      print("lwd needs two arguments, but more provided. First two used.")
      lwd=lwd[1:2]
    }

    ##check line type:
    if(length(lty)==1){
      print("lty needs two arguments, but only one provided. Single argument applied to all lines plotted.")
      lty=rep(lty,2)
    }else if (length(lwd)>2){
      print("lty needs two arguments, but more provided. First two used.")
      lty=lty[1:2]
    }

    ##check colour:
    if(length(col)!=2){

      #####randomcol - returns random combination of colours of length 2
      randomcol <- function()
      {
        r.r <- runif(2)
        while(round(r.r[1],1)==round(r.r[2],1)) {#don't want colours too similar so add check
          r.r <- runif(2)
        }
        g.r <- runif(2)
        while(round(g.r[1],1)==round(g.r[2],1)) {#don't want colours too similar so add check
          g.r <- runif(2)
        }
        b.r <- runif(2)
        while(round(b.r[1],1)==round(b.r[2],1)) {#don't want colours too similar so add check
          b.r <- runif(2)
        }
        col<-rgb(runif(2),runif(2),runif(2))
        return(col)
      } #end randomcol

      ###clashcol function - returns clashing (opposite) combination of colours
      clashcol <- function()
      {
        #using http://forum.processing.org/one/topic/the-opposite-of-a-color.html formula for opposite colour
        r.1 <- runif(1)
        g.1 <- runif(1)
        b.1 <- runif(1)
        b.2 <- min(r.1,min(g.1,b.1)) + max(r.1,max(g.1,b.1))
        col <- rgb(c(r.1,b.2-r.1),c(g.1,b.2-g.1),c(b.1,b.2-b.1))
        return(col)
      } #end clashcol

      ##if random:
      if(col[1]=="random") {
        col <- randomcol()
      }else if(col[1]=="awful.clash") {
        col <- clashcol()
      }else{
        print("Wrong number of colours specified; using random set of colours instead.")
        col<-randomcol()
      }
    }

    ##get fitted and actual values, and dependent variable name
    fitted <- x$mean.fit
    actual <- zoo(x$aux$y, order.by=x$aux$y.index)
    residuals <- x$resids.std
    actual.name <- x$aux$y.name

    ##get current par-values, set new ones:
    def.par <- par(no.readonly=TRUE)
    par(mar=c(2,2,0.5,0.5))

    ##if isat or coef-path:
    if( (x$gets.type=="isat" || coef.path==TRUE) && length(x$ISnames)!=0 ){
      par(mfrow=c(3,1))
      is.x <- cbind(x$aux$mX[,x$aux$mXnames %in% x$ISnames])
      is.coef.ests <- coef.isat(x)[x$ISnames]
      coef.path.0 <- zoo(is.x %*% is.coef.ests, order.by=x$aux$y.index)
    }else{
      par(mfrow=c(2,1))
    }

    ##comment?
    if(is.regular(actual)) {
      plot(actual, main = "",ylim=range(min(actual,fitted,na.rm=TRUE),max(actual,fitted,na.rm=TRUE)),
           type="l",ylab="",xlab="",col=col[2])
    }else{
      plot(as.Date(index(actual)),coredata(actual), main = "",ylim=range(min(actual,fitted,na.rm=TRUE),max(actual,fitted,na.rm=TRUE)),
           type="l",ylab="",xlab="",col=col[2])
    }
    if(is.regular(fitted)) {
      lines(fitted,col=col[1])
    }else{
      lines(as.Date(index(fitted)),coredata(fitted),col=col[1])
    }

    legend("topleft",lty=lty,lwd=lwd,ncol=2,col=col[c(2,1)],legend=c(actual.name,"fitted"),bty="n")
    if(is.regular(residuals)) {
      plot(residuals,type="h",col=col[1])
    }
    else{
      plot(as.Date(index(residuals)),coredata(residuals),type="h",col=col[1])
    }

    abline(0,0)
    legend("topleft",lty=1,col=col[1],legend="standardised residuals",bty="n")
#    legend("topleft",lty=1,col=col[1],legend=c(paste(actual.name,"standardised residuals",sep=": ")),bty="n")

    ##coefficient path
    if( (x$gets.type=="isat" | coef.path==TRUE) & length(x$ISnames)!=0 ) {
      ## we only get standard error bars if TIS *not* run


      ###if tis is there and it is not null, then don't plot, else


      if(!is.null(as.list(x$call)$tis) && as.list(x$call)$tis==TRUE){
        cat("\nNB: Because TIS selected, coefficient standard errors invalid hence not plotted\n", sep="")
        ylim.values <- range(coef.path.0)
        if(is.regular(coef.path.0)) {
          ylim.values <- range(coef.path.0)
          plot(coef.path.0,type="l",col=col[1],ylim=ylim.values)
        }else{
          plot(as.Date(index(coef.path.0)),coredata(coef.path.0),type="l",col=col[1],ylim=ylim.values)
        }
      }   else {

        coef.path.v <- isatvar(x)
        if(is.regular(coef.path.0)) {
          ylim.values <- range(min(coef.path.0-qt(0.975, NROW(coef.path.0))*coef.path.v$const.se),
                               max(coef.path.0+qt(0.975, NROW(coef.path.0))*coef.path.v$const.se))
          plot(coef.path.0,type="l",col=col[1],
               ylim=ylim.values)
          lines(coef.path.0+qt(0.975, NROW(coef.path.0))*coef.path.v$const.se,type="l",col=col[1],lty=3)
          lines(coef.path.0-qt(0.975, NROW(coef.path.0))*coef.path.v$const.se,type="l",col=col[1],lty=3)
        }else{
          plot(as.Date(index(coef.path.0)),coredata(coef.path.0),type="l",col=col[1],ylim=ylim.values)
          lines(as.Date(index(coef.path.0)),coredata(coef.path.0)+qt(0.975, NROW(coef.path.0))*coef.path.v$const.se,type="l",col=col[1],lty=3)
          lines(as.Date(index(coef.path.0)),coredata(coef.path.0)-qt(0.975, NROW(coef.path.0))*coef.path.v$const.se,type="l",col=col[1],lty=3)

        }

      }

      abline(0,0,lty=3)
      legend("topleft",lty=1,col=col[1],legend=c(paste(actual.name,"Coefficient Path",sep=": ")),bty="n")
    }

    #return to old par-values:
    par(def.par)
  }
}
