plot.arx <-
function(x, spec=NULL, col=c("red","blue"),
  lty=c("solid","solid"), lwd=c(1,1), ...)
{
  ##one is always non-NULL?:
  if(!is.null(x$mean.fit) || !is.null(x$var.fit)){

    ##lwd argument:
    if(length(lwd)==1){
      print("lwd needs two arguments, but only one provided. Single argument applied to all lines plotted.")
      lwd=rep(lwd,2)
    }else if (length(lwd)>2){
      print("lwd needs two arguments, but more provided. First two used.")
      lwd=lwd[1:2]
    }

    ##lty argument:
    if(length(lty)==1){
      print("lty needs two arguments, but only one provided. Single argument applied to all lines plotted.")
      lty=rep(lty,2)
    }else if (length(lwd)>2){
      print("lty needs two arguments, but more provided. First two used.")
      lty=lty[1:2]
    }

    ##col argument:
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
    } #if length(col!=2)

    ##spec argument:
    if(is.null(spec)){
      if(!is.null(x$mean.results)){
        spec <- "mean"
      }
      if(is.null(x$mean.results)
         && !is.null(x$variance.results) ){
        spec <- "variance"
      }
      if(!is.null(x$mean.results)
         && !is.null(x$variance.results) ){
        spec <- "both"
      }
    }else{
      spec.type <- c("mean", "variance", "both")
      which.type <- charmatch(spec, spec.type)
      spec <- spec.type[which.type]
    }

    ##plot if spec is not NULL:
    if(!is.null(spec)){

      ##if variance modelled, plot square root of fitted variance and absolute residuals against time
      if(spec=="variance" || spec=="both"){
        vfitted <- sqrt(x$var.fit)
        vactual <- abs(x$resids)
      }

      ##if mean modelled, plot fitted and actual values against time
      if(spec=="mean" || spec=="both"){
        mfitted <- x$mean.fit
        mactual <- zoo(x$aux$y, order.by=x$aux$y.index)
      }
      actual.name <- x$aux$y.name
      residsStd <- x$resids.std

      ##get current par-values:
      def.par <- par(no.readonly=TRUE)

      ##set new par values for plot
      if(spec=="both") {##if both mean and variance modelled, plot both
        par(mfrow=c(3,1))
      }else {##else just plot the one specified
        par(mfrow=c(2,1))
      }

      ##what is happening here (add comment)?:
      par(mar=c(2,2,0.5,0.5))
      if(spec=="mean" || spec=="both") {##plotting mean variables

        ##check whether ?? zoo object is regular, then plot:
        if(is.regular(mactual)) {
          plot(mactual, main = "",
             ylim=range(min(mactual,mfitted),max(mactual,mfitted)),
             type="l",ylab="",xlab="",col=col[2])
        } else {##if irregular, plot manually
          plot(as.Date(index(mactual)),coredata(mactual), main = "",
             ylim=range(min(mactual,mfitted),max(mactual,mfitted)),
             type="l",ylab="",xlab="",col=col[2])
        }

        ##check whether ?? zoo object is regular, then plot:
        if(is.regular(mfitted)) {
          lines(mfitted,col=col[1])
        } else {
          lines(as.Date(index(mfitted)),coredata(mfitted),col=col[1])
        }
        legend("topleft",lty=lty,lwd=lwd,ncol=2,col=col[c(2,1)],legend=c(actual.name,"fitted"),bty="n")

      } #close ##comment? what is happening here?

      ##plotting variance parts:
      if(spec=="variance" || spec=="both") {

        ##add comment?
        if(is.regular(vactual)) {
          plot(vactual, main = "",
             ylim=range(min(vactual,vfitted,na.rm=TRUE),max(vactual,vfitted,na.rm=TRUE)),
             type="l",ylab="",xlab="",col=col[2])
        } else {
          plot(as.Date(index(vactual)),coredata(vactual), main = "",
             ylim=range(min(vactual,vfitted,na.rm=TRUE),max(vactual,vfitted,na.rm=TRUE)),
             type="l",ylab="",xlab="",col=col[2])
        }

        ##add comment?
        if(is.regular(vfitted)) {
          lines(vfitted,col=col[1])
        } else {
          lines(as.Date(index(vfitted)),coredata(vfitted),col=col[1])
        }
        legend("topleft",lty=lty,lwd=lwd,ncol=2,col=col[c(2,1)],
          legend=c("abs(residuals)","fitted sd"),bty="n")

      } #close plotting variance parts

      ##add comment?
      if(is.regular(residsStd)) {
        plot(residsStd,type="h",col=col[1])
      } else {
        plot(as.Date(index(residsStd)),coredata(residsStd),type="h",col=col[1])
      }
      abline(0,0)
      legend("topleft",lty=1,col=col[1],legend=c("standardised residuals"),bty="n")

      #return to old par-values:
      par(def.par)

    } #close if(!is.null(spec))

  } #close if(!is.null(x$mean.fit))

}
