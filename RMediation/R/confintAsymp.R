confintAsymp <- function(mu, Sigma, quant=NULL, alpha=0.05, type="asymp", plot=FALSE, plotCI=FALSE){
  #Sigma <- vech.reverse(Sigma)
  q1 <- quant
  fx <- deriv(quant,names(mu),func=TRUE)
  muList <- as.list(mu)
  grad1 <- do.call("fx",muList)
  grad <- as.vector(attr(grad1,"gradient"))
  quantSE <- sqrt(t(grad)%*%Sigma%*%grad)
  quant <- parse(text=sub("~","",quant))
  quantMean <- eval(quant,muList)
  q.a <- qnorm(alpha/2,0,1)
  CI <- c(quantMean+q.a*quantSE,quantMean-q.a*quantSE)
  res <- list(CI,quantMean,quantSE)
  names(res) <- c(paste((1-alpha/2)*100,"% ", "CI",sep=""), "Estimate", "SE")
  attr(res,"quant")  <- q1
  outer <- FALSE #outer position
  mcex <- .8
  
  if (plot & quantSE>40*.Machine$double.eps)
  {
    if(type=="all"){     
      f <- Vectorize(function(x) dnorm(x,quantMean,quantSE), "x")
      curve(f,lwd=2, col="blue", lty=2, add=TRUE)
      if(plotCI){
        max1 <- par("usr")[2]
        min1 <- par("usr")[1]
        yci<-par("usr")[3]-.1*diff(par("usr")[3:4])/25
        arrows(CI[1],yci,CI[2],yci,length=0,angle=90,code=3,cex=1.5,lwd=2, col="blue",lty=2)
        points(quantMean,yci,pch=19,cex=1.5, col="blue")
        #New
        mtext(paste("LL=",round(CI[1],3)),side=3,line=1,outer=outer,at=max1-1*(max1-min1)/9,col="blue", cex=mcex, adj=0)
        mtext(paste("UL=",round(CI[2],3)),side=3,line=2,outer=outer,at=max1-1*(max1-min1)/9, col="blue" , cex=mcex, adj=0)
        mtext("Asymptotic Normal",side=3,line=3,outer=outer, at=max1-1*(max1-min1)/9, col="blue", cex=mcex, adj=0)
      }                    
    }
    
    if(type %in% c("asymp", "Asymp")){
      f <- Vectorize(function(x) dnorm(x,quantMean,quantSE),"x")
    curve(f,from=quantMean-2.58*quantSE,to=quantMean+2.58*quantSE, lwd=2, xlab="", ylab="", axes=FALSE)
    range1 <- c(quantMean-2.58*quantSE, quantMean+2.58*quantSE)
    max1 <- range1[2]
    min1 <- range1[1]
    xrange<-seq(min1,max1,length=7)
    xrange<-pretty(xrange,n=9)
    mtext(quant,1,4.5)
    axis(1,xrange,xrange, line=2.5);
    axis(2,line=1.1)
    mtext(paste("LL=",round(CI[1],3)),side=3,line=1,outer=outer,at=max1-4*(max1-min1)/9,col="black", cex=mcex)
    mtext(paste("UL=",round(CI[2],3)),side=3,line=2,outer=outer,at=max1-4*(max1-min1)/9, col="black", cex=mcex)
    mtext("Asymptotic",side=3,line=3,outer=outer, at=max1-4*(max1-min1)/9, col="black", cex=mcex)

    if(plotCI){
      yci<-par("usr")[3]-1.2*diff(par("usr")[3:4])/25
      arrows(CI[1],yci,CI[2],yci,length=0,angle=90,code=3,cex=1.5,lwd=2)
      points(quantMean,yci,pch=19,cex=1.5)
    } #if
  } 
  } #if
  
  if (plot & quantSE<=40*.Machine$double.eps){
    if(type %in% c("asymp", "Asymp")){
#     range1 <- c(quantMean-2.58*quantSE, quantMean+2.58*quantSE)
#     max1 <- range1[2]
#     min1 <- range1[1]
#     xrange<-round(seq(min1,max1,length=7),1)
#     axis(1,c(xrange, round(c(quantMean,CI),2)));axis(2)
#     smidge <- par("cin")*abs(par("tcl"))
#     text(max1-0*(max1-min1)/9,(par("usr")[4]),adj=0, bquote(mu== .(round(quantMean,3)) ))
#     text(max1-0*(max1-min1)/9,(par("usr")[4]-1*par("cxy")[2]),adj=0, bquote(sigma== .(round(quantSE,3)) ))
    cat("
        \n ******************************************* 
Note: Because the SE =",quantSE,", the plot of the sampling distribution using aymptotic normal is not available.
        \n *******************************************\n")
    
    }
    
    if(type=="all"){
      max1 <- par("usr")[2]
      min1 <- par("usr")[1]
      yci<-par("usr")[3]+2.5*diff(par("usr")[3:4])/25
      #arrows(CI[1],yci,CI[2],yci,length=0,angle=90,code=3,cex=1.5,lwd=2, col="blue",lty=2)
      #points(quantMean,yci,pch=19,cex=1.5, col="blue")
#       text(max1+0*(max1-min1)/9,(par("usr")[4]-2.5*par("cxy")[2]),adj=0, paste("LL=",round(CI[1],3)), col="blue")
#       text(max1+0*(max1-min1)/9,(par("usr")[4]-3.5*par("cxy")[2]),adj=0, paste("UL=",round(CI[2],3)), col="blue")
      mtext(paste("LL=",round(CI[1],3)),side=3,line=2,outer=outer,at=max1-1*(max1-min1)/9,col="blue", cex=mcex)
      mtext(paste("UL=",round(CI[2],3)),side=3,line=3,outer=outer,at=max1-1*(max1-min1)/9, col="blue", cex=mcex)
      mtext("Asymptotic",side=3,line=4,outer=outer, at=max1-1*(max1-min1)/9, col="blue", cex=mcex)
      cat("Because the SE,",quantSE,", is near zero, the plot of the sampling distribution using aymptotic normal is not available.\n")
      
    }
    
    }
  
  return(res)
}
