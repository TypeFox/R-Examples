confintMC <- function(mu, Sigma, quant=NULL, alpha=0.05, type="MC", plot=FALSE, plotCI=FALSE, n.mc = 1e+06, H0=FALSE, mu0, Sigma0, ...){
  q1 <- quant
  quant <- parse(text=sub("~","",quant))
  df <- data.frame(mvrnorm(n.mc,mu,Sigma))
  colnames(df) <-names(mu)
  quant.vec <- eval(quant,df) # MC Vector
  CI <- quantile(quant.vec,c(alpha/2,1-alpha/2))
  names(CI) <- c(paste((alpha/2*100),"%"),paste((1-alpha/2)*100,"%"))
  quantMean <- mean(quant.vec)
  quantSE <- sd(quant.vec)
  quantError <- quantSE/n.mc
  pMinH1 <- mean(quant.vec > 0) 
  pMinH1 <- 2*(min(pMinH1, 1-pMinH1))
  # This is to determine min and max of MC Samples--Do Not Delete this
 # if(H0) xrange <- c(min(-4*quantSE,quantMean-4*quantSE), max(4*quantSE ,quantMean+4*quantSE) )
 # else xrange <- c(quantMean-4*quantSE, quantMean+4*quantSE)
  xrange <- c(quantMean-4*quantSE, quantMean+4*quantSE)
  max1 <- xrange[2]
  min1 <- xrange[1]
  

  # Added for min null test 6-14-14
  # Calculates Acceptance Region
  if(H0)
  {
    if(missing(Sigma0)| is.null(Sigma0) ) Sigma0 <- Sigma
    if(!is.matrix(Sigma0)){
      if(length(mu)!= (sqrt(1 + 8 * length(Sigma0)) - 1)/2) stop(
        paste("Please check the length of", sQuote("Sigma0"),"and",sQuote("mu"),". If the length(dimension) of the", sQuote("mu"),"vector (",length(mu),") is correct, the stacked lower triangle matrix", sQuote("Sigma0"), "must have ",((2*length(mu)+1)^2-1)/8, "elements, instead of", length(Sigma0)) 
      )
      
      Sigma0 <- lav_matrix_vech_reverse(Sigma0) #converts to a symmetric matrix
    }
    
    #If mu0 is not specified, we use conservative min approach
    if(is.null(mu0) | missing(mu0) ){
      mu0 <- mu
      mu0s <- mu0/sqrt(diag(Sigma0)); #Srandardized
      mu0[which(mu0s==min(mu0s))] <- 0 # setting the smallest z value to 0
    }
    names(mu0) <- names(mu)
    df0 <- data.frame(mvrnorm(n.mc,mu0,Sigma0))
    colnames(df0) <-names(mu0)
    H0quant.vec <- eval(quant,df0)
    H0Mean <- mean(H0quant.vec)
    H0SE <- sd(H0quant.vec)
    H0CI <- quantile(H0quant.vec, c(alpha/2,1-alpha/2) )
    names(H0CI) <- c(paste((alpha/2*100),"%"),paste((1-alpha/2)*100,"%"))
    #yciH0<-par("usr")[3] 
    pMinH0 <- mean(H0quant.vec > quantMean) 
    pMinH0 <- 2*(min(pMinH0, 1-pMinH0))
    H0xy <- H0quant.vec + quantMean
    H0xy <- H0xy[H0xy>min1 & H0xy<max1] #This is used to make the plot prettier
#    H0xy <- H0xy+quantMean # We add mean to make it comparable with the CI
    H0xyDens <- density(H0xy)
    H0CI <- H0CI + quantMean # We add mean to make it comparable with the CI
    H0res <- list(CI=H0CI, Estimate= H0Mean,SE= H0SE,p= pMinH0, mu0) #Results
    #names(H0res) <- c( paste( (1-alpha)*100, "% ", "AC",sep="" ) , "Mean", "SE", "p", "mu0")   
    #attr(H0res,"quant")  <- q1
  }

  ###########   Plot ##########################
  
  
  if (plot){
    
    outer <- FALSE #outer position
    mcex <- .8

    if(type=="all")
    {
      res.asymp <- confintAsymp(mu=mu, Sigma=Sigma, quant=quant, type=type, alpha=alpha)
      range.asymp <- c(res.asymp$Estimate-5*res.asymp$SE, res.asymp$Estimate+5*res.asymp$SE)
      max1 <- max(max1,range.asymp)
      min1 <- min(min1,range.asymp)      
    }
    xy <- quant.vec[quant.vec>min1 & quant.vec<max1]
    xyDens <- density(xy)
    #smidge <- par("cin")*abs(par("tcl"))
    # Added for min null test 6-14-14
    # To get a more reasonable y range
    if(H0){ 
      yrange <- range(quantile(xyDens$y,c(alpha/10,1-alpha/10)),quantile(H0xyDens$y,c(alpha/10,1-alpha/10)))
    }
    else{
      yrange <- quantile(xyDens$y,c(alpha/10,1-alpha/10))
    }
    
    plot(xyDens,xlab="", ylab="", axes=FALSE, xlim=xrange, ylim=yrange, main="", lwd=2)
    mtext(quant,1,5)
    xrange<-pretty(xrange,n=9)
    axis(1,xrange,xrange, line=2.5);
    axis(2,line=1.1)
    
    # This adds legends, bars etc
    if(H0)
    {
      lines(H0xyDens, lty=2, col="blue", lwd=2) #Reference Dist
      #arrows(H0CI[1],yciH0,H0CI[2],yciH0,length=0,angle=90,code=3,cex=1.5,lwd=5,lty=1, col="blue")
      segments(H0CI[1],par("usr")[3],H0CI[1], par("usr")[4]/2, cex=1.5,lwd=2,lty=2, col="blue")
      segments(H0CI[2],par("usr")[3],H0CI[2], par("usr")[4]/2, cex=1.5,lwd=2,lty=2, col="blue")
      
      mtext(paste("P Value=",round(pMinH0,4) ), side=3, line=-1 , outer=outer, at=max1-1*(max1-min1)/9, cex=mcex, adj=0, col="blue")
      
      mtext(paste("Kurtosis=",round(kurtosis(H0quant.vec, type=2),3)), side=3, line=0, outer=outer, at=max1-1*(max1-min1)/9, cex=mcex, adj=0, col="blue")
      mtext(paste("Skewness=",round(skewness(H0quant.vec, type=2),3)), side=3, line=1, outer=outer, at=max1-1*(max1-min1)/9, cex=mcex, adj=0, col="blue")
      mtext(paste("Critical Value=", round(H0CI[1],3)), side=3,line=2, outer=outer, at=max1-1*(max1-min1)/9, cex=mcex, adj=0, col="blue")
      mtext(paste("Critical Value=", round(H0CI[2],3)), side=3,line=3, outer=outer, at=max1-1*(max1-min1)/9, cex=mcex, adj=0, col="blue")
      mtext(paste("H0:",quant,"=0"), side=3,line=4, outer=outer, at=max1-1*(max1-min1)/9, cex=mcex, adj=0, col="blue",font=2)
      text(H0CI[1], par("usr")[4]/2, labels="Lower Critical Value",adj=c(.5,0), cex=mcex )
      text(H0CI[2], par("usr")[4]/2, labels="Upper Critical Value", adj=c(.5,0), cex=mcex)
    }
    
    
    if(type %in% c("mc", "MC")){
      #New- 1/24/14-DT
      mtext(paste("P Value=",round(pMinH1,4)), side=3, line=-1, outer=outer, at=max1-3*(max1-min1)/9, cex=mcex,adj=0)
      mtext(paste("Kurtosis=",round(kurtosis(quant.vec, type=2),3)), side=3, line=0, outer=outer, at=max1-3*(max1-min1)/9, cex=mcex, adj=0)
      mtext(paste("Skewness=",round(skewness(quant.vec, type=2),3)), side=3, line=1, outer=outer, at=max1-3*(max1-min1)/9, cex=mcex, adj=0)      
      mtext(paste("LL=",round(CI[1],3)),side=3,line=2,outer=outer, at=max1-3*(max1-min1)/9, cex=mcex, adj=0)
      mtext(paste("UL=",round(CI[2],3)),side=3,line=3,outer=outer,at=max1-3*(max1-min1)/9, cex=mcex, adj=0)
      mtext("Monte Carlo CI",side=3,line=4,outer=outer,at=max1-3*(max1-min1)/9, cex=mcex, font=2, adj=0)
    }
    
    if(type=="all"){
      #New- 1/24/14-DT
      mtext(paste("Kurtosis=",round(kurtosis(quant.vec, type=2),3)), side=3,line=1,outer=outer, at=max1-3*(max1-min1)/9, cex=mcex, adj=0)
      mtext(paste("Skewness=",round(skewness(quant.vec, type=2) ,3)) ,side=3,line=2,outer=outer, at=max1-3*(max1-min1)/9, cex=mcex, adj=0)
      mtext(paste("LL=", round(CI[1],3)), side=3,line=3, outer=outer, at=max1-3*(max1-min1)/9, cex=mcex, adj=0)
      mtext(paste("UL=", round(CI[2],3)), side=3,line=4, outer=outer, at=max1-3*(max1-min1)/9, cex=mcex, adj=0)
      mtext("Monte Carlo", side=3,line=5, outer=outer, at=max1-3*(max1-min1)/9, cex=mcex, adj=0)
      if (res.asymp$SE>40*.Machine$double.eps) legend(x=max1,y= par("usr")[4],c("Monte Carlo","Asymptotic Normal"),col=c("black","blue"),lty=c(1,2),lwd=c(2:2),bty="n",title="", cex=mcex, y.intersp=.5, xpd=FALSE, xjust=.5) 
    }
    
    if(plotCI){
      yci<-par("usr")[3]-1.2*diff(par("usr")[3:4])/25
      arrows(CI[1],yci,CI[2],yci,length=0,angle=90,code=3,cex=1.5,lwd=2)
      points(quantMean,yci,pch=19,cex=1.5)
    }
    
  } #end of if plot
  
  res <- list(CI,Estimate=quantMean,SE=quantSE,MCError=quantError,p= pMinH1) #Results
  attr(res,"quant")  <- q1
  #names(res) <- c( paste((1-alpha)*100,"% ", "CI",sep=""), "Estimate", "SE","MC Error", "p")
  if(H0) return(list(res,H0res))
  else return(res)
}
