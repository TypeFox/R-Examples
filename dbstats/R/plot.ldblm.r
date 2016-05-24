

 ####################
 #### plot.ldblm ####
 ####################

 ## Description: 
 ##    generic function plot for a ldblm or ldbglm object.  
 ##    
 ##       - a plot of fitted values vs response
 ##       - a plot of residual vs fitted values.
 ##       - a plot of bandwidth vs method statisticals.
 ##
        

plot.ldblm <-function(x,which=c(1L,2L),id.n=3,main="",...){
     
     # stop if the object is not a dblm object.
     if (!inherits(x, "ldblm")&&!inherits(x, "ldbglm")) 
      stop("use only with \"ldblm\" or \"ldbglm\" objects")
     
     # auxiliar boolean vector variable with 4 items. 
     # True if the plot is selected by which or false if not. 
     show <- rep(FALSE, 3)
     show[which] <- TRUE
    
     # change the panel like lm plot (with mouse clik or with enter) 
     # if only one plot is selected ask=FALSE --> not need to change the panel.
     ask = prod(par("mfcol")) < length(which)  
      one.fig <- prod(par("mfcol")) == 1
     if (ask) {
        oask <- devAskNewPage(TRUE)  # ask for new page
        on.exit(devAskNewPage(oask)) # exit the format plot (if true) 
     }
    
    yhat<-x$fitted.values 
    iid <- 1L:id.n  
    y<-x$y
    
    if(show[1L]){ 
     # Fitted values vs response
      if ((main==""&&length(which)==1)||length(which)>1)
       main<-"Fitted values vs response" 
     
     ylim <- range(yhat, na.rm = TRUE)
     if (id.n > 0) 
       ylim <- extendrange(r = ylim, f = 0.08) #extended the range of y if id.n > 0  
     
     plot(y,yhat,main=main,ylim=ylim,ylab="fitted.values",xlab="y")
     abline(0,1,lty = 3, col = "red")
     
      if (id.n > 0) {
            show.y <- sort.list(abs(x$residuals), decreasing = TRUE)[iid]
            y.id <- yhat[show.y]
            # highlighting the id.n higher residuals.
            text(y[show.y],pos=3,cex=0.7, y.id, show.y) 
       }
   }
   
   if(show[2L]){ 
     # Residuals vs fitted
     if ((main==""&&length(which)==1)||length(which)>1)
      main<-"Residuals vs fitted"
     
     ylim <- range(x$residuals, na.rm = TRUE)
     if (id.n > 0) 
       # extended the range of y if id.n > 0
       ylim <- extendrange(r = ylim, f = 0.08) 
     
     plot(yhat,x$residuals,main=main,ylim=ylim,ylab="Residuals",
        xlab="fitted values")
     abline(h=0,lty = 3, col = "red")
     
      if (id.n > 0) {
            show.y <- sort.list(abs(x$residuals), decreasing = TRUE)[iid]
            y.id <- x$residuals[show.y]
            # highlighting the id.n higher residuals. 
            text(yhat[show.y],pos=3,cex=0.7, y.id, show.y)      
       }
    }
    
      
    if(show[3L]){
      # bandwidth vs method statisticals
      method.h<-attr(x,"method.h") # method of ldblm call
      noh<-attr(x,"noh")       # number of possible bandwidths
      h_vec<-attr(x,"h_vec")   # bandwidth used in ldblm
       
      # Ordinary Cross validation method
      if(method.h=="OCV"){
          
         ocvs<-attr(x,"OCV")              # ocv's for each bandwidth
         ocv_opt<-attr(x,"OCV_opt")       # optimal ocv
         
         # color vector, only the optimal ocv becomes red
         color<-c(rep("black",(noh)))   
         color[which(ocv_opt==ocvs)]<-"red"
         
         if ((main==""&&length(which)==1)||length(which)>1)
          main="Bandwidth h of OCV method.h"
         ylim<-c(0,max(ocvs))
         
         # plot with the optimal bandwidth of OCV method.h
         plot(h_vec,ocvs,type="h", main = main,col=color,  
            xlab = "bandwidth h", ylab = "Ordinary Cross-validation")
         # identify the selected h
         text(x$h.opt,pos=3,col="red",cex=0.7,ocv_opt,round(x$h.opt,4)) 
      }
     
      # Generalized Cross validation method.h
      if(method.h=="GCV"){
        
         gcvs<-attr(x,"GCV")
         gcv_opt<-attr(x,"GCV_opt")
         
         color<-c(rep("black",(noh))) 
         color[which(gcv_opt==gcvs)]<-"red"
         
        if ((main==""&&length(which)==1)||length(which)>1)
          main="Bandwidth h of GCV method.h"
     
         
         # plot with the optimal bandwidth of GCV method.h                        
         plot(h_vec,gcvs,type="h", main = main,col=color,     
            xlab = "bandwidth h", ylab = "Generalized Cross-validation")
         # identify the selected h
         text(x$h.opt,pos=3,col="red",cex=0.7,gcv_opt,round(x$h.opt,4))   
      }                                                                                       
      
        
      if(method.h=="AIC"){
         aics<-attr(x,"AIC")
         aic_opt<-attr(x,"AIC_opt")
         
         color<-c(rep("black",(noh))) 
         color[which(aic_opt==aics)]<-"red"
        
         if ((main==""&&length(which)==1)||length(which)>1)
          main="Bandwidth h of AIC method.h"
         
         
         # plot with the optimal bandwidth of AIC method.h
         plot(h_vec,aics,type="h", main = main,col=color,      
            xlab = "bandwidth h", ylab = "Aikaike Information Criterium")
         # identify the selected h
         text(x$h.opt,pos=3,col="red",cex=0.7,aic_opt,round(x$h.opt,4))     
        }
        
       if(method.h=="BIC"){
         bics<-attr(x,"BIC")
         bic_opt<-attr(x,"BIC_opt")
         
         color<-c(rep("black",(noh))) 
         color[which(bic_opt==bics)]<-"red"
         
         if ((main==""&&length(which)==1)||length(which)>1)
          main="Bandwidth h of BIC method.h"
   
         
         # plot with the optimal bandwidth of BIC method.h
         plot(h_vec,bics,type="h", main = main,col=color,       
            xlab = "bandwidth h", ylab = "Bayesian Information Criterium")
         # identify the selected h 
         text(x$h.opt,pos=3,col="red",cex=0.7,bic_opt,round(x$h.opt,4))     
        }
        if(method.h=="user.h")
         warning("there is not plot of bandwidth selection if method.h=user.h" )     
      }
      return (invisible())
 }