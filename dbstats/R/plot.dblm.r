

 ###################
 #### plot.dblm ####
 ###################

 ## Description: 
 ##    generic function plot for a dblm object. Six plots (selected by which) 
 ##    are available:
 ##       - a plot of residual vs fitted values.
 ##       - the Q-Qplot of normality.
 ##       - a Scale-Location plot of residuals against fitted values.  
 ##       - cook.distance
 ##       - Leverage vs residuals
 ##        - minimum effective rank of OCV, GCV, AIC or BIC method.
 ##
        

plot.dblm<-function(x,which=c(1L:3L, 5L),id.n=3,main="",cook.levels = c(0.5, 1),
              cex.id = 0.75,type.pred=c("link","response"),...){
                          
     # stop if the object is not a dblm object.
     if (!inherits(x, "dblm")&&!inherits(x, "dbglm")) 
        stop("use only with \"dblm\" or \"dbglm\" objects")
     
     # auxiliar boolean vector variable with 4 items. True if the plot is 
     # selected by which or false if not. 
     if (!is.numeric(which) || any(which < 1) || any(which > 6)) 
        stop("'which' must be in 1:6")
             
     show <- rep(FALSE, 6)
     show[which] <- TRUE
    
     # change the panel like lm plot (with mouse clik or with enter).
     # if only one plot is selected ask=FALSE --> not need to change the panel.  
     ask = prod(par("mfcol")) < length(which) 
      one.fig <- prod(par("mfcol")) == 1
     if (ask) {
        oask <- devAskNewPage(TRUE)  # ask for new page
        on.exit(devAskNewPage(oask)) # exit the format plot (if true) 
     }


     # parameters shared by all graphics: 
     iid <- 1L:id.n       # vector 1:(number of observations to hightlihgt)  
    
      
     # residual values
     if (inherits(x, "dblm")) 
      r<- x$residuals      
     if (inherits(x, "dbglm"))
      r<- summary(x)$deviance.resid
    
     # fitted values
     if (inherits(x, "dblm"))
      yh<- x$fitted.values
     if (inherits(x, "dbglm")){
      type.pred <- match.arg(type.pred)
      if (type.pred=="response")
        yh <- x$fitted.values
      else 
        yh <- attr(x,"eta")
     }
     
     # weights
     if (inherits(x, "dblm"))
       w<-x$weights
     if (inherits(x, "dbglm"))
      w<-attr(x,"ori_weights")
    
     # weighting the residuals  
     r.w <- if (is.null(w))
            r
       else sqrt(w) * r
     
     # hat values         
     hii<-diag(x$H)+x$weights/sum(x$weights)   
     
     # residual standard desviation
     if (inherits(x, "dblm")) 
      s<-summary(x)$sigma  
     if (inherits(x, "dbglm"))
      s<-sqrt(summary(x)$dispersion)
     
     # standaritzed residuals 
     rs <- r.w/(s * sqrt(1-hii)) 
     
     # takes the elements to be highlighted in the standaritzed residuals plots
     show.rs <- sort.list(abs(rs), decreasing = TRUE)[iid] 
     
     # "Residual vs fitted plot"
     if(show[1L]){
       # y dimensions 
       ylim <- range(r, na.rm = TRUE)
       if (id.n > 0) 
          # extended the range of y if id.n > 0
          ylim <- extendrange(r = ylim, f = 0.08)
      
      if ((main==""&&length(which)==1)||length(which)>1)
        main="Residuals vs fitted"
       
       # plot with yh (fitted values) and r (residual values)
       plot(yh,r,ylim=ylim,ylab="Residuals",xlab="Fitted values",main=main) 
       abline(h = 0, lty = 3, col = "gray")
       lines(lowess(yh,r),col="red")   
       if (id.n > 0) {
            ro<-sort(r) 
            show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
            y.id <- r[show.r]
            y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
            # highlighting the id.n higher residuals.
            text(yh[show.r],pos=3,cex=0.7, y.id, show.r)      
       }
      }
     
     
     # "Q-Qnorm plot"
     if(show[2L]){
      
      ylim <- range(rs, na.rm = TRUE)
      ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
      
      if ((main==""&&length(which)==1)||length(which)>1)
       main<-"Normal Q-Q"
     
      ylab23 <- "Standardized residuals"
      
      # call qqnorm function. See the normality of the standaritzed residuals.        
      qq <- qqnorm(rs, main = main, ylab = ylab23, ylim = ylim, 
            ...) 
      qqline(rs, lty = 3, col = "gray50")
      if (id.n > 0){
           ro<-sort(r) 
           show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
           y.id <- r[show.r]
           y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
           # highlighting the id.n higher standard residuals.
           text(x=qq$x[show.rs],y=qq$y[show.rs],pos=3,cex=0.7,adj=y.id,
                  labels=show.rs) 
     }
    }

     # "Scale-location plot"   
      if(show[3L]){
        ylab23 <- "Standardized residuals"   
        sqrtabsr <- sqrt(abs(rs))
        # y dimensions (>0)
        ylim <- c(0, max(sqrtabsr, na.rm = TRUE)) 
        if (id.n > 0) 
          # extended the range of y if id.n > 0
          ylim <- extendrange(r = ylim, f = 0.08)
        #y label. sqrt(abs(rs))
        yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name(ylab23)))) 
        
        if ((main==""&&length(which)==1)||length(which)>1)
         main<- "Scale-location"
         
       # fitted values (weight>0) 
        yhn0 <- if (is.null(w))  
            yh
        else yh[w != 0]
        # plot between fitted values and the square root of the absolute(stand.resiudals)
        plot(yhn0, sqrtabsr, xlab ="Fitted values",ylab=yl, main = main,  
            ylim = ylim)
        lines(lowess(yhn0,sqrtabsr),col="red") 
        if (id.n > 0) 
            # highlighting the id.n higher standard residuals.
            text(yhn0[show.rs],pos=3,cex=0.7,sqrtabsr[show.rs],show.rs)
    }    

     # "obs. number vs cooks distance" 
     if(show[4L]){
     
      num2<-hii/(1-hii)^2
      p=x$eff.rank+1

      if (inherits(x, "dbglm")){
         r.w <- summary(x)$pears.resid  
        cook<-(r.w^2/(p*s^2)*num2) # cook's distance (s'ha de multiplicar per p)
      }else
        cook<-((x$residuals)^2/(p*s^2)*num2) 
      
      if ((main==""&&length(which)==1)||length(which)>1)
       main="Cook's distance"
      
      if (id.n > 0) {
            show.r <- order(-cook)[iid]   # the most extreme points.
            ymx <- cook[show.r[1L]] * 1.075  # extended ylim
      } else ymx <- max(cook, na.rm = TRUE)
      
      plot(cook, type = "h", ylim = c(0, ymx), main = main, 
            xlab = "Obs. number", ylab = "Cook's distance", ...)
      if (id.n > 0)
        #laballed the extreme points  
        text(show.r,pos=3,cex=0.7,cook[show.r],show.r) 
      }
      
     
      # "leverage vs standaritzed residuals (cooks distance)" 
     if(show[5L]){
      
       if (inherits(x, "dbglm")){
         r.w <- summary(x)$pears.resid # standaritzed residuals
         rsp<- r.w/(s * sqrt(1 - hii))            
         num2<-hii/(1-hii)^2
         cook<-(r.w^2/s^2*num2)  # cook's distance
         ylab = "Std. Pearson resid"
       }
       else{
         rsp <- rs         # standaritzed residuals
         num2<-hii/(1-hii)
         cook<-((x$residuals)^2/s^2*num2) 
         ylab = "Standaritzed residuals"
       } 
       
       ylim <- range(rsp, na.rm = TRUE) # y lim 
       if (id.n > 0) {
           ylim <- extendrange(r = ylim, f = 0.08) # extended the y lim
             show.rsp <- order(-cook)[iid]         # extreme points
        }

       if ((main==""&&length(which)==1)||length(which)>1)
        main="Residuals vs Leverage"
       
       # plot of leverage vs residual standaritzed 
       plot(hii,rsp, main = main,ylim=ylim, xlim = c(0, max(hii, na.rm = TRUE)),      
            ylab = ylab, xlab = "Leverage", ...) 
       abline(h = 0, v = 0, lty = 3, col = "gray") # lines in the axes
       lines(lowess(hii,rsp),col="red") 
       usr <- par("usr")
       r.hat <- range(hii, na.rm = TRUE)
       hh <- seq.int(min(r.hat[1L], r.hat[2L]/100),
                  usr[2L], length.out = 101)
       for (crit in cook.levels) {
                  cl.h <- sqrt(crit * (x$eff.rank+1) * (1 - hh)/hh)
                  lines(hh, cl.h, lty = 2, col = 2)
                  lines(hh, -cl.h, lty = 2, col = 2)
                }  
       text(hii[show.rsp],pos=3,cex=0.7,rsp[show.rsp],show.rsp)
       
       legend("bottomleft", legend = "Cook's distance",
                  lty = 2, col = 2, bty = "n")
                xmax <- min(0.99, usr[2L])
                ymult <- sqrt((x$eff.rank+1) * (1 - xmax)/xmax)
                aty <- c(-sqrt(rev(cook.levels)) * ymult, sqrt(cook.levels) *
                  ymult)
                axis(4, at = aty, labels = paste(c(rev(cook.levels),
                  cook.levels)), mgp = c(0.25, 0.25, 0), las = 2,
                  tck = 0, cex.axis = cex.id, col.axis = 2)
     
      }
       
    # "rank effective against ocv, gcv, aic or bic estimator"  
     if(show[6L]){
    
        # range.eff.rank
        if (inherits(x, "dbglm"))
         range.eff.rank <-  c(attr(x,"range.eff.rank")[1]:attr(x,"range.eff.rank")[2])
        else{
         if (inherits(x, "dblm")) 
          range.eff.rank <-  1: attr(x,"threshold")    
         }      

      #if (inherits(x,"dbglm")) stop("this plot is not provided for a dbglm object")
      if (!attr(x,"full.search"))  stop("this plot is not provided for a dblm without 'full.search=TRUE'")
      else{
       method<-attr(x,"method")
       # if the user method is eff.rank or epsilonis not posible to makes this plot.
       if (method=="eff.rank"||method=="epsilon"){  
            warning("the effective rank election plot can only be done, if the object dblm has been called with the parameter method different that 'eff.rank' or 'epsilon'")
       }else{
        
        if (method=="OCV"){
         ocvs<-attr(x,"ocvs")
         # red for the selected eff.rank
         color<-c(rep("black",(x$eff.rank-1)),"red",rep("black",(length(ocvs)-x$eff.rank)))[range.eff.rank]
         ocvs <-  ocvs[range.eff.rank]
         
        if ((main==""&&length(which)==1)||length(which)>1)
          main="OCV vs Effective rank"
        
        ylim<-c(0,max(ocvs))
         
         if (id.n > 0) 
           # extended the range of y if id.n > 0
           ylim <- extendrange(r = ylim, f = 0.08) 
         # plot with the optimal effective rank of OCV method   
         plot(range.eff.rank, ocvs, type="h",ylim = ylim, main = main,col=color,
            xlab = "effective rank", ylab = "Ordinary Cross-validation")
         # identify the selected eff.rank   
         text(x$eff.rank,pos=3,col="red",cex=0.7,x$ocv,x$eff.rank)           
        }
       
        if (method=="GCV"){
         gcvs<-attr(x,"gcvs")
         color <-c(rep("black",(x$eff.rank-1)),"red",rep("black",(length(gcvs)-x$eff.rank)))[range.eff.rank]
         gcvs <-  gcvs[range.eff.rank]
          
         
        if ((main==""&&length(which)==1)||length(which)>1)
          main="GCV vs Effective rank"
         ylim<-c(0,max(gcvs))
         
         if (id.n > 0)
          # extended the range of y if id.n > 0 
          ylim <- extendrange(r = ylim, f = 0.08)
         
         # plot with the optimal effective rank of GCV method 
         plot( range.eff.rank, gcvs, type="h",ylim = c(0,max(gcvs)), main = main,col=color,  
            xlab = "effective rank", ylab = "Generalized Cross-validation")
         text(x$eff.rank,pos=3,col="red",cex=0.7,x$gcv,x$eff.rank)
        }
       
        if (method=="AIC"){
         aics<-attr(x,"aics")
         color<-c(rep("black",(x$eff.rank-1)),"red",rep("black",(length(aics)-x$eff.rank)))[range.eff.rank]
         aics<-aics[range.eff.rank]
         
         
         if ((main==""&&length(which)==1)||length(which)>1)
          main="AIC vs Effective rank"
         ylim <- range(aics, na.rm = TRUE)
         
         if (id.n > 0)
          # extended the range of y if id.n > 0 
          ylim <- extendrange(r = ylim, f = 0.08)
         # plot with the optimal effective rank of AIC method 
         plot(range.eff.rank, aics, type="h", main = main,col=color,         
            xlab = "effective rank",ylim=ylim, ylab = "AIC criteria")
         text(x$eff.rank,pos=1,col="red",cex=0.7,aics[x$eff.rank],x$eff.rank)  
        }
        
        if (method=="BIC"){
         bics<-attr(x,"bics")
         color<-c(rep("black",(x$eff.rank-1)),"red",rep("black",(length(bics)-x$eff.rank)))[range.eff.rank]
         bics<-bics[range.eff.rank]
        
         if ((main==""&&length(which)==1)||length(which)>1)
          main="BIC vs Effective rank"
         ylim <- range(bics, na.rm = TRUE)
        
         if (id.n > 0) 
          # extended the range of y if id.n > 0
          ylim <- extendrange(r = ylim, f = 0.08) 
          
         # plot with the optimal effective rank of BIC method  
         plot(range.eff.rank, bics, type="h", main = main,col=color,     
            xlab = "effective rank",ylim=ylim, ylab = "BIC criteria")
         text(x$eff.rank,pos=1,col="red",cex=0.7,bics[x$eff.rank],x$eff.rank)  
        }        
      }    
     }
    }
  return (invisible())
}