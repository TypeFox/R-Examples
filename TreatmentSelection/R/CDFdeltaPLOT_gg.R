CDFdeltaPLOT_gg <-
function(x, ci, ci.bounds, get.F, fixed.values,conf.bands,  rho, xlab, ylab, xlim, ylim, main){ 

  trt.effect <- x$derived.data$trt.effect

  event <- x$derived.data$event
  trt <- x$derived.data$trt
  n = length(trt.effect)

  F.D <- get.F(trt.effect, event, trt, rho = rho)*100

  mydata = data.frame(trt.effect, F.D )
  mydata = mydata[with(mydata, order(F.D)),]

  

    cen <- mean(c(min(trt.effect, na.rm=TRUE), max(trt.effect, na.rm=TRUE)))
    ran <- max(trt.effect, na.rm=TRUE) - min(trt.effect, na.rm=TRUE)
    
    if(substr(ci, 1,1) =="h"){ cen <- mean(c(min(c(trt.effect, ci.bounds), na.rm=TRUE), max(c(trt.effect, ci.bounds), na.rm=TRUE))); ran <- max(c(trt.effect, ci.bounds), na.rm=TRUE) - min(c(trt.effect, ci.bounds), na.rm=TRUE)}
    
    
    ran <- ran*1.1
    mylim <- c(cen-ran/2, cen+ran/2)
    
    if(is.null(ylab)) ylab <- "% population below treatment effect"
    if(is.null(xlab)) xlab <- "treatment effect"
    if(is.null(ylim)) ylim <- c(0,100)
    if(is.null(xlim)) xlim <- mylim
    if(is.null(main)) main <- "Treatment effect distribution"
    p <- ggplot(mydata)     

 # breaks = seq(xlim[1], xlim[2], length.out = 5)
  
  
  if(!is.null(ci.bounds)){
    
    ci.bounds <- matrix(ci.bounds, ncol=length(fixed.values), nrow = 2)
    
    if(substr(ci, 1,1)=="h"){
      width = 5
      index.fix  <- (fixed.values<= max(F.D) & fixed.values >= min(F.D)) 
    }else{
      width = .05
      index.fix  <- (fixed.values<= max(trt.effect) & fixed.values >= min(trt.effect)) 
    }
    
    p <- shade_gg(p, ci.bounds[,index.fix], fixed.values[index.fix], type = substr(ci, 1, 1), bands = conf.bands, width = width)
  }
  
  
      
    p <- p+geom_step(data = mydata, aes(y = F.D, x = trt.effect), direction = "vh", size = 1)
    
    #add x/y labels and main
    p <- p + xlab(xlab) + ylab(ylab) + ylim(ylim[1], ylim[2])+ ggtitle(main) 
    
    #change the names for the legend
      
    p <- p + theme( text = element_text(size=14)) #, 
    
    
    p <- p + scale_x_continuous( limits = xlim)
    
  
  
  
  print(p)
 list(p=p)
}
