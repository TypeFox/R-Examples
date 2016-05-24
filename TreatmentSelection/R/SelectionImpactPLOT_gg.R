SelectionImpactPLOT_gg <-
function(x, ci, ci.bounds, get.F, fixed.values,conf.bands,  rho, xlab, ylab, xlim, ylim, main){ 
  
  risk.t0 <- x$derived.data$fittedrisk.t0
  risk.t1 <- x$derived.data$fittedrisk.t1
  trt.effect <- x$derived.data$trt.effect

  event <- x$derived.data$event
  trt <- x$derived.data$trt
  n = length(trt)

  F.D <- get.F(trt.effect, event, trt, rho = rho)*100
  theta.curve <- EventRateVec(risk.t0, risk.t1, F.D, rho, event, trt)
  
  lty = 1 #define the lty for the main curve
  mydata = data.frame(theta.curve, F.D, lty )
  mydata = mydata[with(mydata, order(F.D)),]

  ## need the mean risk given t=0 and t=1. We need to account for subsampling in this. 
  allMeasures <- x$functions$get.summary.measures( x$derived.data, rho, x$model.fit$thresh)
  
  
  avglines <- cbind(allMeasures$ER.trt0.mod, sort(F.D), 4)
  avglines <- rbind(avglines, 
              cbind(allMeasures$ER.trt1.mod, sort(F.D), 3))
  
  avglines = data.frame(avglines); names(avglines) = names(mydata)
  mydata <- rbind(mydata, avglines)
  
  
    cen <- mean(c(min(theta.curve, na.rm=TRUE), max(theta.curve, na.rm=TRUE)))
    ran <- max(theta.curve, na.rm=TRUE) - min(theta.curve, na.rm=TRUE)
    
    if(substr(ci, 1,1) =="v"){ cen <- mean(c(min(c(theta.curve, ci.bounds), na.rm=TRUE), max(c(theta.curve, ci.bounds), na.rm=TRUE))); ran <- max(c(theta.curve, ci.bounds), na.rm=TRUE) - min(c(theta.curve, ci.bounds), na.rm=TRUE)}
    
    
    ran <- ran*1.1
    mylim <- c(cen-ran/2, cen+ran/2)
    
    if(is.null(xlab)) xlab <- "d = % population below treatment effect"
    if(is.null(ylab)) ylab <- "Event rate given trt rule: T = 1 if F(v) > d"
    if(is.null(xlim)) xlim <- c(0,100)
    if(is.null(ylim)) ylim <- mylim
    if(is.null(main)) main <- "Selection impact curve"
    p <- ggplot()     
    
    #add x/y labels and main
    p <- p + xlab(xlab) + ylab(ylab) + ylim(ylim[1], ylim[2])+ ggtitle(main) 
    
    p <- p + theme( text = element_text(size=18)) #, 
    
    p <- p + scale_x_continuous(limits = xlim)
    
  
  
  if(!is.null(ci.bounds)){
    
    ci.bounds <- matrix(ci.bounds, ncol=length(fixed.values), nrow = 2)
    
    if(substr(ci, 1,1)=="v"){
      
      index.fix  <- (fixed.values<= max(F.D) & fixed.values >= min(F.D))
      width = 5
    }else{
      width = .05
      index.fix  <- (fixed.values<= max(theta.curve) & fixed.values >= min(theta.curve)) 
    }
    
    p <- shade_gg(p, ci.bounds[,index.fix], fixed.values[index.fix], type = substr(ci, 1, 1), bands = conf.bands, lty, width = width)
  }
  
  
  
  
  
  p <- p+geom_step(data =  mydata[(1:(n)),], aes(x = F.D, y = theta.curve), size = 1, direction = "vh")
  p <- p+geom_line(data =  mydata[-c(1:(n)),], aes(x = F.D, y = theta.curve, linetype = factor(lty)), size = 0.5)
  
  
  p <- p+ scale_linetype_manual(name = "Event Rate", breaks = c( "3", "4"), values = c( 3, 4), labels = c("treat all", "treat none"))
  
  
  print(p) 
  
  return(list(p=p))

}
