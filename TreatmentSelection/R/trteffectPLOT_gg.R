trteffectPLOT_gg <-
function(x, ci, ci.bounds, get.F, fixed.values, conf.bands,  rho, xlab, ylab, xlim, ylim, main, markerTWO=FALSE, lty = 1,  p=NULL){ 
  

  trt.effect <- x$derived.data$trt.effect

  event <- x$derived.data$event
  trt <- x$derived.data$trt
  n = length(trt.effect)
  
  F.D <- get.F(trt.effect, event, trt, rho = rho)*100
  
  mydata = data.frame(trt.effect, F.D, lty )
  mydata = mydata[with(mydata, order(F.D)),]

  
  
  ## need to adjust these for scc and cc sample designs. 
  allMeasures <- x$functions$get.summary.measures( x$derived.data, rho, x$model.fit$thresh)
  
  
  avglines <- cbind(0, sort(F.D), 4)
  avglines <- rbind(avglines,
                    cbind(allMeasures$ER.trt0.mod - allMeasures$ER.trt1.mod,
                          sort(F.D), 3))
  avglines = data.frame(avglines); names(avglines) = names(mydata)
  mydata <- rbind(mydata, avglines)
  

 if(!markerTWO){
   cen <- mean(c(min(trt.effect, na.rm=TRUE), max(trt.effect, na.rm=TRUE)))
   ran <- max(trt.effect, na.rm=TRUE) - min(trt.effect, na.rm=TRUE)
   
   if(substr(ci, 1,1) =="v"){ cen <- mean(c(min(c(trt.effect, ci.bounds), na.rm=TRUE), max(c(trt.effect, ci.bounds), na.rm=TRUE))); ran <- max(c(trt.effect, ci.bounds), na.rm=TRUE) - min(c(trt.effect, ci.bounds), na.rm=TRUE)}
   
   
   ran <- ran*1.1
   mylim <- c(cen-ran/2, cen+ran/2)
   
   if(is.null(xlab)) xlab <- "% population below treatment effect"
   if(is.null(ylab)) ylab <- "treatment effect"
   if(is.null(xlim)) xlim <- c(0,100)
   if(is.null(ylim)) ylim <- mylim
   if(is.null(main)) main <- "Treatment effect distribution"
   p <- ggplot()     

     #add x/y labels and main
     p <- p + xlab(xlab) + ylab(ylab) + ylim(ylim[1], ylim[2])+ ggtitle(main) 
     
     p <- p + theme( text = element_text(size=18)) #, 
     
     p <- p + scale_x_continuous(limits = xlim)
   
  }

  if(!is.null(ci.bounds)){
 
  ci.bounds <- matrix(ci.bounds, ncol=length(fixed.values), nrow = 2)

  if(substr(ci, 1,1)=="v"){

    index.fix  <- (fixed.values<= max(F.D) & fixed.values >= min(F.D))
    width = 5
  }else{
    width = .05
    index.fix  <- (fixed.values<= max(trt.effect) & fixed.values >= min(trt.effect)) 
  }
  
  p <- shade_gg(p, ci.bounds[,index.fix], fixed.values[index.fix], type = substr(ci, 1, 1), bands = conf.bands, lty, width = width)
  }

  


 
    p <- p+geom_step(data =  mydata[(1:(n)),], aes(x = F.D, y = trt.effect), size = 1, direction = "vh")
    p <- p+geom_line(data =  mydata[-c(1:(n)),], aes(x = F.D, y = trt.effect, linetype = factor(lty)), size = 0.5)
    

  p <- p+ scale_linetype_manual(name = "Treatment Effect", breaks = c( "3", "4"), values = c( 3, 4), labels = c("Mean", "Zero"))
  

  print(p) 
    
  return(list(p=p))
}
