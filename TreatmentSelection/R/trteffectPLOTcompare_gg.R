trteffectPLOTcompare_gg <-
function(x1, x2, ci, ci.bounds, get.F, fixed.values, conf.bands,  rho, xlab, ylab, xlim, ylim, main, markerTWO=FALSE, lty = 1,  p=NULL){ 
  
  mylim <- NULL
  trt.effect <- x1$derived.data$trt.effect

  event <- x1$derived.data$event
  trt <- x1$derived.data$trt
  
  trt.effect2 <- x2$derived.data$trt.effect

  event2 <- x2$derived.data$event
  trt2 <- x2$derived.data$trt
  
  n = length(trt.effect)
  
  F.D <- get.F(trt.effect, event, trt, rho = rho)*100
  F.D2 <- get.F(trt.effect2, event2, trt2, rho = rho)*100

  mydata = data.frame(trt.effect, F.D, lty = 1, size =1)
  mydata = mydata[with(mydata, order(F.D)),]
  
  mydata2 = data.frame(trt.effect = trt.effect2, F.D = F.D2,lty =  2, size = 1 )
  mydata2 = mydata2[with(mydata2, order(F.D)),]
  mydata <- rbind(mydata, mydata2)
 # mydata <- rbind(mydata, mydata2, c(-100, -100, 3,.5), c(-100, 100, 4, .5))
  
  
  ## need to adjust these for scc and cc sampling designs
  allMeasures <- x1$functions$get.summary.measures( x1$derived.data, rho, x1$model.fit$thresh)
  
  avglines <- cbind(0, sort(F.D), 4, .5)
  avglines <- rbind(avglines, 
                    cbind(allMeasures$ER.trt0.emp-allMeasures$ER.trt1.emp, 
                          sort(F.D), 3, .5))
  avglines = data.frame(avglines); names(avglines) = names(mydata)
  mydata <- rbind(mydata, avglines)

   if(is.null(xlab)) xlab <- "% population below treatment effect"
   if(is.null(ylab)) ylab <- "treatment effect"
   if(is.null(xlim)) xlim <- c(0,100)
   if(is.null(ylim)) ylim <- mylim
   if(is.null(main)) main <- "Treatment effect distribution"
   p <- ggplot(mydata) 
  
       
     #add x/y labels and main
     p <- p + xlab(xlab) + ylab(ylab) + ylim(ylim[1], ylim[2])+ ggtitle(main) 
     
    # p <- p + stat_hline(yintercept  = mean(trt.effect), aes(linetype = factor(3), size = factor(3)), show.guide = FALSE)+
    #   stat_hline(yintercept = 0, aes( linetype = factor(4), size = factor(4)), show.guide = FALSE)
     
     p <- p + theme( text = element_text(size=14)) #, 
     
     
     p <- p + scale_x_continuous(limits = xlim)
     
     
     
     
   
  
  


  if(!is.null(ci.bounds)){
 

  fixed.values1 <- fixed.values[1,]
  fixed.values2 <- fixed.values[2,]
  
  
  if(substr(ci, 1,1)=="v"){
    
    index.fix1  <- (fixed.values1<= max(F.D) & fixed.values1 >= min(F.D)) 
    index.fix2 <- (fixed.values2<= max(F.D2) & fixed.values2 >= min(F.D2))
    
    width = 5
    
  }else{
    index.fix1  <- (fixed.values1<= max(trt.effect) & fixed.values1 >= min(trt.effect)) 
    index.fix2  <- (fixed.values2<= max(trt.effect2) & fixed.values2 >= min(trt.effect2)) 
    width = .05
  }

  p <- shade_gg(p, ci.bounds[1:2,index.fix1], fixed.values1[index.fix1], type = substr(ci, 1, 1), bands = conf.bands, lty=1, width = width)
  p <- shade_gg(p, ci.bounds[3:4,index.fix2], fixed.values2[index.fix2], type = substr(ci, 1, 1), bands = conf.bands, lty=2, width = width)

  }


    p <- p+geom_step(data =  mydata[(1:(2*n)),], aes(x = F.D, y = trt.effect, linetype = factor(lty), size = factor(lty)), direction = "vh")
    p <- p+geom_line(data =  mydata[-c(1:(2*n)),], aes(x = F.D, y = trt.effect, linetype = factor(lty), size = factor(lty)))
       


  return(list(p=p))
}
