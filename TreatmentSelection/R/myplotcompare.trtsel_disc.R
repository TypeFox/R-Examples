myplotcompare.trtsel_disc <-
function(x, bootstraps = 500, alpha = .05,
                           ci   = "horizontal", marker.names = c("Marker 1", "Marker 2"),  
                           xlab = NULL, 
                           ylab = NULL, 
                           xlim = NULL, 
                           ylim = NULL, 
                           main = NULL, offset = offset, conf.bands,  annotate.plot)
{
  quantile <- NULL #appease check
  ts1 <- x$trtsel1
  ts2 <- x$trtsel2
  fittedrisk.t0.y1 <- ts1$derived.data$fittedrisk.t0
  fittedrisk.t1.y1 <- ts1$derived.data$fittedrisk.t1
  marker1 <- ts1$derived.data$marker
  delta.y1 <- ts1$derived.data$trt.effect
  link <- ts1$model.fit$link
  fittedrisk.t0.y2 <- ts2$derived.data$fittedrisk.t0
  fittedrisk.t1.y2 <- ts2$derived.data$fittedrisk.t1
  marker2 <- ts2$derived.data$marker
  delta.y2 <-  ts2$derived.data$trt.effect
  
  n<-length(marker1)
  
  rho  <- ts1$model.fit$cohort.attributes
  study.design <- ts1$model.fit$study.design
  trt <- ts1$derived.data$trt
  event <- ts1$derived.data$event
  
  boot.sample <- ts1$functions$boot.sample

  
  
  
  
  one.boot.plot_disc <-
    function(event, trt, marker1, marker2, rho = rho,  obp.boot.sample){
      
      myboot.sample <- obp.boot.sample( event, trt, rho)
      
      rho.b <- myboot.sample[1:7]
      ind   <- myboot.sample[-c(1:7)]
      
      event.b  <- event[ind]
      trt.b  <- trt[ind]
      marker1.b  <- marker1[ind] 
      marker2.b  <- marker2[ind] 
      mval1 <- sort(unique(marker1))
      mval2 <- sort(unique(marker2))
      
      
      c(   trteff.mkr10 =mean(event.b[trt.b==0 & marker1.b ==mval1[1]]) - mean(event.b[trt.b==1 & marker1.b ==mval1[1]]), 
           trteff.mkr11 =mean(event.b[trt.b==0 & marker1.b ==mval1[2]]) - mean(event.b[trt.b==1 & marker1.b ==mval1[2]]),
           trteff.mkr20 =mean(event.b[trt.b==0 & marker2.b ==mval2[1]]) - mean(event.b[trt.b==1 & marker2.b ==mval2[1]]), 
           trteff.mkr21 =mean(event.b[trt.b==0 & marker2.b ==mval2[2]]) - mean(event.b[trt.b==1 & marker2.b ==mval2[2]])
           
           )
      
    }
  
 if(conf.bands){
  boot.data <- replicate(bootstraps, one.boot.plot_disc( event, trt, marker1, marker2, rho,obp.boot.sample = boot.sample))
  mval1 <- sort(unique(marker1))
  mval2 <- sort(unique(marker2))
  
  row.names(boot.data) = c(paste("trteffect.1mkr", mval1[1], sep = ""), 
                           paste("trteffect.1mkr", mval1[2], sep = ""),
                           paste("trteffect.2mkr", mval2[1], sep = ""), 
                           paste("trteffect.2mkr", mval2[2], sep = ""))
  
  #horizontal
  if(substr(ci, 1,1) =="h") { warning("Horizontal CI bands are not allowed for treatment effect plots with a discrete marker. Vertical bands will be computed"); ci <- "vertical";}

  #vertical
  myconf.ints <- apply(boot.data, 1, quantile, probs = c(alpha/2, 1-alpha/2))

  ci = "vertical"
 }else{
   myconf.ints = NULL
 }
  

    ts1.curves <- trteffectPLOTcompare_gg_disc(x1=ts1, x2 = ts2, ci.bounds = myconf.ints, conf.bands = conf.bands, offset = offset, xlab=xlab, ylab = ylab, xlim=xlim, ylim = ylim, main = main,marker.names = marker.names,  annotate.plot = annotate.plot) 
  
    p <- ts1.curves[[1]]
  print(p)
    #if(is.null(xlim)) xlim = c(0,100)
    #legend(x=xlim[2]+diff(xlim)/15, y = quantile(ylim, prob = .5), legend = marker.names, lty = c(1,2), col = c("black", "black"), bty="n", cex = 1, xpd = TRUE, lwd=c(2,2))


   result <- list("plot" = p, 
                  "ci.bounds" = ts1.curves[[2]])
 


   invisible(result)
#par(old.par)

}
