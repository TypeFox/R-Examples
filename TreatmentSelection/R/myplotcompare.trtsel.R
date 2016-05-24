myplotcompare.trtsel <-
function(x, bootstraps = 500, alpha = .05,
                           ci   = "horizontal", marker.names = c("Marker 1", "Marker 2"),
                           fixeddeltas.y1 =  NULL, fixeddeltas.y2,  
                           xlab = NULL, 
                           ylab = NULL, 
                           xlim = NULL, 
                           ylim = NULL, 
                           main = NULL, offset = offset, conf.bands)
{


  ts1 <- x$trtsel1
  ts2 <- x$trtsel2

  if(substr(ci, 1, 4) =="hori") {
           ## trt selcurve plot with horizontal ci bands
           fix.ind = 2    # fix delta
           out.ind = 1    # vary F.marker
        }else if(substr(ci, 1, 4) =="vert"){
           ## trt sel curve plot with vertical ci bands
           fix.ind = 1    #fix F.event
           out.ind = 2    #vary delta
        }

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
  get.F <- ts1$functions$get.F

###
#old.par <- par(no.readonly = TRUE)


  if( length(fixeddeltas.y1) > 0 & bootstraps >0 ) {


    #bootstrapping done by fixing response and strapping marker 

   # if(is.null(fixed.values) & ci =="horizontal"){
   # fixeddeltas.y1 <- seq(from = min(delta.y1), to = max(delta.y1), length.out = 100)
   # fixeddeltas.y2 <- seq(from = min(delta.y2), to = max(delta.y2), length.out = 100)
   # }else if(is.null(fixed.values) & ci =="vertical"){
   # fixeddeltas.y1 <- 1:100/100
   # fixeddeltas.y2 <- 1:100/100
   
   # }else{
    #fixeddeltas.y1 <- fixed.values
    #fixeddeltas.y2 <- fixed.values
   # }
    if(link == "risks_provided"){
      
      provided_risk <- cbind(fittedrisk.t0.y1, 
                         fittedrisk.t1.y1, 
                         fittedrisk.t0.y2,
                         fittedrisk.t1.y2)
    }else{
      
      provided_risk = NULL
    }

    boot.dat <- replicate( bootstraps, one.boot.plot.compare(event, trt, marker1, marker2, ci, 
                                                              fixeddeltas.y1 = fixeddeltas.y1, fixeddeltas.y2 = fixeddeltas.y2,
                                                              rho = rho, study.design = study.design,  obp.boot.sample = boot.sample, obp.get.F = get.F, fix.ind, out.ind, link = link, 
                                                             provided_risk = provided_risk))
    

    if(length(fixeddeltas.y1)==1){
      bounds.delta.y1<- quantile(boot.dat[1,,], probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
      bounds.delta.y2<- quantile(boot.dat[2,,], probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)

      bounds.delta.y1 = t(t(bounds.delta.y1))
      bounds.delta.y2 = t(t(bounds.delta.y2))
    }else{
      bounds.delta.y1<- apply(boot.dat[1,,], 1, function(x, ...){quantile(unlist(x), ...)}, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
      bounds.delta.y2<- apply(boot.dat[2,,], 1, function(x, ...){quantile(unlist(x), ...)}, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
   

    }


    }else{
    bounds.delta.y1 <- NULL
    bounds.delta.y2 <- NULL
    }
    ## end bootstrap
    
    #set up the ylimit of the treatment effect plot
    if(is.null(ylim)){

    if(substr(ci, 1,1) %in% c("v", "V")){
    min.delta <- min(c(delta.y1, delta.y2, bounds.delta.y1, bounds.delta.y2), na.rm=TRUE)
    max.delta <- max(c(delta.y1, delta.y2, bounds.delta.y1, bounds.delta.y2), na.rm=TRUE)

    cen <- mean(c(min.delta, max.delta), na.rm=TRUE)
    }else{

    min.delta <- min(c(delta.y1, delta.y2), na.rm=TRUE)
    max.delta <- max(c(delta.y1, delta.y2), na.rm=TRUE)

    cen <- mean(c(min.delta, max.delta), na.rm=TRUE)
    
    }

    ran <- max.delta - min.delta
    ran <- ran*1.1
    ylim <- c(cen-ran/2, cen+ran/2)
    }
  


    ts1.curves <- trteffectPLOTcompare_gg(x1=ts1, x2 = ts2, ci = ci, ci.bounds = rbind(bounds.delta.y1, bounds.delta.y2), get.F = get.F, fixed.values = rbind(fixeddeltas.y1, fixeddeltas.y2+offset), conf.bands = conf.bands, rho=rho, xlab=xlab, ylab = ylab, xlim=xlim, ylim = ylim, main = main) 

    p <- ts1.curves[[1]]
    p <- p + scale_linetype_manual(name = "", breaks = c("1","2", "3", "4"), values = c(1, 2, 3, 4), labels = c(marker.names, "Mean", "Zero"))
    p <- p + scale_size_manual(name = "", breaks = c("1", "2", "3", "4"), values = c(1, 1, .5, .5), labels = c(marker.names, "Mean", "Zero"))
    
  print(p)
    #if(is.null(xlim)) xlim = c(0,100)
    #legend(x=xlim[2]+diff(xlim)/15, y = quantile(ylim, prob = .5), legend = marker.names, lty = c(1,2), col = c("black", "black"), bty="n", cex = 1, xpd = TRUE, lwd=c(2,2))

   
   if(bootstraps >0 & length(fixeddeltas.y1) > 0 ){
   conf.ints.y1 <- as.data.frame(cbind(fixeddeltas.y1, t(bounds.delta.y1)))
   names(conf.ints.y1) <- c("fixed", "lower", "upper")
   conf.ints.y2 <- as.data.frame(cbind(fixeddeltas.y2, t(bounds.delta.y2)))
   names(conf.ints.y2) <- c("fixed", "lower", "upper")


   result <- list("plot" = p, 
                  "trtsel1" = list( "conf.intervals" = conf.ints.y1), 
                  "trtsel2" = list( "conf.intervals" = conf.ints.y2))
 
   }else{

   result <- list("plot" = p)
 

   }

   invisible(result)
#par(old.par)

}
