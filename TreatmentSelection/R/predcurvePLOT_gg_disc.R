predcurvePLOT_gg_disc <-
function(x, ci, ci.bounds, get.F,  xlab, ylab, xlim, ylim, main, trt.names){ 

  
  fittedrisk.t0 <- x$derived.data$fittedrisk.t0
  fittedrisk.t1 <- x$derived.data$fittedrisk.t1
  marker <- x$derived.data$marker
  event <- x$derived.data$event
  trt <- x$derived.data$trt
  
  mydata <- data.frame(risk = fittedrisk.t0*(1-trt)+fittedrisk.t1*trt, trt = trt,  marker)

  mydata <- unique(mydata)
  mydata$lower <- rep(NA, nrow(mydata))
  mydata$upper <- rep(NA, nrow(mydata))
  mval = sort(unique(marker))
   
  #appease check 
  risk <- lower <- upper <- NULL

  if(!is.null(ci.bounds)){
    
  

  #order matters here!
  mydata[mydata$trt==0 & mydata$marker==mval[1], 4:5] <- ci.bounds[,1]
  mydata[mydata$trt==1 & mydata$marker==mval[1], 4:5] <- ci.bounds[,2]
  mydata[mydata$trt==0 & mydata$marker==mval[2], 4:5] <- ci.bounds[,3]
  mydata[mydata$trt==1 & mydata$marker==mval[2], 4:5] <- ci.bounds[,4]

  
    p <- ggplot(mydata, aes(x = factor(marker), y = risk, ymin = lower, ymax = upper, group = factor(trt), shape = factor(trt), linetype = factor(trt) ))
    p <- p + geom_errorbar(size = 1, width = .1) + geom_point(size = 4)
  
  
  }else{

      p <- ggplot(mydata, aes(x = factor(marker), y = risk, group = factor(trt), shape = factor(trt), linetype = factor(trt) ))
      p <- p + geom_point(size = 4)


  }

  if(is.null(xlab)) xlab <- "marker value"
  if(is.null(ylab)) ylab <- "risk given marker"
  if(is.null(xlim)) xlim <- c(0,100)
  if(is.null(ylim)) ylim <- c(0,1)
  if(is.null(main)) main <- "Risk curves by treatment"

  #add x/y labels and main
  p <- p + xlab(xlab) + ylab(ylab) + ylim(ylim[1], ylim[2]) + ggtitle(main) 
  #change the names for the legend
  p <- p + scale_shape_discrete(breaks=factor(c(1,0)), labels = trt.names) +
    scale_linetype_discrete(breaks=factor(c(1,0)), labels = trt.names)+
    theme(legend.title = element_blank(),  text = element_text(size=14))
  #legend.text = element_text(size = 16))
  p <- p + scale_x_discrete(labels = c(paste(mval[1], "\n(", round(mean(marker==mval[1])*100, 1),"%)", sep = ""), 
                                       paste(mval[2], "\n(", round(mean(marker==mval[2])*100, 1),"%)", sep = "")))
  print(p)
  return(list(p, mydata))
}
