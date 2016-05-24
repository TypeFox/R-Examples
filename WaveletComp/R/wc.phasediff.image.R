wc.phasediff.image <-
function(WC, use.sAngle = F, 
         plot.coi = T, 
         plot.contour = T, which.contour = "wp", siglvl = 0.1, col.contour = "white",
         n.levels=100, color.palette = "rainbow(n.levels, start=0, end=.7)",
         useRaster = T, max.contour.segments = 250000,
         plot.legend = T,
         legend.params = list(width=1.2, shrink=0.9, mar=5.1, n.ticks=6, label.digits=2, label.format="f", lab=NULL, lab.line=2.5),
         label.time.axis = T, show.date = F, date.format = NULL, timelab = NULL, 
         label.period.axis = T, periodlab = NULL,                           
         main = NULL,
         graphics.reset = T,
         verbose = F){
 
  ##########################################
  
  if(verbose == T){
       out <- function(...){ cat(...) }
    }
    else{
       out <- function(...) { }
    }  
    
  default.options = options()   
 
  options(max.contour.segments = as.integer(max.contour.segments))
  
  ##########################################
  
  if (class(WC) == 'analyze.wavelet') { stop("Your object class is 'analyze.wavelet' --- please use wt.phase.image!") }
  
  series.data = WC$series
  
  ##########################################
                          
  axis.1 <- WC$axis.1
  axis.2 <- WC$axis.2 
  
  # which version of Angle?
  if (use.sAngle == T) {
      Angle <- WC$sAngle
  }
  if (use.sAngle == F) {
      Angle <- WC$Angle
  }

  phasediff.levels = seq(-pi,pi,length.out=n.levels+1)
  key.cols = rev(eval(parse(text=color.palette)))
  
   # legend parameters  
  
  if (is.null(legend.params$width))         legend.params$width = 1.2
  if (is.null(legend.params$shrink))        legend.params$shrink = 0.9
  if (is.null(legend.params$mar)) legend.params$mar = ifelse(is.null(legend.params$lab), 5.1, 6.1)   
  if (is.null(legend.params$n.ticks))       legend.params$n.ticks = 6
  if (is.null(legend.params$label.digits))  legend.params$label.digits = 1
  if (is.null(legend.params$label.format))  legend.params$label.format = "f"
  if (is.null(legend.params$lab.line))      legend.params$lab.line = 2.5
  
  #######################################################################################
  ## start plotting
  #######################################################################################
  
  op = par(no.readonly = TRUE)
    
  image.plt  = par()$plt
  legend.plt = NULL
  
  if (plot.legend == T) {
  
      # construct plot regions for image and legend
   
      legend.plt = par()$plt
      
      char.size = par()$cin[1]/par()$din[1]
    
      hoffset       = char.size * par()$mar[4]
      legend.width  = char.size * legend.params$width
      legend.mar    = char.size * legend.params$mar
      
      legend.plt[2] = 1 - legend.mar
      legend.plt[1] = legend.plt[2] - legend.width
  
      vmar = (legend.plt[4] - legend.plt[3]) * ((1 - legend.params$shrink)/2)
  
      legend.plt[4] = legend.plt[4] - vmar
      legend.plt[3] = legend.plt[3] + vmar

      image.plt[2] = min(image.plt[2], legend.plt[1] - hoffset)
      
      # plot legend first
               
      par(plt = legend.plt)

      key.marks  = round(seq(from = 0, to = 1, length.out=legend.params$n.ticks)*n.levels)
      key.labels = formatC(as.numeric(phasediff.levels), digits = legend.params$label.digits, format = legend.params$label.format)[key.marks+1]
      
      image(1, seq(from = 0, to = n.levels), matrix(phasediff.levels, nrow=1), col = key.cols, breaks = phasediff.levels, useRaster=T, xaxt='n', yaxt='n', xlab='', ylab='')
      axis(4, lwd=.25, at=key.marks, labels=NA, tck=0.02, tcl=(par()$usr[2]-par()$usr[1])*legend.params$width-0.04)
      mtext(key.labels, side = 4, at = key.marks, line = 0.5, las=2)
      text(x = par()$usr[2] + (1.5+legend.params$lab.line)*par()$cxy[1], y=n.levels/2, labels=legend.params$lab, xpd=NA, srt = 270)
  
      box(lwd = .25)    
    
      par(new=TRUE, plt = image.plt)  
      
  }    
  
  #######################################################################################
  ## plot image of phase differences
  #######################################################################################
  
  image(axis.1, axis.2, t(Angle), col = key.cols, breaks = phasediff.levels, useRaster = useRaster,
        ylab = "", xlab = '', axes = FALSE, main = main)                
       
  # plot contour lines?     
  if  (plot.contour == T) {
       
       if (which.contour == 'wc') {
           W.pval = WC$Coherence.pval
       }  
       if (which.contour == 'wp') {
           W.pval = WC$Power.xy.pval
       } 
       
       if (is.null(W.pval) == F) {      
           contour(axis.1, axis.2, t(W.pval) < siglvl, levels = 1, lwd = 2, # lwd = 0.5 before
                   add = TRUE, col = col.contour, drawlabels = FALSE)
       } 
  }     

  # plot cone of influence?
  if (plot.coi == T) {
      polygon(WC$coi.1, WC$coi.1, border = NA, col = rgb(1, 1, 1, 0.5))
  }  
  
  box(lwd = .25)
  
  
  # label period axis ?  
  if (label.period.axis == T) {
  
      if (is.null(periodlab)) {periodlab='period'}
    
      period.tick <- unique(trunc(axis.2))
      period.tick[period.tick<log2(WC$Period[1])] = NA
      period.tick = na.omit(period.tick)
      period.tick.label <- 2^(period.tick)   
      axis(2, lwd = .25, at = period.tick, labels = NA, tck=0.02, tcl=0.5)
      axis(4, lwd = .25, at = period.tick, labels = NA, tck=0.02, tcl=0.5)
      mtext(period.tick.label, side = 2, at = period.tick, las = 1, line = 0.5)

      mtext(periodlab, side = 2, line = 2.5,font = 1)
          
  }
  
  # label time axis ?   
  if (label.time.axis == T) {
  
      if (is.null(timelab)) {timelab='time'}
      
      if (show.date == F) {
               A.1 = axis(1, lwd = .25, labels=NA, tck = 0.0)
               mtext(A.1, side = 1, at = A.1, line = 0.5)
               mtext(timelab, side = 1, line = 2)
      }
      if (show.date == T) {  
      
               if (is.element('date',names(series.data))) { my.date = series.data$date }  
               else { my.date = rownames(series.data) }  
      
               if (is.null(date.format)) { my.date = as.Date(my.date) }
               else { my.date = as.POSIXct(my.date, format=date.format) }
               par(new=TRUE)
               #empty plot, but calendar
               plot(my.date, seq(min(axis.2),max(axis.2), length.out=WC$nc), type="n", xaxs = "i", yaxs ='i', yaxt='n', xlab="", ylab="",
                    lwd = .25, tck=0.02, tcl=0.5)
               mtext(timelab, side = 1, line = 2.5)      
      }    
  }  
    
   
  #######################################################################################
  ## apropos graphical parameters
  #######################################################################################
  
  # reset contour line options
  options(default.options)
  
  # reset graphical parameters?
  if (graphics.reset == T) {
      par(op)
  }     
  
  # output of graphical parameters
  
  output = list(op = op, image.plt = image.plt, legend.plt=legend.plt)
  class(output) = "graphical parameters"

  out("Class attributes are accessible through following names:\n")
  out(names(output), "\n")    
    
  return(invisible(output))   
  
}
