wt.sel.phases <-
function(WT, 
         sel.period = NULL, sel.lower = NULL, sel.upper = NULL,   
         only.coi = F,
         only.sig = T, siglvl=0.05, 
         show.avg.phase = F, phase.avg.col = 'black',
         label.time.axis = T, 
         show.date = F, date.format = NULL, 
         timelab = NULL, 
         label.phase.axis = T, phaselab = NULL,
         main = NULL, sub = NULL, 
         verbose = F) {
         
   ##########################################      
                     
   if(verbose == T){
       out <- function(...){ cat(...) }
    }
    else{
       out <- function(...) { }
    }   
    
   nc = WT$nc    
   nr = WT$nr

   series.data = WT$series
   
  
   #######################################################################################
   ## Select periods
   #######################################################################################
   
   Period = WT$Period

   # sel.period available?
   if (length(sel.period) != 0) {
   
       # nearest available period:
       sel.rnum = which(abs(Period-sel.period) == min(abs(Period-sel.period)))
       
   }    
   
   # in case, sel.period is not available, refer to sel.upper, sel.lower
   if (length(sel.period) == 0) {
   
       # in case, sel.lower (sel.upper) is not available, use minimum (maximum) period 
       if (length(sel.lower) == 0) {sel.lower = min(Period)}
       if (length(sel.upper) == 0) {sel.upper = max(Period)}
             
       # in case, sel.lower > sel.upper, interchange
       if (sel.lower > sel.upper) { 
           sel.lower.h = sel.lower
           sel.lower = sel.upper
           sel.upper = sel.lower.h
       }
      
       sel.rnum = which(((Period >= sel.lower) & (Period <= sel.upper)))
       
   }    
       
   # selected band / range of periods   
   sel.period.band  = Period[sel.rnum]  
   sel.period.range = as.character(round(sel.period.band, 3))
    
   #######################################################################################
   ## Select corresponding phases and phase differences
   ####################################################################################### 
   
   Phase = WT$Phase 
    
   # only cone of influence for reconstruction? 
   if (only.coi == T) {
       for (i in (1:nc)) {
            for (s.ind in seq_len(nr)) {
                 if (WT$Scale[s.ind] > 2^WT$coi.2[i]) {Phase[s.ind, i] = NA}
            }
       } 
   }  
    
   # use significant parts only?
   if (only.sig == T) {
       if (!is.null(WT$Power.pval)) { 
       
          ind.sig = (WT$Power.pval < siglvl)
          Phase[ind.sig == 0] = NA
          
       }      
   }  
   if (only.sig == F) { siglvl = NA }
  
   Phase = Phase[sel.rnum,]
   phase.col = 1:length(sel.rnum)
    
   # in case, several periods have been selected, compute column averages if desired
   if (length(sel.rnum) > 1) {
       sel.period.range = paste(round(range(sel.period.band), 3), collapse=' - ')
       if (show.avg.phase == T) {
           Phase = colMeans(Phase, na.rm=T)  
           phase.col = phase.avg.col
       }  
       if (show.avg.phase == F) {
           Phase = t(Phase)        
       }    
   }
   
   #######################################################################################
   ## start plotting
   #######################################################################################
   
   if (is.null(sub)) { sub=paste('selected period: ', sel.period.range, sep='') }
      
   matplot(WT$axis.1, Phase, type='l', lty=1, col=phase.col, 
           xaxs = 'i', xaxt='n', ylim=c(-pi,pi), yaxt='n', 
           xlab='', ylab='', 
           main = main,
           sub=sub)
   
   
   # label phase axis ?  
   if (label.phase.axis == T) {
   
       if (is.null(phaselab)) {
           phaselab='phase'
           if (show.avg.phase == T) {
               phaselab='average phase'
           }
       }    
       
       axis(2, lwd = .25, at=seq(-pi, pi, pi/2), labels=NA, tck=0.02, tcl=0.5)
       axis(4, lwd = .25, at=seq(-pi, pi, pi/2), labels=NA, tck=0.02, tcl=0.5)
       mtext(c(expression(-pi),expression(-pi/2),0,expression(pi/2),expression(pi)), at=seq(-pi, pi, pi/2), side = 2, line = 0.5, las=1, font = 1)
       mtext(phaselab, side = 2, line = 2.5, font = 1)
       
   }
   
   my.date = NULL
   
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
               plot(my.date, rep(0,WT$nc), ylim=c(-pi,pi), type="n", xaxs = "i", yaxt='n', xlab="", ylab="",
                    lwd = .25, tck=0.02, tcl=0.5)
               mtext(timelab, side = 1, line = 2.5) 
               
      }    
   }  
   
  
  ####################################
  ## Output
  ####################################
    
  output <- list(Period = sel.period.band,
                 Phase = Phase,
                 only.coi = only.coi,
                 only.sig = only.sig, siglvl = siglvl, 
                 date = my.date,
                 time.axis = WT$axis.1
                 )
                 
  class(output) = "sel.phases"

  out("Class attributes are accessible through following names:\n")
  out(names(output), "\n")             
     
  return(invisible(output))  

}
