reconstruct <-
function(WT, my.series = 1, lvl = 0,
         only.coi = F,
         only.sig = T, siglvl=0.05,         
         only.ridge = F, 
         sel.period = NULL, sel.lower = NULL, sel.upper = NULL,   
         rescale = T,
         plot.waves = F, 
         plot.rec = T, lty = 1, lwd = 1, col = 1:2, ylim = NULL,
         show.legend = T, legend.coords = 'topleft', legend.horiz=F,
         legend.text = NULL,
         label.time.axis = T, show.date = F, date.format = NULL, timelab = NULL,
         main.waves = NULL, main.rec = NULL, main = NULL,
         verbose = T) {
                         
    if(verbose == T){
       out <- function(...){ cat(...) }
    }
    else{
       out <- function(...) { }
    } 
    
    lwd.axis = 0.25
    
    series.data = WT$series
    
    ####################################
    ## Identify the scenario
    ####################################
    
    if (class(WT) == 'analyze.wavelet') {
    
        out("Your input object class is 'analyze.wavelet'...\n") 
        
        my.series = ifelse(names(series.data)[1] == 'date', names(series.data)[2], names(series.data)[1]) 
        
        Wave  = WT$Wave
        Power = WT$Power
        Power.pval = WT$Power.pval
        Ridge = WT$Ridge     
        
    } 
    if (class(WT) == 'analyze.coherency') {   
    
        out("Your input object class is 'analyze.coherency'...\n")   
        
        if (is.numeric(my.series)) { 
            if (!is.element(my.series,c(1,2))) { stop("Please choose either series number 1 or 2!") }
            my.series = ifelse(names(series.data)[1] == 'date', names(series.data)[my.series+1], names(series.data)[my.series])  
        }
        
        ind = which( names(series.data) == my.series ) 
        which.series.num = ifelse(names(series.data)[1] == 'date', ind-1, ind)
        if (!is.element(which.series.num, c(1,2))) { stop("Your series name is not available, please check!") }
      
        if (which.series.num == 1) {
            Wave  = WT$Wave.x
            Power = WT$Power.x
            Power.pval = WT$Power.x.pval
            Ridge = WT$Ridge.x
        }
        if (which.series.num == 2) {
            Wave  = WT$Wave.y
            Power = WT$Power.y
            Power.pval = WT$Power.y.pval
            Ridge = WT$Ridge.y
        }      
        
    }   
           
    out(paste("Your time series '", my.series, "' will be reconstructed...", sep=''), '\n')   
   

    ####################################
    ## Prepare reconstruction components
    ####################################
    
    out("Starting the reconstruction process...\n")
    
    nc = WT$nc    
    nr = WT$nr

    dt = WT$dt
    dj = WT$dj
    
    Scale = WT$Scale
    Period = WT$Period   
       
    loess.span = WT$loess.span
   
   
    rec.waves = matrix(0, nrow=nr, ncol=nc)
    for (s.ind in seq_len(nr)) {
         rec.waves[s.ind,] = (Re(Wave[s.ind,])/sqrt(Scale[s.ind]))*dj*sqrt(dt)/(pi^(-1/4)*0.776)
    }
    
    # select minimum level?
    rec.waves = rec.waves * (Power >= lvl)
    comment.lvl = paste('minimum power level: ',lvl, sep='')
    
    # select ridge?
    if (only.ridge == T) {      
        rec.waves = rec.waves * Ridge  
        rec.waves[Ridge == 0] = NA
    }  
    comment.ridge = paste('only ridge: ', only.ridge, sep='')
    

    # use significant parts only?
    if (only.sig == T) {
        if (!is.null(Power.pval)) { 
            rec.waves = rec.waves * (Power.pval < siglvl)
            rec.waves[Power.pval >= siglvl] = NA
        }      
    }  
    if (only.sig == F) { siglvl = NA }
    comment.sig = paste('significance level: ',siglvl,sep='')
    
    # only cone of influence for reconstruction? 
    if (only.coi == T) {
        for (i in (1:nc)) {
             for (s.ind in seq_len(nr)) {
                  if (Scale[s.ind] > 2^WT$coi.2[i]) {rec.waves[s.ind, i] = NA}
             }
        } 
    }  
    comment.coi = paste('only coi: ', only.coi, sep='') 
    
    # so far
    rnum.used = which(rowSums(rec.waves, na.rm=T)!=0)
    comment.periods = 'period: all relevant'
    
       
    # sel.period available?
    if (length(sel.period) != 0) {
    
        sel.rnum = numeric()
    
        # nearest available period:
        for (i in (1:length(sel.period))) {
             sel.rnum = union(sel.rnum, which(abs(Period-sel.period[i]) == min(abs(Period-sel.period[i])))) 
        }  
        
        rec.waves = rec.waves[sel.rnum,]
        comment.periods = paste('period: ',paste(as.character(round(Period[sel.rnum], 1)), collapse=', '), sep='')
           
        if (length(sel.rnum) == 1) {
            rec.waves = t(rec.waves)
        } 
        
        rnum.used = intersect(rnum.used, sel.rnum)
    } 
    
    # in case, sel.period is not available, refer to sel.upper, sel.lower
    if (length(sel.period) == 0 & ((length(sel.lower) != 0) | (length(sel.upper) != 0))) {   
    
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
       rec.waves = rec.waves[sel.rnum,]
       
       # selected band / range of periods   
       sel.period.band  = Period[sel.rnum]  
       sel.period.range = as.character(round(sel.period.band, 1))
       if (length(sel.rnum) > 1) {
           sel.period.range = paste(round(range(sel.period.band), 1), collapse=' - ')
       }   
       if (length(sel.rnum) == 1) {rec.waves = t(rec.waves)}
       comment.periods = paste('period: ', sel.period.range, sep='') 
       
       rnum.used = intersect(rnum.used, sel.rnum)
       
    }
         
        
    
    ####################################
    ## Compute reconstructed series
    ####################################   
       
    # reconstructed time series
    x.r  = colSums(rec.waves, na.rm=T)
    
    # retrieve original time series
    x    = series.data[[my.series]]
    
    # rescale the reconstructed time series?
    if (rescale == T) {
        x.r  = (x.r-mean(x.r))*sd(x)/sd(x.r) + mean(x) 
    }   
    
    
    ####################################
    ## Plottings
    ####################################
    
    ## plot of reconstruction waves
        
    if (plot.waves == T) {
    
        out("Reconstruction waves are being plotted...\n")
        
        if (is.null(main) == F) {main.waves = main}
        
        range.rec.waves = range(rec.waves, na.rm =T)
        matplot(WT$axis.1, t(rec.waves), type = 'l', ylim = range.rec.waves,
                main = main.waves, sub = paste(comment.lvl, ', ', comment.sig,', ',comment.coi,', ',comment.ridge,', ', comment.periods,  sep=''),
                xaxs = 'i', xaxt = 'n', 
                xlab = '', ylab = '')

       # label time axis ?   
       if (label.time.axis == T) {
  
           if (is.null(timelab)) {timelab='time'}
      
           if (show.date == F) {
               A.1 = axis(1, lwd = lwd.axis, labels=NA, tck = 0.0)
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
               plot(my.date, seq(range.rec.waves[1], range.rec.waves[2],length.out=WT$nc), type="n", ylim = range.rec.waves,
               xaxs = "i", yaxt='n', xlab="", ylab="",
                    lwd = lwd.axis, tck=0.02, tcl=0.5)
               mtext(timelab, side = 1, line = 2.5)  
               
           }    
       }                   
                
       
    } #end if (plot.waves == T)
    
    x.rec = cbind(series.data, x.r=x.r)
    my.rec.series = paste(my.series,'.r',sep='')
    colnames(x.rec) = c(colnames(series.data),my.rec.series)
    rownames(x.rec) = rownames(series.data)
        
    if (plot.waves & plot.rec) { 
        par(ask = T)
    }
    
    
    
    ## plot of reconstructed series
    
    if (plot.rec == T) {
    
          
       out("Original (detrended) and reconstructed series are being plotted...\n")
       
       if (is.null(main) == F) { main.rec = main }
       
       range.rec = range(x.rec[,c(my.series,my.rec.series)], na.rm=T)
       if (is.null(ylim)) {ylim = range.rec} 
       matplot(WT$axis.1, x.rec[,c(my.series,my.rec.series)], type = 'l', ylim = ylim, lty = lty, lwd = lwd, col = col,
               main = main.rec, sub = paste(comment.lvl, ', ', comment.sig,', ',comment.coi,', ',comment.ridge,', ', comment.periods,  sep=''),
               xaxs = 'i', xaxt = 'n', 
               xlab = '', ylab = '')
               
       par(ask = F)
       
       if (show.legend == T) {      
       
           if (is.null(legend.text)) { 
               legend.text = c(paste('original', ifelse(loess.span!=0, paste(' (detrended, span: ', loess.span, ')', sep=''),''), sep=''), 'reconstructed')
               }
               
           legend(legend.coords, horiz=legend.horiz, legend = legend.text, lty = lty, lwd = lwd, col = col)        
       }
       
       
       # label time axis ?   
       if (label.time.axis == T) {
  
           if (is.null(timelab)) {timelab='time'}
      
           if (show.date == F) {
               A.1 = axis(1, lwd = lwd.axis, labels=NA, tck = 0.0)
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
               plot(my.date, seq(ylim[1], ylim[2], length.out=WT$nc), type="n", 
               ylim = ylim, xaxs = "i", yaxt='n', xlab="", ylab="",
                    lwd = lwd.axis, tck=0.02, tcl=0.5)
               mtext(timelab, side = 1, line = 2.5)  
               
           }    
       }                             
              
  
    }
    
    
    
    ####################################
    ## Output
    ####################################
    
    output <- list(series = x.rec,
                 rec.waves = rec.waves,
                 loess.span = loess.span,
                 lvl = lvl,
                 only.coi = only.coi,
                 only.sig = only.sig, siglvl = siglvl, 
                 only.ridge = only.ridge,                
                 rnum.used = rnum.used,   
                 rescale = rescale,
                 dt = dt, dj = dj,
                 Period = Period, Scale = Scale,
                 nc = nc, nr = nr, 
                 axis.1 = WT$axis.1,
                 axis.2 = WT$axis.2
                 )
                 
    class(output) = "reconstruct"

    out("Class attributes are accessible through following names:\n")
    out(names(output), "\n")             
     
    return(invisible(output))
       
      
}
