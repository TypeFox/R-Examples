wc.avg <-
function(WC, which.avg = 'wp', 
         show.siglvl = T, siglvl = c(0.05, 0.1), sigcol = c('red','blue'), sigpch = 20, 
         label.avg.axis = T, averagelab = NULL, label.period.axis = T, periodlab = NULL, 
         show.legend = T, legend.coords='topright',
         main = NULL, lwd = 0.5, 
         verbose = F){
                   
                   
  ##########################################       
 
  if(verbose == T){
       out <- function(...){ cat(...) }
    }
  else{
       out <- function(...) { }
    }  
           
  ##########################################           
           

  if (class(WC) == 'analyze.wavelet') { stop("Your object class is 'analyze.wavelet' --- please use wt.avg!") }                 
                   
  if (which.avg == 'wc') {
      W.avg = WC$Coherence.avg
      W.avg.pval = WC$Coherence.avg.pval
  }  
  if (which.avg == 'wp') {
      W.avg = WC$Power.xy.avg
      W.avg.pval = WC$Power.xy.avg.pval
  }  
  
  out(paste("Averages across time of '", which.avg, "' will be plotted...", sep=''), '\n') 

  plot(W.avg, log2(WC$Period), lwd = lwd, type = "l", axes = FALSE, ylab = "", xlab = '', yaxs = 'i', main = main)
  
  if ((show.siglvl == T) & (is.null(W.avg.pval) == F)) { 
     
     P.dat = data.frame(Pvalue = W.avg.pval, Log.period = log2(WC$Period), Average = W.avg)
     
     n.sig = length(siglvl)
     if (length(sigpch) == 1) {sigpch = rep(sigpch,n.sig)}
     if (n.sig != length(sigpch)) {sigpch = rep(20,n.sig)}
     if (n.sig != length(sigcol)) {sigcol = 1:n.sig}
       
     siglvl.order = order(siglvl,decreasing=T)
     sig.params = data.frame(siglvl = siglvl[siglvl.order], sigcol=sigcol[siglvl.order], sigpch = sigpch[siglvl.order])
     for (i in (1:length(siglvl))) {
         with(P.dat[P.dat$Pvalue < sig.params$siglvl[i],], points(Average, Log.period, pch = sig.params$sigpch[i], col = as.character(sig.params$sigcol[i])))
     } 
      
     
     if (show.legend == T) {
         legend(legend.coords, legend=siglvl, pch=sigpch, col=sigcol, horiz=F, box.lwd=0.25)
     }
  }
  
  box(lwd = 0.25)
  
  
  # label average axis ?  
  if (label.avg.axis == T) {
  
    if (is.null(averagelab)) {
        if (which.avg=='wc') {averagelab='average coherence'}
        if (which.avg=='wp') {averagelab='average cross-wavelet power'}
    }    
  
      
    A.1 = axis(1, lwd = .25, labels=NA, tck = 0.0)
    mtext(A.1, side = 1, at = A.1, line = 0.5)
    mtext(averagelab, side = 1, line = 2, font = 1)  
     
  }

  # label period axis ?  
  if (label.period.axis == T) {
  
    if (is.null(periodlab)) {periodlab='period'}
    
    period.tick <- unique(trunc(WC$axis.2))
    period.tick[period.tick<log2(WC$Period[1])] = NA
    period.tick = na.omit(period.tick)
    period.tick.label <- 2^(period.tick)   
    axis(2, lwd = .25, at = period.tick, labels = NA, tck=0.02, tcl=0.5)
    axis(4, lwd = .25, at = period.tick, labels = NA, tck=0.02, tcl=0.5)
    mtext(period.tick.label, side = 2, at = period.tick, las = 1, line = 0.5)

    mtext(periodlab, side = 2, line = 2.5,font = 1)
  
  }
  
}
