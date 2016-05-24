analyze.coherency <-
function(my.data, my.pair = c(1,2), 
         loess.span = 0.75, dt = 1, dj=1/20, 
         lowerPeriod = 2*dt, upperPeriod = floor(nrow(my.data)/3)*dt,
         window.type.t=1, window.type.s=1, window.size.t=5, window.size.s=1/4,         
         make.pval = T, method = "white.noise", params = NULL, n.sim = 100,
         verbose = T) {
                             
  if(verbose == T){
     out <- function(...){ cat(...) }
  }
  else{
     out <- function(...) { }
  }    

###################################################################################################
## The following function smoothes the series in a data frame.
## Input: a data frame with dates as row names
## Output: a data frame with the same row and column names, but with smoothed series
###################################################################################################

  loess.data.frame = function(x, loess.span)  {
    x.smoothed = x
    for (i in 1:ncol(x))  {
      day.index = 1:nrow(x)
      my.loess.x = loess(x[, i] ~ day.index, span = loess.span)
      # smoothed series = fitted values:
      x.loess = as.numeric(predict(my.loess.x, data.frame(x = 1:nrow(x))))
      x.smoothed[, i] = x.loess
    }
    return(x.smoothed)
  }

###################################################################################################
## Select the pair of time series to be analyzed
###################################################################################################

  xy = my.data[,my.pair]
  my.pair = colnames(xy)
  
  if (length(my.pair) != 2) { stop('Please select two series for analysis!\n') }
  if (is.element('date', my.pair)) { stop('Please review your selection of series!\n') }

###################################################################################################
## Smooth the data (if requested)
###################################################################################################

  if (loess.span != 0) {
     out("Smoothing the time series...\n")
     xy.trend = loess.data.frame(xy, loess.span)
     xy = xy-xy.trend
     xy = cbind(xy, xy.trend)
     colnames(xy) = c(my.pair, paste(my.pair,'.trend',sep=''))    
  }
  
###################################################################################################
## Add date column if available
###################################################################################################  
  
  if (is.element('date',names(my.data))) {xy = cbind(date=my.data$date, xy)}  
  

###################################################################################################
## Start the analysis of wavelet coherency
###################################################################################################

  out("Starting wavelet transformation and coherency computation...\n")
  if (make.pval == T) { out("... and simulations... \n") }
  my.wc = wc(x=xy[[my.pair[1]]], y=xy[[my.pair[2]]], start = 1, 
             dt = dt, dj = dj, 
             lowerPeriod = lowerPeriod, upperPeriod = upperPeriod,
             window.type.t=window.type.t, window.type.s=window.type.s, window.size.t=window.size.t, window.size.s=window.size.s,
             make.pval = make.pval, method = method, params = params, n.sim = n.sim)
              
##################################################################################################
## Compute the ridges
################################################################################################## 
 
  Ridge.xy = ridge(my.wc$Power.xy)
  Ridge.co = ridge(my.wc$Coherence)

  Ridge.x = ridge(my.wc$Power.x)
  Ridge.y = ridge(my.wc$Power.y)
  
##################################################################################################
## Prepare the output
##################################################################################################  


  output = list(series = xy, loess.span = loess.span, dt = dt, dj = dj,
                Wave.xy = my.wc$Wave.xy, Angle = my.wc$Angle,
                sWave.xy = my.wc$sWave.xy, sAngle = my.wc$sAngle,     
                Power.xy = my.wc$Power.xy, Power.xy.avg = my.wc$Power.xy.avg, 
                Power.xy.pval = my.wc$Power.xy.pval, Power.xy.avg.pval = my.wc$Power.xy.avg.pval, 
                Coherency = my.wc$Coherency,
                Coherence = my.wc$Coherence, Coherence.avg = my.wc$Coherence.avg,
                Coherence.pval = my.wc$Coherence.pval, Coherence.avg.pval = my.wc$Coherence.avg.pval,
                Wave.x = my.wc$Wave.x, Wave.y = my.wc$Wave.y,
                Phase.x = my.wc$Phase.x, Phase.y = my.wc$Phase.y,
                Ampl.x = my.wc$Ampl.x, Ampl.y = my.wc$Ampl.y,             
                Power.x = my.wc$Power.x, Power.y = my.wc$Power.y, 
                Power.x.avg = my.wc$Power.x.avg, Power.y.avg = my.wc$Power.y.avg,
                Power.x.pval = my.wc$Power.x.pval, Power.y.pval = my.wc$Power.y.pval, 
                Power.x.avg.pval = my.wc$Power.x.avg.pval, Power.y.avg.pval = my.wc$Power.y.avg.pval,
                sPower.x = my.wc$sPower.x, sPower.y = my.wc$sPower.y,
                Ridge.xy = Ridge.xy,
                Ridge.co = Ridge.co,
                Ridge.x = Ridge.x, Ridge.y = Ridge.y,
                Period = my.wc$Period, Scale = my.wc$Scale,
                nc = my.wc$nc, nr = my.wc$nr,
                coi.1 = my.wc$coi.1, coi.2 = my.wc$coi.2, 
                axis.1 = my.wc$axis.1, axis.2 = my.wc$axis.2)
                            
  
  class(output) = "analyze.coherency"

  out("Class attributes are accessible through following names:\n")
  out(names(output), "\n")

  return(invisible(output))
}
