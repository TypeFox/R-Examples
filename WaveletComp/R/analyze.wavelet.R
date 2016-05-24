analyze.wavelet <-
function(my.data, my.series = 1, loess.span = 0.75, dt = 1, dj=1/20, 
         lowerPeriod = 2*dt, upperPeriod = floor(nrow(my.data)/3)*dt,
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
## Output: a data frame with the same row and column names, with smoothed series
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
## Select the time series to be analyzed
###################################################################################################

   if (is.numeric(my.series)) { 
          my.series = names(my.data)[my.series] 
      }
      
   if (length(my.series) != 1) { stop('Please select (only) one series for analysis!\n') }   
   if (is.element('date', my.series)) { stop('Please review your selection of series!\n') }  
   
   ind = which( names(my.data) == my.series )
   x = data.frame(my.data[,ind])   
   colnames(x) = my.series
   rownames(x) = rownames(my.data)

###################################################################################################
## Smooth the data (if requested)
###################################################################################################

  if (loess.span != 0) {
     out("Smoothing the time series...\n")
     x.trend = loess.data.frame(x, loess.span)
     x = x-x.trend
     x = cbind(x, x.trend)
     colnames(x) = c(my.series, paste(my.series,'.trend',sep=''))
  }
  
###################################################################################################
## Add date column if available
###################################################################################################  
  
  if (is.element('date',names(my.data))) {x = cbind(date=my.data$date, x)}  
  
###################################################################################################
## Start the analysis of wavelets
###################################################################################################

  out("Starting wavelet transformation...\n")
  if (make.pval == T) { out("... and simulations... \n") }
  my.wt = wt(x=x[[my.series]], start = 1, 
             dt = dt, dj = dj, 
             lowerPeriod = lowerPeriod, upperPeriod = upperPeriod,
             make.pval = make.pval,
             method = method,
             params = params,
             n.sim = n.sim, save.sim = F)
              
##################################################################################################
## Compute the power ridge
##################################################################################################

  Ridge = ridge(my.wt$Power)
  
##################################################################################################  
## Prepare the output  
##################################################################################################


  output <- list(series = x, loess.span = loess.span, dt = dt, dj = dj,
                 Wave = my.wt$Wave, Phase = my.wt$Phase, Ampl = my.wt$Ampl,
                 Power = my.wt$Power, Power.avg = my.wt$Power.avg,
                 Power.pval = my.wt$Power.pval, Power.avg.pval = my.wt$Power.avg.pval,  
                 Ridge = Ridge,     
                 Period = my.wt$Period, Scale = my.wt$Scale,                      
                 nc = my.wt$nc, nr = my.wt$nr,      
                 coi.1 = my.wt$coi.1, coi.2 = my.wt$coi.2,
                 axis.1 = my.wt$axis.1, axis.2 = my.wt$axis.2               
                )

  
  class(output) = "analyze.wavelet"

  out("Class attributes are accessible through following names:\n")
  out(names(output), "\n")

  return(invisible(output))

}
