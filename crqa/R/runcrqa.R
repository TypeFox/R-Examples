## GPL > 3: written by Moreno I. Coco (moreno.cocoi@gmail.com) 
#  Translating some Matlab codes by Rick Dale (rdale@ucmerced.edu)

## mother function to compute cross-recurrence between 
## two time series (numerical and categorical)

## arguments: 

## Two ways (type: 1, 2) are implemented to compute recurrence:

## type = 1; a computationally fast solution that only calculates recurrence of the diagonal
## diagonal recurrence can be computed in two methods (profile, window)
## method = "profile": compute the recurrence profile over the all time series for different lags
## method = "window": compute recurrence over time by sliding a window
## of a given size, calculating recurrence for a maximum
## number of lags smaller than half the size of the window.
## I would have to explain the above better.

## pro: very fast, ideal if the interest of the researcher
## is on properties of the diagonal profile
## cons: it does not return the full spectrum of measures
## from the recurrence plots

## type = 2; a recurrence matrix is calculated and several measures of
## recurrence: meanline, maxline, det, entropy ... are extracted from the matrix
## again recurrence can be computed both on the whole profile
## and based on the window.

## pro:  rich output about the characteristics of the recurrence profile.
## cons: very slow, and computationally intensive.

## par: list of arguments to be passed according to type
## of analysis and the method used. 

## For type 1:

## with method = "profile"
## ws = the width -/+ timestamps to use to lag the series

## with method = " window"  
## step = the sampling jumps over which windows are rolled.
## window_size = the size of the window of analysis
## lag_width = the number of lags to be analyzed.

## for type 2:

## delay = the delay introduced to the time series
## embed = the embedding dimensions; e.g., 1
## rescale (1 = mean, 2 = max) =  the method used to rescale the distance matrix
## radius = the radius when calculating euclidean distance, 
## NOTE: set it very small (e.g., 0.00001) when dealing  with categorical data.
##  normalize (0 = do not normalize, 1 = center around the mean, 2 = z-score transformation) the timeseries:
## minline = the minimum length of off-diagonal recurrent lines in the recurrence plot (e.g., 1)
## together with this arguments add those associated with
## the method used ("profile", "window") described above

.packageName <- 'crqa'

runcrqa <- function(ts1, ts2, par){

    datatype = thrshd = type = method = ws = radius = windowstep = windowsize =
        step = embed =  delay = rescale = normalize = mindiagline =
        minvertline = lagwidth = tw = whiteline = recpt = pad = NULL
    ## stupid initialization to please CRAN
    
    for (v in 1:length(par)) assign(names(par)[v], par[[v]])
    ## assign parameters

    tryCatch({
        ## set up a tryCatch to detach parameters values if errors occur

        res = checkts(ts1, ts2, datatype, thrshd, pad)
        ## first check that sequences
        ## have the same length
 
        if ( res[[2]] == TRUE ){

            tsnorm = res[[1]]
            ts1 = tsnorm[,1]; ts2 = tsnorm[,2]
    
            switch(type,
                   ##   Ways of calculating recurrence
    
                   {1
                    ## Quick Recurrence Profile (only diagonal)
                    
                    if (method == "profile"){
                 
                        res = drpdfromts(ts1, ts2, ws, datatype, radius)
              
                    }

                    
                    if (method == "window"){
                 
                        res = windowdrp(ts1, ts2, step, windowsize,
                            lagwidth, datatype, radius)
              
                    }
              
                }, ## close type one 
       
                   {2  ## CRQA in-depth measures (maxline, determinims, etc.)
              
              
                    if (method == "profile"){
                
                        res = crqa(ts1, ts2, delay, embed,
                            rescale, radius, normalize, mindiagline,
                            minvertline, tw, whiteline, recpt)
                        
                    }

                    if (method == "window"){
                        
                        res = wincrqa(ts1, ts2, windowstep,
                            windowsize, delay, embed, rescale,
                            radius, normalize, mindiagline, minvertline,
                            tw, whiteline, trend = F)
              
                    }
               
                } ## close type 2
             
                   ) ## close the switch option
      
        }

        else { print (paste ("Sequences differ by", res[[1]], "units") )
           }
  
#        print( paste( "Finished computing!") )
  
        return(res)
  
    }, warning = function(war) { ## here exception can be handled
    
    # warning handler picks up where error was generated
    # maybe an handler with Restarts point would be better
        print(paste("WARNING:  ", war))

    }, error = function(err) {
 
     # warning handler picks up where error was generated
        print(paste("ERROR:  ",err))

    })

}
