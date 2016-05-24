
## written by Moreno I. Coco 2013 (moreno.cocoi@gmail.com)
## original Matlab code by Rick Dale

## calculate recurrence over different sized windows
## arguments: step = interval considered on the serie;
##            windowstep = the step of the window. 
##            windowsize =  the size of the window
##           
##          
## other arguments to pass are the same as crqa

# windowsize = 200; windowstep = 50;
# type = 1; delay = 1; embed = 1;
# rescale = 1; radius = 0.0001;
# normalize = 0; minline = 2

# source("crqa.R")

# tS = simts(0.25, 0.05, 0.2, 0.2, 0.25, 1000)
# ts1 = tS[1,]; ts2 = tS[2,]

.packageName <- 'crqa'

wincrqa <- function(ts1, ts2, windowstep, windowsize, delay, embed,
                    rescale, radius, normalize, mindiagline,
                    minvertline, tw = 0, whiteline = F, trend = FALSE){

    ## we do not expect as input a recurrent plot
    ts1 = as.vector(as.matrix(ts1));   ts2 = as.vector(as.matrix(ts2))
    points = seq(1, (length(ts1) - (windowsize)-1), windowstep)
    
    TREND = crawin = vector()
    tsp = 0 ## set a counter with the win at which rec was found    

    i = 1
    for (i in points){
        tsp = tsp +1
        
        ts1win = ts1[i:(i+windowsize)];
        ts2win = ts2[i:(i+windowsize)];
        
        ans = crqa(ts1win, ts2win, delay, embed, rescale,
            radius, normalize, mindiagline, minvertline, tw,
            whiteline)
        
        RP = ans$RP
        if (length(RP) == 1) RP = vector() ## a trick for cases
                                           ## with empty recurrence plot
        ans = as.numeric( unlist(ans[1:9]) )
        ## if trend needs to be calculated do it here

        if (trend == TRUE){

            if (length(RP) > 0){
            
                NX = ncol(RP)
                T = vector("numeric", length = NX-1)

                k = 1
                for (k in 1:(NX-1)){
                    T[k] = length(
                         which(diag(RP, k) != F)) / (NX-k)*100;
                }
            
                Ntau = NX - 1 - round(0.1*NX);

                ## last 10% of the RP will be skipped
                p = polyfit(2:(Ntau+1),T[1:Ntau], 1) # slope
                TREND = c(TREND, 1000 * p[1])
            ## Webber's definition includes factor 1000

            } else { TREND = NA}
        } else {
            TREND = c(TREND, NA)
        }
        
        crawin = rbind(crawin, c(ans, tsp), deparse.level = 0)
        
    }
 
    return(list(crqwin = crawin, TREND = TREND))
    
}

