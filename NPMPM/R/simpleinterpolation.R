simpleinterpolation <-
function(psteptime, series, cfu_in){

# GET THE MEASURED VALUES FOR ONE SERIES OF MEASURED VALUES
mvalue <- measured_values[(measured_values$raw_data_key==as.vector(series)),]

cfu <- log10(cfu_in)
tprevious <- max(mvalue[mvalue$time<=psteptime,]$time)
cfuprevious <- mvalue[mvalue$time==tprevious,]$logc
tlater<- min(mvalue[mvalue$time>=psteptime,]$time)
cfulater <- mvalue[mvalue$time==tlater,]$logc

if (tprevious==tlater) {
cfu_out <-  cfuprevious

} else { 
bla <- (cfulater - cfuprevious) / (tlater - tprevious)
cfu_out <- cfuprevious + bla* (psteptime-tprevious)
} # end else

return(cfu_out)

} # END FUNCTION

