cfuinterpolation <-
function(psteptime, series, cfu_in){

# GET THE MEASURED VALUES FOR ONE SERIES OF MEASURED VALUES
mvalue <- measured_values[(measured_values$raw_data_key==as.vector(series)),]
bla <- lm(mvalue$logc~mvalue$time)
n0 <- mvalue[mvalue$time==0,]$logc
n0 <- log10(cfu_in)-(n0[1]-coef(bla)[1])
cfu_out <- n0 + (coef(bla)[2])*psteptime
return(cfu_out)
} # END FUNCTION

