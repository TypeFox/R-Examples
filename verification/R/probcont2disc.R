probcont2disc <- function(x, bins = seq(0,1,0.1) ){
  ## converts continuous prob forecasts into a range of discrete
  ## probability forecsast assigned the value at the midpoint of their
  ## bin.
 
  if(prod(x >= 0 & x <= 1) != 1)
   {stop("Are you sure x is a probability? \n
     Values must be between 0 and 1 \n")}
if(max(x) >  max(bins)  | min(x) < min(bins) ){stop("
  Bins must span the interval of predictions.") }

  mids <- bins[-length(bins)] +  0.5* diff(bins)
xx<- cut(x, breaks = bins, include.lowest = TRUE)
# new <- mids[xx]  

  return(list(new = mids[xx], mids = mids))
} ## close function
