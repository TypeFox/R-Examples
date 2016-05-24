.ROCnosp <- function(index, boolean, mask=NULL, nthres=NULL, thres=NULL, NAval=0, progress=FALSE, 
                     ones.bool=NULL, zeros.bool=NULL) {

# calculate basic roc table
tocd0 <- roctable(index, boolean, maskval=mask, nthres=nthres, thres=thres, NAval=NAval, progress=progress, 
                 ones.bool=ones.bool, zeros.bool=zeros.bool)
tocd <- tocd0$tocd
minInd <- tocd0$minInd

tocd1 <- tocd
tocd1[nrow(tocd1), "Threshold"] <- paste("<= ", minInd)

id <- order(tocd$falseAlarms1)
AUC <- sum(tocd$Model1[-length(tocd$Model1)] * diff(tocd$falseAlarms1)) + 
  sum(diff(tocd$falseAlarms1[id])*diff(tocd$Model1[id]))/2 

# calculate uncertainty
if(!is.null(mask)) {
  mask[mask == NAval] <- NA
  index <- index*mask
}
uncertain <- uncertainty(index, tocd)

# output
output <- new("Roc", table=tocd1[, c(1,6,7)], AUC=AUC, maxAUC = AUC + uncertain/2, minAUC = AUC - uncertain/2)

return(output) 

}
