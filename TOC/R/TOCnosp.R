.TOCnosp <- function(index, boolean, mask=NULL, nthres=NULL, thres=NULL, NAval=0, P=NA, Q=NA, progress=FALSE,  
                     ones.bool=NULL, zeros.bool=NULL, population=NULL, units=character(0)) {

  boolval <- boolean
  # calculate population if not given
  if (is.null(population)){
    if(!is.null(mask)){
      mask[mask == NAval] <- NA
      boolval <- boolean*mask
    }
    # extract total number of no NA cells
    boolvals <- boolval[!is.na(boolval)]
    population <- length(boolvals)
    if(!is.na(P) & !is.na(Q)){
      population <- P + Q 
    }
  }
  
# calculate basic roc (toc) table
tocd0 <- roctable(index, boolean, maskval=mask, nthres=nthres, thres=thres, NAval=NAval, progress=progress, 
                 ones.bool=ones.bool, zeros.bool=zeros.bool)
tocd <- tocd0$tocd
minInd <- tocd0$minInd

# calculate additional columns for toc
maxA <- tocd$A[nrow(tocd)]
maxB <- tocd$B[nrow(tocd)]
tocd$m <- (maxA - tocd$A)/(maxA + maxB) 
tocd$h <- tocd$A/(maxA + maxB)
tocd$f <- tocd$B/(maxA + maxB) 
tocd$c <- 1 - tocd$m - tocd$h -tocd$f
prevalence <- tocd$h[nrow(tocd)]

tocd$Hits <- tocd$Model1 * prevalence * population
tocd$hitsFalseAlarms <- tocd$Hits + tocd$falseAlarms1*(1-prevalence)*population
tocd$hitsMisses <- prevalence*population
tocd$maximum <- pmin(tocd$hitsMisses, tocd$hitsFalseAlarms)
tocd$minimum <- pmax(0, tocd$hitsFalseAlarms + tocd$hitsMisses - population)
tocd$Uniform1 <- tocd$hitsMisses * tocd$hitsFalseAlarms / population

tocd1 <- tocd
tocd1[nrow(tocd1), "Threshold"] <- paste("<= ", minInd)

tocd2 <- tocd1[, c("Threshold", "hitsFalseAlarms", "Hits")]

# adjustment to population if user provides P and Q
if(!is.na(P) & !is.na(Q)){
tocd1$hitsFalseAlarmsP <- P * tocd1$Model1 + Q * tocd1$falseAlarms1
tocd1$HitsP <- P * tocd1$Model1
tocd2 <- tocd1[, c("Threshold", "hitsFalseAlarms", "Hits", "hitsFalseAlarmsP", "HitsP")]
}

# calculate totalAUC in data units and AUC as a proportion 
id <- order(tocd2$hitsFalseAlarms)
totalAUC <- sum(tocd2$Hits[-length(tocd2$Hits)] * diff(tocd2$hitsFalseAlarms)) + 
  sum(diff(tocd2$hitsFalseAlarms[id])*diff(tocd2$Hits[id]))/2 - ((prevalence * population)^2)/2 

AUC <- totalAUC/(population * prevalence * population - (prevalence * population)^2)

colnames(tocd2)[2] <- "Hits+FalseAlarms"
if (any(colnames(tocd2) == "hitsFalseAlarmsP")) colnames(tocd2)[4] <- "Hits+FalseAlarmsP"

# calculate uncertainty
if(!is.null(mask)) index <- index*mask
uncertain <- uncertainty(index, tocd)

# output
output <- new("Toc", table=tocd2, prevalence=prevalence*population, population=population, units=units, AUC=AUC,  
               maxAUC = AUC + uncertain/2, minAUC = AUC - uncertain/2)
return(output) 

}
