DelLessData <- function(data, ncases = 0) {
#This function deletes cases of a missing pattern with less than or equal to ncases 
  if(length(data)==0)
 {
   cat("Warning: data is empty")
   return
 }
  if (is.matrix(data)) {
  data <- OrderMissing(data)
  }
  n <- nrow(data$data)
  p <- ncol(data$data)
  ind <- which(data$patcnt <= ncases)
  spatcntz <- c(0, data$spatcnt)
  rm <- c()
  removedcases <- c()
  if(length(ind) != 0){
    #cat("cases with insufficient number of observations were removed")
    for(i in 1:length(ind))
    {
      rm <- c(rm, seq(spatcntz[ind[i]] + 1, spatcntz[ind[i] + 1]));
    }
    y <- data$data[-1 * rm, ]
    removedcases <- data$caseorder[rm]
    patused <- data$patused[-1 * ind, ]
    patcnt <- data$patcnt[-1 * ind]
    caseorder <- data$caseorder[-1 * rm]
    spatcnt <- cumsum(patcnt)
  }else {
  patused <- data$patused
  patcnt <- data$patcnt
  spatcnt <- data$spatcnt
  caseorder <- data$caseorder
  y <- data$data
  }
  
  newdata <- list(data = y, patused = patused, patcnt = patcnt,
                  spatcnt = spatcnt, g = length(patcnt), caseorder = caseorder, 
                  removedcases = removedcases)
#  colnames(newdata)<-colnames(data)
  class(newdata) <- "orderpattern"
  newdata
}