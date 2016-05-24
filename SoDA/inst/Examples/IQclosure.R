newIQ <- function(nData = 1000, probs = seq(0, 1, 0.25)) {
   dataBuf <- numeric(0)
   qBuf <- numeric(0)

   addData <- function(newdata) {
       n <- length(newdata);
       if(n + length(dataBuf) > nData)
         recompute(newdata)
       else
         dataBuf <<- c(dataBuf, newdata)
   }

   recompute <- function(newdata = numeric(0)) {
       qBuf <<- doQuantile(qBuf, c(dataBuf, newdata), probs)
       dataBuf <<- numeric(0)
   }

   getQ <- function() {
       if(length(dataBuf) > 0)
         doQuantile(qBuf, dataBuf, probs)
       else
         qBuf
   }
   list(addData = addData, getQ = getQ)
}
