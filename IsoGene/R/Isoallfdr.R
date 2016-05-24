# Arguments :  qqstat, the output from Isoqqstat
#              ddelta contains various values for delta's
# Return FDR for various values of delta

# input :
#    qqstat: the output from function Isoqqstat
#    ddelta: can be missing
#    stat: "E2", "Williams", "Marcus", "M","ModifM"
# output :
#    delta table:
#    col 1: delta values
#    col 2: FP 50%
#    col 3: FP 90%
#    col 4: Number of significant genes
#    col 5: FDR 50%
#    col 6: FDR 90%

Isoallfdr <- function(qqstat, ddelta, stat) {

  switch(stat,
      E2 = {
        qstat <- qqstat[[1]]
        dperm <- qqstat[[2]]},
      Williams = {
        qstat <- qqstat[[3]]
        dperm <- qqstat[[4]]},
      Marcus = {
        qstat <- qqstat[[5]]
        dperm <- qqstat[[6]]},
      M = {
        qstat <- qqstat[[7]]
        dperm <- qqstat[[8]]},
      ModifM = {
        qstat <- qqstat[[9]]
        dperm <- qqstat[[10]]  
      })

  k1 <- nrow(qstat)
   
  # Use percentiles if ddelta is not specified
  if (missing(ddelta)) {
      qd <- round(quantile(abs(qstat[,3]), c(0.01, 0.999)), 2)
      bb <- (qd[2] - qd[1]) / 100
      b0 <- 10^-floor(log(bb, 10))
      bb <- trunc(bb * b0) / b0
      ddelta <- seq(from = qd[1], to = qd[2], by = bb)
  }

  k2 <- length(ddelta)
  # Suppress warning messages if uupcut is empty
  # Determine the number of significant genes for each value in ddelta
  low.point <- up.point <- clow <- cup <- sn <- NULL
  nsn <- array(0, c(length(ddelta), 2))
  for (i in 1:k2){
      if (sum(qstat[qstat[,1] < 0, 3] < - ddelta[i])>0){
        low <- max(qstat[qstat[,1] < 0, 1][qstat[qstat[,1] < 0, 3] < - ddelta[i]])
      } else {
        low <- -Inf
      }
      if (low != "-Inf") { # use  !is.finite(low) combined with negative ??
         which.low <- which(qstat[qstat[,1] < 0, 1] == low)
         low.point[i] <- which.low[length(which.low)]
         # max(sort(as.numeric(grep("TRUE",as.character(test.low)))))
      }
      if (low == "-Inf"){
         low <- min(qstat[,1])
         low.point[i] <- 0
      }
      clow[i] <- low

      if (sum(qstat[qstat[,1] > 0,3] > ddelta[i])>0){ 
        up <- min(qstat[qstat[,1] > 0, 1][qstat[qstat[,1] > 0,3] > ddelta[i] ]) 
      } else {
        up <- Inf
      }

      if (up != "Inf") {
        which.up <- which(qstat[qstat[,1]>0,1] == up)
        up.point[i] <- which.up[1]
      }
      if (up == "Inf") {
         up <- max(qstat[qstat[,1] > 0,1])
         up.point[i] <- sum(qstat[,1] > 0) + 1
      }
      cup[i] <- up
      sn[i] <- sum(qstat[,1] > 0) + 1 - sum(up.point[i]) + sum(low.point[i])
      qq <- c(rep(1, nrow(dperm)) %*% ((dperm > cup[i]) | (dperm < clow[i])))
      nsn[i, ] <- quantile(qq, c(0.5, 0.9))
             
  }
  prun <- quantile(dperm, c(0.25, 0.75))
  p <- min(sum(qstat[,1] < prun[2] & qstat[,1] > prun[1]) / (nrow(qstat) *0.5), 1)
  nsn <- nsn * p
  dtable = round(cbind(Ddelta = ddelta, "FalsePositive50%" = nsn[,1],
      "FalsePositive90%" = nsn[,2], Called = sn, "FDR50%" = nsn[,1]/sn,
      "FDR90%" = nsn[, 2]/sn), 4)
  
  return(dtable)
}
