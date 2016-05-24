env_dif_logistic <-
function(formula, cluster = 'default', data, ...){
  formula0 <- formula
  data0 <- data[,all.vars(formula0)]
    if(nrow(data0) %% 2 != 0){print("Please remove all unmatched observations and try again using paired data")}
  convert <- all.vars(formula0)[2:length(all.vars(formula0))]
    if(cluster=='default'){cluster0 <- rep(1:(nrow(data0)/2), each=2)}
    else {cluster0 <- cluster}
      newnamesBb <- rep("Bb.", length(convert))
      newnamesBw <- rep("Bw.", length(convert))  
      data0[,(1+ncol(data0)):(ncol(data0)+2*length(convert))] <- NA
    for(j in 1:length(convert)){
    newnamesBb[j] <- paste(newnamesBb[j], convert[j], sep="")
    newnamesBw[j] <- paste(newnamesBw[j], convert[j], sep="")
      for(k in seq(1, nrow(data0), by=2)){
      data0[k, j + ncol(data0)-2*length(convert)] <- mean(c(data0[k,convert[j]], data0[k+1,convert[j]]))
      data0[k+1, j + ncol(data0)-2*length(convert)] <- data0[k, j + ncol(data0)-2*length(convert)]
      data0[k, j + ncol(data0)-length(convert)] <- data0[k,convert[j]] - data0[k, j + ncol(data0)-2*length(convert)]
      data0[k+1, j + ncol(data0)-length(convert)] <- data0[k+1,convert[j]] - data0[k+1, j + ncol(data0)-2*length(convert)]
      }
    }
    colnames(data0)[(1+ncol(data0)-2*length(convert)):(ncol(data0))] <- c(newnamesBb, newnamesBw)
  formula1 <- paste( "cbind(cluster0,", colnames(data0)[1], ") ~", "0 + ", paste(colnames(data0)[(ncol(data0)-length(convert)+1):ncol(data0)], collapse=" + "), sep="")
  modelfit <- mclogit(as.formula(formula1), data = data0, ...)
}
