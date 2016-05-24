env_dif <-
function(formula, data, ...){
  formula0 <- formula
  data0 <- data[,all.vars(formula0)]
    if(nrow(data0) %% 2 != 0){print("Please remove all unmatched observations and try again using paired data")}
  data1 <- data0
  data1[seq(1, nrow(data0), by=2),] <- data0[seq(1, nrow(data0), by=2),] - data0[seq(2, nrow(data0), by=2),]
  data1 <- data1[seq(1, nrow(data1), by=2),]
    colnames(data1) <- paste("Bw.", colnames(data1), sep="")
  formula1 <- paste( colnames(data1)[1], "~", "0 + ", paste(colnames(data1)[2:ncol(data1)], collapse=" + "))
  modelfit <- lm(formula = as.formula(formula1), data = data1, x=T)
}
