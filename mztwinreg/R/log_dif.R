log_dif <-
function(formula, data, ...){
  formula0 <- formula
  data0 <- data[,all.vars(formula0)]
    if(nrow(data0) %% 2 != 0){print("Please remove all unmatched observations and try again using paired data")}
  convert <- all.vars(formula0)[1]
    data0[,ncol(data0)+1] <- NA
    data0[seq(1, nrow(data0), by=2), ncol(data0)] <- log(data0[seq(1, nrow(data0), by=2),convert] / data0[seq(2, nrow(data0), by=2),convert])
    data0 <- data0[seq(1, nrow(data0), by=2), c(all.vars(formula0)[2:length(all.vars(formula0))], colnames(data0)[ncol(data0)])]
    colnames(data0)[ncol(data0)] <- as.character(convert)
  modelfit <- ols(as.formula(formula0), data = data0, ...)
}
