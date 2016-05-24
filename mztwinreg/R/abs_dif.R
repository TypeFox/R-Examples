abs_dif <-
function(formula, regression = 'linear', data, ...){
  formula0 <- formula
  data0 <- data[,all.vars(formula0)]
    if(nrow(data0) %% 2 != 0){print("Please remove all unmatched observations and try again using paired data")}
  convert <- all.vars(formula0)[1]
    data0[,ncol(data0)+1] <- NA
    colnames(data0)[ncol(data0)] <- paste("absDiff.", convert, sep="")
    data0[seq(1, nrow(data0), by=2), ncol(data0)] <- abs(data0[seq(1, nrow(data0), by=2),convert] - data0[seq(2, nrow(data0), by=2),convert])
    data0 <- data0[seq(1, nrow(data0), by=2), c(all.vars(formula0)[2:length(all.vars(formula0))], colnames(data0)[ncol(data0)])]
  formula1 <- paste( paste("absDiff.", convert, sep=""), "~", paste(all.vars(formula0)[2:length(all.vars(formula0))], collapse=" + ") )
  if(regression=='linear'){modelfit <- ols(as.formula(formula1), data = data0, ...)
  } else if(regression=='logistic'){modelfit <- lrm(as.formula(formula1), data = data0, ...)}
}
