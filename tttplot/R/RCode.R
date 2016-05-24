### Time-to-Target Plot - version in R

tttPlot <- function(timeValue = NULL, tGraph = "TTTPlot", snTheorical = FALSE) {
      if(is.null(timeValue)) return(list("xSortVal" = c(-1), "probTV" = c(-1)));
      time_value = timeValue;
      nTV        = length(time_value);
      
    # Sort values
      sorted_time_value = sort(time_value);

    # Calcule the AVG
      avgTV      = sum(time_value) / nTV;

    # Compute probabilities for distribution plot
      probTV     = c(1:nTV) * 0;
      np1        = nTV+1;
      for(nn in 1:nTV){
        probTV[nn] = nn + 0.5;
        probTV[nn] = probTV[nn] / np1;
      }

    # Compute distribution parameters
      fq         = (np1 * 0.25);
      tq         = (np1 * 0.75);
      fq         = floor(np1 * 0.25);
      tq         = floor(np1 * 0.75);
      y          = probTV[fq];
      zl         = sorted_time_value[fq];
      ql         = -log(1-y);
      y          = probTV[tq];
      zu         = sorted_time_value[tq];
      qu         = -log(1-y);
      lambda     = (zu - zl)/(qu - ql);
      mu         = zl - (lambda * ql);
      shifted_mean = mu+lambda;

    # Compute theoretical plot (400 points)
      tmax       = sorted_time_value[nTV-1];
      inv_lambda = 1/lambda;
      eps        = tmax/400;
      theory_t   = c(1:400) * 0;
      theory_p   = theory_t;
      for(nn in 1:400){
        theory_t[nn] = eps * nn;
        theory_p[nn] = 1-exp(-inv_lambda*(eps * nn - mu));
      }

    # Compute theoretical time values
      theoretical_time = c(1:nTV) * 0;
      for(nn in 1:nTV){
        theoretical_time[nn] = -log(1-probTV[nn]);
      }

    # Compute qqplot line, lower and upper error lines.
      x              = c(1:nTV) * 0;
      qq_err         = x;
      lo_error_point = x;
      up_error_point = x;
      for(nn in 1:nTV){
        piTV               = probTV[nn];
        x[nn]              = -log(1-piTV);
        qq_err[nn]         = lambda * x[nn] + mu;
        dev                = lambda * (sqrt(piTV/((1-piTV)*np1)));
        lo_error_point[nn] = qq_err[nn] - dev;
        up_error_point[nn] = qq_err[nn] + dev;
      }

      if(tGraph == "QQPlot") {
        plot(sorted_time_value, theoretical_time, col = "black", type="l", pch = 20);
        lines(x,qq_err,col="blue");
        lines(x,up_error_point,col="red");
        lines(x,lo_error_point,col="green");
      }

      if(tGraph == "TTTPlot") {
        plot(sorted_time_value, probTV, col = "black", type="l", pch = 20);
        if(snTheorical) lines(theory_t,theory_p,col="blue");
      }

      return(list("xSortVal" = sorted_time_value, "probTV" = probTV));
}

tttPlotCompare <- function(timeValue1 = NULL, timeValue2 = NULL, tGraph = "TTTPlot", 
                           snTheorical = FALSE, xLab = "Time", yLab = "Accum. Prob.", legendTT = NULL, 
                           snReturn = TRUE, posLegend = "topleft") {
  if(is.null(timeValue1) || is.null(timeValue2)) return(list("xSortVal1" = c(-1), "probTV1" = c(-1), "xSortVal2" = c(-1), "probTV2" = c(-1)));
  tVal1 = tttPlot(timeValue1, "None", FALSE);
  tVal2 = tttPlot(timeValue2, "None", FALSE);
  nVal1 = length(tVal1$xSortVal);
  nVal2 = length(tVal2$xSortVal);
  if(is.null(legendTT)) legendTT = c("Model1", "Model2");
  matplot(cbind(tVal1$xSortVal,tVal2$xSortVal), cbind(tVal1$probTV, tVal2$probTV), pch = 1, type = "l", xlab = xLab, ylab = yLab)
  legend(posLegend, legendTT, col = c("black", "red"), lty = c(1,2), lwd = c(2,2), bty = 'n', cex = 0.5);
  
  if(snReturn) return(list("xSortVal1" = tVal1$xSortVal, "probTV1" = tVal1$probTV, "xSortVal2" = tVal2$xSortVal, "probTV2" = tVal2$probTV));
}
