deseasonalized_trend <- function(df, w=NULL) {
  dow = strftime(as.POSIXct(df$timestamp, origin="1970-01-01"), format="%a")
  pdf = cbind(df, dow=as.factor(dow))
  if (!is.null(w) && sum(w!=1) > 0) {
    assumed_family="binomial"
  } else {
    w = rep(1, nrow(df))
    # if all are positive integers, assume poisson; otherwise, gaussian
    if (min(df$value) >= 0 && sum(round(df$value) != df$value) == 0) {
      assumed_family = "poisson"
    } else {
      assumed_family = "gaussian"
    }
  }

  my_gam = gam(value ~ dow + s(timestamp), family=assumed_family, data=pdf, weights=w)

  pval = summary(my_gam)$anova["s(timestamp)", "Pr(F)"]
  if (length(pval) == 0) {
    pval = summary(my_gam)$anova["s(timestamp)", "P(Chi)"]
  }
  if (length(pval) == 0) {
    pval = NA
  }

  smoothed_prediction = predict(my_gam, newdata=data.frame(timestamp=pdf$timestamp, dow=as.factor(rep("Wed",nrow(pdf)))), type="response")

  # should probably be looking for a p-val of 0.01 or lower
  return(list(pval=pval, smoothed_prediction=as.vector(smoothed_prediction)))
}
