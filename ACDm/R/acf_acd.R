acf_acd <- function(fitModel = NULL, conf_level = 0.95, max = 50, min = 1){
  lag <- acf <- NULL
  
  if("acdFit" %in%  class(fitModel)){
    dur <- fitModel$durations$durations
    adjDur <- fitModel$durations$adjDur
    resi <- fitModel$residuals
  } else if("data.frame" %in%  class(fitModel)){
    dur <- fitModel$durations
    adjDur <- fitModel$adjDur
    resi <- fitModel$residuals
  } else if(is.vector(fitModel)){
    dur <- fitModel
    adjDur <- NULL
    resi <- NULL
  } else stop("fitModel is not of the correct object type")
  
  df <- data.frame()
  if(length(dur) != 0){
    temp_acf <- stats::acf(dur, plot = FALSE, lag.max = max)
    df <- rbind(df, data.frame(acf = temp_acf$acf[-(1:min)], lag = temp_acf$lag[-(1:min)], data = "durations"))
    conf <- stats::qnorm((1 - conf_level)/2)/sqrt(length(dur))
  }
  if(length(adjDur) != 0){
    temp_acf <- stats::acf(adjDur, plot = FALSE, lag.max = max)
    df <- rbind(df, data.frame(acf = temp_acf$acf[-(1:min)], lag = temp_acf$lag[-(1:min)], data = "adj. durations"))
    conf <- stats::qnorm((1 - conf_level)/2)/sqrt(length(adjDur))
  }
  if(length(resi) != 0){
    temp_acf <- stats::acf(resi, plot = FALSE, lag.max = max)
    df <- rbind(df, data.frame(acf = temp_acf$acf[-(1:min)], lag = temp_acf$lag[-(1:min)], data = "residuals"))
    conf <- stats::qnorm((1 - conf_level)/2)/sqrt(length(resi))
  }
  
  g <- ggplot2::ggplot(df, ggplot2::aes(x = lag, y = acf))
  g <- g + ggplot2::geom_bar(stat = "identity", position = "identity") + ggplot2::ylab("autocorrelation") + ggplot2::xlab("lag")
  g <- g + ggplot2::geom_hline(yintercept = -conf, color = "blue",size = 0.2) + ggplot2::geom_hline(yintercept = conf, color = "blue",size = 0.2)
  g <- g + ggplot2::geom_hline(yintercept = 0, color = "red", size = 0.3) + ggplot2::theme_bw(base_size=20) + ggplot2::facet_wrap(~data) 
  
  print(g)
  acf_acd <- df
}