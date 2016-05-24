boot.plots <- function(object, comp = object$ncomp, parm = NULL, 
                       type = c("coefs", "weights", "loadings")) {
  if ((object$val.method == "none" | object$val.method == "loo")) {
    stop("No bootstrapping was done for this model")
  }
  type <- match.arg(type)
  if(is.null(parm)) {
    parm <- sample(names(object$Xdata), 1)
    hist.dat <- switch(type, 
      coefs = subset(do.call("rbind", object$validation$coefficients), 
        row.names(do.call("rbind", object$validation$coefficients)) == parm)[, comp],
      weights = subset(do.call("rbind", object$validation$weights), 
        row.names(do.call("rbind", object$validation$weights)) == parm)[, comp], 
      loadings = subset(do.call("rbind", object$validation$loadings), 
        row.names(do.call("rbind", object$validation$loadings)) == parm)[, comp])
    df <- data.frame(x = hist.dat)  
  } else {
    if(length(parm) != 1) {
      stop("Only one parameter at a time")
    }
    if(is.numeric(parm)) {
      stop("We need an actual variable name")
    }
    hist.dat <- switch(type, 
      coefs = subset(do.call("rbind", object$validation$coefficients), 
        row.names(do.call("rbind", object$validation$coefficients)) == parm)[, comp],
      weights = subset(do.call("rbind", object$validation$weights), 
        row.names(do.call("rbind", object$validation$weights)) == parm)[, comp], 
      loadings = subset(do.call("rbind", object$validation$loadings), 
        row.names(do.call("rbind", object$validation$loadings)) == parm)[, comp])
    df <- data.frame(x = hist.dat)
  }
  print(with(df, ggplot(df, aes(sample = x)) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    stat_qq() + 
    ylab("Bootstrap Sample") + 
    xlab("Quantiles") + 
    ggtitle(paste("QQPlot for Bootstrap Distribution for", type, "\nparameter no.", parm, "\ncomponent no.", comp))))
  print(with(df, ggplot(df, aes(x = x)) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_histogram(aes(fill = ..count..)) + 
    xlab(type) + 
    ggtitle(paste("Bootstrap Distribution for", type, "\nparameter no.", parm, "\ncomponent no.", comp))))
}