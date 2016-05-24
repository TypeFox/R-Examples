jk.after.boot <- function(object, ncomp = object$ncomp, 
                          type = c("coefficients", "loadings", "weights"),  
                          parm = NULL) {
  if ((object$val.method == "none" | object$val.method == "loo")) {
    stop("No bootstrapping was done for this model")
  }
  if ((length(ncomp) != 1)) {
    stop("Only one component at a time for now")
  }
  if(is.null(parm)) {
    parm <- sample(coefficients(object)[[1]]$variable, 1)
  } else {
    parm <- parm
  }
  if ((length(parm) != 1)) {
    stop("Only one variable at a time for now")
  }
Obj <- nrow(object$Xdata)
jab <- llply(1:Obj, function(this.obs) {
  boot.length <- 1:length(object$validation$in.bag)
  A <- llply(boot.length, function(x) {
    this.obs %in% object$validation$in.bag[[x]]
  })
  B <- do.call("rbind", llply(boot.length, function(x) !(FALSE %in% A[[x]])))
  B <- data.frame(In = B, Row = 1:nrow(B))
  IN <- B[B$In == TRUE, 2]
  OUT <- B[B$In == FALSE, 2]
  IN.Samp <- do.call("rbind", object$validation[names(object$validation) == type][[1]])
  IN <- subset(IN.Samp, row.names(IN.Samp) %in% parm)[, ncomp]
  OUT <- IN[OUT]
  coefficients.boot.in.dat <- IN
  coefficients.boot.out <- OUT
  coefficients.boot.out.cis <- mean(coefficients.boot.out, na.rm = T)
  coefficients.boot.in.cis <- mean(coefficients.boot.in.dat, na.rm = T)
  inf <- coefficients.boot.out.cis
  Empirical.Quantiles.in <- quantile(coefficients.boot.in.dat, 
                                     prob = c(0.05, 0.1, 0.16, 0.5, 0.84, 0.9, 0.95), na.rm = T)
  Empirical.Quantiles.out <- quantile(coefficients.boot.out, 
                                      prob = c(0.05, 0.1, 0.16, 0.5, 0.84, 0.9, 0.95), na.rm = T)
  Out <- data.frame(Empirical.Quantiles.in, Empirical.Quantiles.out)
  Out$Percentile <- row.names(Out)
  Out$Obs <- this.obs
  Out$inf <- inf
  Out
})
  df <- ldply(jab)
  infl.bt <- (Obj - 1) * (mean(unique(df$inf), na.rm = T) - unique(df$inf))
  infl.values <- (round(infl.bt / sqrt(sum(infl.bt^2 / (Obj - 1))), 2))
  Infl.values <- data.frame(infl.values, Obs = 1:Obj)
  df <- merge(df, Infl.values, by = "Obs")
  # df.m <- melt(df, measure = c(2, 3))
print(with(df, ggplot(df, aes(infl.values, Empirical.Quantiles.out)) + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        geom_text(aes(label = as.factor(Obs), col = as.factor(Obs))) + 
        theme(legend.position = "none") + 
        geom_line(aes(group = Percentile)) + 
        geom_hline(data = df, aes(yintercept = Empirical.Quantiles.in), lty = 2) + 
        xlab("Estimate of Relative Influence") + 
        ylab("Empirical Quantiles") + 
        ggtitle(paste("Jackknife-After-Boot (",  type,")\nfor Variable =", parm, "\nncomp =", ncomp)) +     
        theme(plot.title = element_text(size = 20)) + 
        theme(axis.title.x = element_text(size = 20)) +
        theme(axis.title.y = element_text(size = 20, angle = 90)) + 
        theme(axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, face = "bold")) + 
        theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold"))))
}


