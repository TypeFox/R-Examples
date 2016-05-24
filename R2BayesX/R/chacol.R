chacol <- function(x)
{
  if(!is.null(cn <- colnames(x))) {
    for(i in 1L:length(cn)) {
      if(cn[i] == "deviance")
        cn[i] <- "Deviance"
      if(cn[i] == "pmean")
        cn[i] <- "Mean"
      if(cn[i] == "pmed")
        cn[i] <- "50%"
      if(cn[i] == "pstddev")
        cn[i] <- "Sd"
      if(cn[i] == "pmin")
        cn[i] <- "Min"
      if(cn[i] == "pmax")
        cn[i] <- "Max"
      if(cn[i] == "pmode")
        cn[i] <- "Estimate"
      if(cn[i] == "ci95lower")
        cn[i] <- "2.5%"
      if(cn[i] == "ci80lower")
        cn[i] <- "10%"
      if(cn[i] == "std")
        cn[i] <- "std. error"
      if(cn[i] == "ci80upper")
        cn[i] <- "90%"
      if(cn[i] == "ci95upper")
        cn[i] <- "97.5%"
      if(cn[i] == "variance")
        cn[i] <- "Variance"
      if(cn[i] == "smoothpar")
        cn[i] <- "Smooth Par."
      if(cn[i] == "stopped")
        cn[i] <- "Stopped"
      if(cn[i] == "loglike")
        cn[i] <- "logLik"
      if(cn[i] == "aic")
        cn[i] <- "AIC"
      if(cn[i] == "bic")
        cn[i] <- "BIC"
      if(cn[i] == "gcv")
        cn[i] <- "GCV"
      if(cn[i] == "dic")
        cn[i] <- "DIC"
      if(grepl("pqu", cn[i], fixed = TRUE)) {
        cn[i] <- paste(gsub("pqu", "", cn[i], fixed = TRUE), "%", sep = "")
        if(grepl("p", cn[i], fixed = TRUE))
          cn[i] <- gsub("p", ".", cn[i], fixed = TRUE)
      }
    }
    colnames(x) <- cn
  }

  return(x)
}

