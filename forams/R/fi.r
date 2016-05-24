`fi` <-
  function (df, groups) {
    Ps <- NA
    Po <- NA
    Ph <- NA
    for(i in 1:ncol(df)) {
      Ps[i] <- as.matrix(sum(df[i][groups=='Ps']) / sum(df[i]))
    }                          
    for(i in 1:ncol(df)) {
      Po[i] <- as.matrix(sum(df[i][groups=='Po']) / sum(df[i]))
    }
    for(i in 1:ncol(df)) {
      Ph[i] <- as.matrix(sum(df[i][groups=='Ph']) / sum(df[i]))
    }
    FI <- new("fi")
    FI@fi <- as.data.frame((10 * Ps) + Po + (2 * Ph))
    FI@fi <- cbind(as.data.frame(1:ncol(df)), FI@fi)
    colnames(FI@fi) <- c('PlotOrder', 'FI')
    rownames(FI@fi) <- colnames(df)
    FI@fi <- as.data.frame(FI@fi)
    class(FI) <- "fi"
    return(FI)
  }

