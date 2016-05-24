xBalance.make.stratwts <- function(stratum.weights,ss.df,goodstrat.df,zz,data,normalize.weights) {
  if (is.function(stratum.weights)) {
    swt.ls <- rep(list(stratum.weights), length(ss.df))
    names(swt.ls) <- names(ss.df)
  }
  if (is.list(stratum.weights) & !all(names(ss.df)%in%names(stratum.weights)))
    stop("list stratum.weights must have entry names matching those of stratifying factors")

  if (!is.list(stratum.weights) & !is.function(stratum.weights) & length(ss.df)>1)
    stop("stratum weights must be specified for each stratifying factor")

  if (!is.list(stratum.weights) & !is.function(stratum.weights)) {
    swt.ls <- list(stratum.weights)
    names(swt.ls) <- names(ss.df)
  }

  if (is.list(stratum.weights) & !is.function(stratum.weights)) {
    swt.ls <- stratum.weights
    names(swt.ls) <- names(ss.df)
  }




  ### change names here!

  wtlist <- list()
  for (nn in names(swt.ls)) {

    if (is.function(swt.ls[[nn]])) {
      sweights <-
        do.call(swt.ls[[nn]],
                args=list(data=data.frame(Tx.grp=zz[goodstrat.df[[nn]]],
                              stratum.code=factor(ss.df[goodstrat.df[[nn]],nn]),
                              data[goodstrat.df[[nn]],,drop=FALSE])),
                envir=parent.frame())
    } else {
      if (!is.numeric(swt.ls[[nn]]))
        stop("stratum.weights must be an expression or numeric vector")

      if (is.null(names(swt.ls[[nn]])))
        stop ("if stratum.weights is a vector, must have names")

      if (!(all(levels(factor(ss.df[[nn]])) %in% names(swt.ls[[nn]])) ))
        stop("if stratum.weights is a vector, must have a name for each stratum")

      sweights <- swt.ls[[nn]][levels(factor(ss.df[[nn]]))]
    }

    if (all(is.na(sweights)))
      stop(paste("All stratum weights NA (strat.",nn,")."))
    if (any(is.na(sweights))) {
      sweights[is.na(sweights)] <- 0
      warning(paste("NAs in stratum weights (",nn," strat.); to be interpreted as 0s.", sep=""))
    }
    if (any(sweights<0))
      stop("stratum weights must be nonnegative")

    if (normalize.weights)
      sweights <- sweights/sum(sweights, na.rm=TRUE)

    if (identical(harmonic, swt.ls[[nn]])) {
      hwts <- sweights
    } else {
      hwts <- harmonic(data.frame(Tx.grp=zz[goodstrat.df[[nn]]],
                                  stratum.code=factor(ss.df[goodstrat.df[[nn]],nn]),
                                  data[goodstrat.df[[nn]],,drop=FALSE]))
    }
    hwts <- hwts/sum(hwts, na.rm=TRUE)

    wtratio <- unsplit(sweights/hwts, ss.df[[nn]], drop=TRUE)[goodstrat.df[[nn]]]
    wtlist[[nn]] <- list(sweights=sweights,wtratio=wtratio)
    NULL
  }
  wtlist
}
