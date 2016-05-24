summary.copula = function(object, ...) {
  summaryIntamap(object, ...)
}

summary.idw = function(object, ...) {
  summaryIntamap(object, ...)
}
summary.automap = function(object, ...) {
  summaryIntamap(object, ...)
}
summary.linearVariogram = function(object, ...) {
  summaryIntamap(object, ...)
}

summary.yamamoto = function(object, ...) {
  summaryIntamap(object, ...)
}



summary.transGaussian = function(object, ...) {
  summaryIntamap(object, ...)
}

summaryIntamap = function(object, ...) {
  object2 = object
  class(object2) = "list"
  cat(paste("The object contains the following elements: \n"))
  print(summary(object2))
  if ("observations" %in% names(object)) {
    cat(paste("\n Observations:\n"))
    print(summary(object$observations, ...))
  }
  if ("variogramModel" %in% names(object)) {
    cat(paste("\n variogramModel:\n"))
    print(object$variogramModel)
  }

  if ("copulaParams" %in% names(object)) {
    cat(paste("\n copulaParams:\n"))
    print((object$copulaParams))
  }
  if ("anisPar" %in% names(object)) {
    cat(paste("\n anisPar: \n"))
    print(object$anisPar)
  }
  if ("predictions" %in% names(object)) {
    cat(paste("\n Predictions:\n"))
    print(summary(object$predictions, ...))
  }

  if ("processDescription" %in% names(object)) {
    cat(paste("\n processDescription:\n"))
    print(object$processDescription)
  }
}




plot.copula = function(x, ...) {
  plotIntamap(x, ...)
}

plot.idw = function(x, ...) {
  plotIntamap(x, ...)
}
plot.automap = function(x, ...) {
  plotIntamap(x, ...)
}
plot.linearVariogram = function(x, ...) {
  plotIntamap(x, ...)
}

plot.yamamoto = function(x, ...) {
  plotIntamap(x, ...)
}


plot.transGaussian = function(x, ...) {
  plotIntamap(x,  ...)
}



plotIntamap = function(object,zcol = "all", sp.layout = NULL, plotMat = c(2,2), ...) {
  shift = 0.03
  plots = list()
  pl = 0
  if ("variogramModel" %in% names(object) && "sampleVariogram" %in% names(object)) {
    pl = pl+1
    x = list(exp_var = object$sampleVariogram, var_model = object$variogramModel)
    class(x) = "autofitVariogram"
    plots[[pl]] = plot(x)
  }
  if (zcol == "all") zcol = expandZcol(object$outputWhat)
  if (all(zcol %in% c(names(object$predictions),"mean","variance"))) {
    predictions = object$predictions
  } else if ("outputTable" %in% names(object)){
    predictions = object$outputTable  
    transp = attr(predictions, "transposed")
    if (!is.null(transp) && transp) predictions = t(predictions)
    predictions = as.data.frame(predictions)
    coordinates(predictions) = ~x+y
  } else {
    predictions = NULL
  }
  if (!is.null(predictions)) {
    try(gridded(predictions) <- TRUE)
    if (!(all(zcol %in% c(names(predictions),"mean","variance")))) {
       zno = zcol[!(zcol %in% c(names(predictions),"mean","variance"))]
       stop(paste("plotIntamap: Column",zno, "does not exist in predictions \n"))
    }
    if (!"mean" %in% names(predictions)) predictions$mean = predictions$var1.pred
    if (!"variance" %in% names(predictions)) predictions$variance = predictions$var1.var
    for (i in 1:length(zcol)) {
      pl = pl+1
      plots[[pl]] = automapPlot(predictions, zcol = zcol[i], main = zcol[i],
        sp.layout = sp.layout, ...)
    }
  }
  nplots = length(plots)
  if (nplots == 0) {
    return()
  } else if (nplots <2) {
    print(plots[[i]]) 
  } else {  
    xinc = 1/plotMat[1]
    yinc = 1/plotMat[2]

    nmat = plotMat[1]*plotMat[2]
    ii = 0
    ij = 1
    if (!par()$ask) {
      par(ask=TRUE)
      achange = TRUE
    } else achange = FALSE
    for (i in 1:nplots) {
      ii = ii + 1
      if (ii > plotMat[1]) {
        ii = 1
        ij = ij + 1
        if (ij > plotMat[2]) ij = 1
      }
      xp = c((ii-1)*xinc,ii*xinc)
      yp = c((plotMat[2]-ij)*yinc,(plotMat[2]-ij+1)*yinc)
      print(plots[[i]],position = c(xp[1],yp[1],xp[2],yp[2]),more = nmat-(ii-1)*ij-ii) 
    }
    if (achange) par(ask = FALSE)
  }
}




expandZcol = function(outputWhat) {
  zcol = c(1:length(outputWhat))
  if ("nsim" %in% names(outputWhat)) {
    nsim = outputWhat$nsim
    zcol = c(1:(length(outputWhat)+outputWhat$nsim-1))
  }  
  i = 0
  for (j in 1:length(outputWhat)) {
    i = i+1
    what = outputWhat[i]
    if (names(what) %in% c("mean","variance")) {
      zcol[i] = names(what)
    } else if (names(what) == "nsim") {
      zcol[i:(i+nsim-1)]
      i = i+nsim-1
    } else {
      zcol[i] = paste(names(what),what[[1]],sep="")
    }
  }
  zcol
}