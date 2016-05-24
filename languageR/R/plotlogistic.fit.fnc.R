`plotlogistic.fit.fnc` <-
function(x, data, method = "cut", where = seq(0, 1, by=0.1), scalesize=NA, ...) {
  require(lme4, quietly = TRUE)
  require(rms, quietly = TRUE)
  if (class(x)[1]=="mer") {
    y=attr(x@frame, "terms")
    depvar = names(attr(terms(y),"dataClasses")[attr(terms(y),"response")])
    probs = fitted(x)   # this used to be the fitted logit, but is now the fitted probability
  } else {
    if (class(x)[1] == "lrm") {
      depvar = as.character(formula(x$call))[2]
      probs = predict(x,type="fitted")
    } else {
      stop("first argument is not an lmer or lrm model")
    }
  }
  if (method == "cut") {
    classes = cut2(probs, where, levels.mean = TRUE)
    classCounts = table(classes)
    means = tapply(as.numeric(data[,depvar])-1, classes, mean)
    means = means[!is.na(means)]
  } else {
    if (method == "shingle") {
      sh = equal.count(probs)
      means = rep(0, length(levels(sh)))
      midpoints = rep(0, length(means))
      for (i in 1:length(levels(sh))) {
        means[i] = mean(probs[probs>levels(sh)[[i]][1] & probs < levels(sh)[[i]][2]])
        midpoints[i] = as.character(mean(levels(sh)[[i]]))
      }
      names(means) = as.character(midpoints)
    }
  }
  plot(as.numeric(names(means)), means, 
     xlab = "mean predicted probabilities", 
     ylab = "observed proportions", type="n")
  abline(0, 1, col = "grey")
  if ((method=="cut") & (!is.na(scalesize))) {
    symbols(as.numeric(names(means)),as.numeric(means), circles=as.numeric(classCounts), 
     inches=scalesize, main=" ", add=T)
  } else {
    points(as.numeric(names(means)), means)
  }

  mtext(paste("R-squared: ", round(cor(as.numeric(names(means)), means)^2,2),
  sep=""), 3, 1.5)
}

