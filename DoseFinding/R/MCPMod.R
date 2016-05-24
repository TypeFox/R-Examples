## wrapper function for MCTtest and fitMod calls
MCPMod <- function(dose, resp, data = NULL, models = NULL, S=NULL,
                   type = c("normal", "general"), 
                   addCovars = ~1, placAdj = FALSE,
                   selModel = c("AIC", "maxT", "aveAIC"),
                   alpha = 0.025, df = NULL, critV = NULL,
                   doseType = c("TD", "ED"), Delta, p,
                   pVal = TRUE, alternative = c("one.sided", "two.sided"), 
                   na.action = na.fail, mvtcontrol = mvtnorm.control(),
                   bnds, control = NULL){

  direction <- attr(models, "direction")
  ## first perform multiple contrast test
  if(!is.null(data)){
    callMCT <- list(deparse(substitute(dose)), deparse(substitute(resp)), data,
                    models, S, type, addCovars, placAdj, alpha, df,
                    critV, pVal, alternative, na.action, mvtcontrol)
    test <- do.call(MCTtest, callMCT)
  } else {
    test <- MCTtest(dose, resp, data, models, S, type, addCovars, placAdj, alpha, df,
                    critV, pVal, alternative, na.action, mvtcontrol)
  }

  
  ## now pre-select models based on contrasts
  tstat <- test$tStat
  pvals <- attr(tstat, "pVal")
  if(!is.null(pvals)){
    tstat <- tstat[pvals < alpha]
  } else {
    tstat <- tstat[tstat > test$critVal]
  }
  if(length(tstat) == 0) ## stop if no model significant
    return(list(MCTtest = test, mods = NULL, modcrit = NULL, selMod = NULL, TD = NULL))

  ## fit models and calculate model selection criteria
  addArgs <- list(off=attr(models, "off"), scal=attr(models, "scal"))
  selModel <- match.arg(selModel)
  builtIn <- c("linlog", "linear", "quadratic", "linInt", "emax",
               "exponential", "logistic", "betaMod", "sigEmax")
  nams <- gsub("[0-9]", "", names(tstat)) ## remove numbers from model-names
  namsU <- unique(nams)
  
  mods <- vector("list", length(namsU));z <- 1
  if(missing(bnds)){
    if(!is.null(data)){
      cal <- as.character(match.call())
      doseVec <- data[, cal[2]]
    } else {
      doseVec <- dose
    }
    bnds <- defBnds(max(doseVec))
  } else {
    if(!is.list(bnds))
      stop("bnds needs to be a list")
  }
  if(selModel %in% c("AIC", "aveAIC")){
    if(type[1] == "normal"){
      modcrit <- AIC
    } else {
      modcrit <- gAIC
    }
  } else {
    modcrit <- function(x)
      max(tstat[attr(x, "model") == nams])
  }
  for(i in 1:length(namsU)){
    if(!is.null(data)){
      callMod <- list(deparse(substitute(dose)), deparse(substitute(resp)), data,
                      namsU[i], S, type, addCovars, placAdj, bnds[[namsU[i]]],
                      df, NULL, na.action, control, addArgs)
      mods[[i]] <- do.call(fitMod, callMod)
    } else {
      mods[[i]] <- fitMod(dose, resp, data, namsU[i], S, type, addCovars,
                          placAdj, bnds[[namsU[i]]], df, NULL,
                          na.action, control, addArgs)
    }
  }
  crit <- sapply(mods, modcrit)
  names(crit) <- names(mods) <- namsU
  attr(crit, "crit") <- selModel

  if(selModel %in% c("maxT", "AIC")){
    if(selModel == "AIC"){
      ind <- which.min(crit)
    }
    if(selModel == "maxT"){
      nam <- names(tstat)[which.max(tstat)]
      ind <- which(gsub("[0-9]", "", nam) == names(mods))
    }
    selMod <- namsU[ind] # name of selected model
  } else {
    aic <- crit-mean(crit)
    selMod <- exp(-0.5*aic)/sum(exp(-0.5*aic)) # model weights
    names(selMod) <- namsU
  }

  ## calculate target dose estimate
  tds <- NULL
  doseType <- match.arg(doseType)
  if(doseType == "TD"){
    if(missing(Delta))
      stop("\"Delta\" needs to be specified for TD estimation")
    tds <- sapply(mods, TD, Delta=Delta, direction = direction)
    attr(tds, "addPar") <- Delta
  }
  if(doseType == "ED"){
    if(missing(p))
      stop("\"p\" needs to be specified for TD estimation")
    tds <- sapply(mods, ED, p=p)
    attr(tds, "addPar") <- p
  }

  out <- list(MCTtest = test, mods = mods, modcrit=crit,
              selMod=selMod, doseEst=tds, doseType = doseType)
  class(out) <- "MCPMod"
  out
}

predict.MCPMod <- function(object,
                           predType = c("full-model", "ls-means", "effect-curve"),
                           newdata = NULL, doseSeq = NULL, se.fit = FALSE, ...){
  lapply(object$mods, function(x) predict(x, predType, newdata, doseSeq, se.fit))
}

print.MCPMod <- function(x, digits=3, eps=1e-03, ...){
  cat("MCPMod\n")

  xx <- x$MCTtest
  cat("\nMultiple Contrast Test:\n")
  ord <- rev(order(xx$tStat))
  if (!any(is.null(attr(xx$tStat, "pVal")))) {
    pval <- format.pval(attr(xx$tStat, "pVal"), digits = digits, 
                        eps = eps)
    dfrm <- data.frame(round(xx$tStat, digits)[ord], pval[ord])
    names(dfrm) <- c("t-Stat", "adj-p")
  }
  else {
    dfrm <- data.frame(round(xx$tStat, digits)[ord])
    names(dfrm) <- c("t-Stat")
  }
  print(dfrm)
  if (!is.null(xx$critVal)) {
    twoSide <- xx$alternative == "two.sided"
    vec <- c(" one-sided)", " two-sided)")
    cat("\n", "Critical value: ", round(xx$critVal, digits), 
        sep = "")
    if (attr(xx$critVal, "Calc")) {
      cat(" (alpha = ", xx$alpha, ",", vec[twoSide + 1], 
          sep = "")
    }
    cat("\n")
  }
  cat("\n")

  cat("Estimated Dose Response Models:")
  for(i in 1:length(x$mods)){
    cat("\n")
    cat(names(x$mods)[i], "model\n")
    cofList <- coef(x$mods[[i]], sep = TRUE)
    cof <- do.call("c", cofList)
    namcof <- c(names(cofList$DRpars), names(cofList$covarPars))
    namcof <- gsub(" ", "", namcof)  # remove white spaces for GUI
    names(cof) <- gsub("doseM", "dose", namcof) # use more obvious names
    print(round(cof, digits))
  }
  if(attr(x$modcrit, "crit") != "aveAIC"){
    cat("\nSelected model (",attr(x$modcrit, "crit"),"): ", x$selMod, "\n", sep="")
  } else {
    cat("\nModel weights (AIC):\n")
    attr(x$selMod, "crit") <- NULL
    print(round(x$selMod, 4))
  }
  
  if(is.null(length(x$doseEst)))
    return()
  if(x$doseType == "TD")
    strn <- ", Delta="
  if(x$doseType == "ED")
    strn <- ", p="
  cat("\nEstimated ",x$doseType,strn,attr(x$doseEst, "addPar"),"\n", sep="")
  attr(x$doseEst, "addPar") <- NULL
  print(round(x$doseEst, 4))
}

summary.MCPMod <- function(object, ...){
  class(object) <- "summary.MCPMod"
  print(object, digits = 3)
}

print.summary.MCPMod <- function(x, ...){
  cat("MCPMod\n\n")

  cat(rep("*", 39), "\n", sep="")
  cat("MCP part \n")
  cat(rep("*", 39), "\n", sep="")
  print(x$MCTtest)
  cat("\n")

  if(length(x$mods) == 0)
    return()
  cat(rep("*", 39), "\n", sep="")
  cat("Mod part \n")
  cat(rep("*", 39), "\n", sep="")
  for(i in 1:length(x$mods)){
    if(i > 1)
      cat("\n")
    if(length(x$mods) > 1)
      cat("** Fitted model", i,"\n")
    summary(x$mods[[i]])
  }
  cat("\n")
  cat(rep("*", 39), "\n", sep="")
  cat("Model selection criteria (",attr(x$modcrit, "crit"),"):\n", sep="")
  cat(rep("*", 39), "\n", sep="")
  crit <- attr(x$modcrit, "crit")
  attr(x$modcrit, "crit") <- NULL
  print(x$modcrit)
  if(crit != "aveAIC"){
    cat("\nSelected model:", x$selMod, "\n")
  } else {
    cat("\nModel weights (AIC):\n")
    attr(x$selMod, "crit") <- NULL
    print(round(x$selMod, 4))
  }
  
  if(is.null(length(x$doseEst)))
    return()
  cat("\n")
  cat(rep("*", 39), "\n", sep="")
  if(x$doseType == "TD")
    strn <- ", Delta="
  if(x$doseType == "ED")
    strn <- ", p="
  cat("Estimated ",x$doseType,strn,attr(x$doseEst, "addPar"),"\n", sep="")
  cat(rep("*", 39), "\n", sep="")
  attr(x$doseEst, "addPar") <- NULL
  print(round(x$doseEst, 4))
}


plot.MCPMod <- function(x, CI = FALSE, level = 0.95,
                        plotData = c("means", "meansCI", "raw", "none"),
                        plotGrid = TRUE, colMn = 1, colFit = 1, ...){
  if(is.null(x$mods))
    stop("No models significant, nothing to plot")
  plotFunc(x, CI, level, plotData, plotGrid, colMn, colFit, ...)
}

