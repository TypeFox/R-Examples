## here the multiple contrast test related functions

## performs multiple contrast test (MCP part of MCPMod)
MCTtest <- function(dose, resp, data = NULL, models, S = NULL,
                    type = c("normal", "general"),
                    addCovars = ~1, placAdj = FALSE, 
                    alpha = 0.025, df = NULL, critV = NULL, pVal = TRUE,
                    alternative = c("one.sided", "two.sided"),
                    na.action = na.fail, mvtcontrol = mvtnorm.control(),
                    contMat = NULL){
  ## perform multiple contrast test
  type <- match.arg(type)
  alternative <- match.arg(alternative)
  ## check for valid arguments
  cal <- as.character(match.call())
  lst <- checkAnalyArgs(dose, resp, data, S, type,
                        addCovars, placAdj, na.action, cal)
  dd <- lst$dd;type <- lst$type;S <- lst$S
  doseNam <- lst$doseNam;respNam <- lst$respNam

  ## calculate optimal contrasts and test-statistics
  doses <- unique(dd[[doseNam]])
  k <- length(doses)
  if(type == "normal"){
    dd[, doseNam] <- as.factor(dd[, doseNam])
    form <- paste(respNam, "~", doseNam, "+", addCovars[2], "-1", sep="")
    lm.fit <- lm(as.formula(form), data = dd)
    est <- coef(lm.fit)[1:k]
    vc <- vcov(lm.fit)[1:k, 1:k]
    df <- lm.fit$df.residual
  } else {
    est <- dd[[respNam]]
    vc <- S
    if(is.null(df))
      df <- Inf
  }
  if(is.null(contMat)){ # calculate optimal contrasts
    contMat <- optContr(models, doses, S=vc, placAdj=placAdj)$contMat
    rownames(contMat) <- doses
  } else { # contrast matrix specified
    if(inherits(contMat, "optContr"))
      contMat <- contMat$contMat
    if(nrow(contMat) != length(est))
      stop("contMat of incorrect dimensions")
  }
  ct <- as.vector(est %*% contMat)
  covMat <- t(contMat) %*% vc %*% contMat
  den <- sqrt(diag(covMat))
  tStat <- ct/den
  
  if(alternative == "two.sided"){
    tStat <- abs(tStat)
  }
  corMat <- cov2cor(covMat)
  
  if(is.null(critV)){
    if(!pVal){
      stop("either p-values or critical value need to be calculated.")
    }
  } else if(is.logical(critV) & critV == TRUE){
    critV <- critVal(corMat, alpha, df, alternative, mvtcontrol)  
    attr(critV, "Calc") <- TRUE # determines whether cVal was calculated
  } else { 
    pVal <- FALSE # pvals are not calculated if critV is supplied
    attr(critV, "Calc") <- FALSE
  }
  if(pVal){
    pVals <- pValues(contMat, corMat, alpha, df,
                     tStat, alternative, mvtcontrol)
  }
  res <- list(contMat = contMat, corMat = corMat, tStat = tStat,
              alpha = alpha, alternative = alternative[1])
  if(pVal)
    attr(res$tStat, "pVal") <- pVals
  res$critVal <- critV
  class(res) <- "MCTtest"
  res
}

print.MCTtest <- function(x, digits = 3, eps = 1e-3, ...){
  cat("Multiple Contrast Test\n")
  cat("\n","Contrasts:","\n", sep="")
  print(round(x$contMat, digits))
  cat("\n","Contrast Correlation:","\n", sep="")
  print(round(x$corMat, digits))
  cat("\n","Multiple Contrast Test:","\n",sep="")
  ord <- rev(order(x$tStat))
  if(!any(is.null(attr(x$tStat, "pVal")))){
    pval <- format.pval(attr(x$tStat, "pVal"),
                        digits = digits, eps = eps)
    dfrm <- data.frame(round(x$tStat, digits)[ord],
                       pval[ord])
    names(dfrm) <- c("t-Stat", "adj-p")
  } else {
    dfrm <- data.frame(round(x$tStat, digits)[ord])
    names(dfrm) <- c("t-Stat")
  }
  print(dfrm)
  if(!is.null(x$critVal)){
    twoSide <- x$alternative == "two.sided"
    vec <- c(" one-sided)", " two-sided)")
    cat("\n","Critical value: ", round(x$critVal, digits), sep="")
    if(attr(x$critVal, "Calc")){
      cat(" (alpha = ", x$alpha,",", vec[twoSide+1], "\n", sep="")
    } else {
      cat("\n")
    }
  }
}

critVal <- function(corMat, alpha = 0.025, df = NULL,
                    alternative = c("one.sided", "two.sided"),
                    control = mvtnorm.control()){
  ## calculate critical value
  alternative <- match.arg(alternative)
  if(missing(corMat))
    stop("corMat needs to be specified")
  if(is.null(df))
    stop("degrees of freedom need to be specified")
  tail <- ifelse(alternative[1] == "two.sided",
                 "both.tails", "lower.tail")
  if (!missing(control)) {
    if(!is.list(control)) {
      stop("when specified, 'control' must be a list")
    }
    ctrl <- do.call("mvtnorm.control", control)
  } else {
    ctrl <- control
  }
  if(!is.finite(df)) # normal case
    df <- 0
  qmvtCall <- c(list(1-alpha, tail = tail, df = df, corr = corMat,
                algorithm = ctrl, interval = ctrl$interval))
  do.call("qmvt", qmvtCall)$quantile
}

checkAnalyArgs <- function(dose, resp, data, S, type,
                           addCovars, placAdj, na.action, cal){
  if(class(addCovars) != "formula")
    stop("addCovars argument needs to be of class formula")
  if(class(placAdj) != "logical")
    stop("placAdj argument needs to be of class logical")
  if(placAdj){
    if(type == "normal")
      stop("\"placAdj == TRUE\" only allowed for type = \"general\"")
    if(any(dose == 0))
      stop("If placAdj == TRUE there should be no placebo group")
  }
  if(!is.null(data)){ # data handed over in data frame
    if(!is.data.frame(data))
      stop("data argument needs to be a data frame")
    nams <- c(cal[2], cal[3], all.vars(addCovars))
    ind <- match(nams, names(data))
    if (any(is.na(ind)))
      stop("variable(s): ", paste(nams[is.na(ind)], collapse= ", "), " not found in ", cal[4])
    dd <- na.action(data[,nams])
  } else { # data handed over via vectors
    if(addCovars != ~1)
      stop("need to hand over data and covariates in data frame")
    if(length(dose) != length(resp))
      stop(cal[2], " and ", cal[3], " not of equal length")
    dd <- na.action(data.frame(dose, resp))
    cal[2:3] <- gsub("\\$", "", cal[2:3])
    cal[2:3] <- gsub("\\[|\\]", "", cal[2:3])
    colnames(dd) <- cal[2:3]
  }
  doseNam <- cal[2];respNam <- cal[3]
  if(any(dd[[doseNam]] < -.Machine$double.eps))
    stop("dose values need to be non-negative")
  if(!is.numeric(dd[[doseNam]]))
    stop("dose variable needs to be numeric")
  if(!is.numeric(dd[[respNam]]))
    stop("response variable needs to be numeric")
  ## check type related arguments
  if(type == "general" & is.null(S))
    stop("S argument missing")
  if(type == "normal" & !is.null(S))
    message("Message: S argument ignored for type == \"normal\"\n")
  if(type == "general" & addCovars != ~1)
    message("Message: addCovars argument ignored for type == \"general\"")
  if(!is.null(S)){
    if(!is.matrix(S))
      stop("S needs to be of class matrix")
    nD <- length(dd[[doseNam]])
    if(nrow(S) != nD | ncol(S) != nD)
      stop("S and dose have non-confirming size")
  }
  ord <- order(dd[[doseNam]])
  dd <- dd[ord, ]
  Sout <- NULL
  if(type == "general")
    Sout <- S[ord, ord]
  return(list(dd=dd, type = type, S = Sout, ord=ord,
              doseNam=doseNam, respNam=respNam))
}

pValues <- function(contMat, corMat, alpha = 0.025, df,
                    tStat, alternative = c("one.sided", "two.sided"),
                    control = mvtnorm.control()){
  ## function to calculate p-values
  nD <- nrow(contMat)
  nMod <- ncol(contMat)
  if(missing(corMat))
    stop("corMat needs to be specified")
  if(missing(df))
    stop("degrees of freedom need to be specified")
  if(length(tStat) != nMod)
    stop("tStat needs to have length equal to the number of models")
  alternative <- match.arg(alternative)
  ctrl <- mvtnorm.control()
  if (!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }
  if(!is.finite(df)) # normal case
    df <- 0
  lower <- switch(alternative[1],
                  one.sided = matrix(rep(-Inf, nMod^2), nrow = nMod),
                  two.sided = matrix(rep(-tStat, each = nMod), nrow = nMod))
  upper <- switch(alternative[1],
                  one.sided = matrix(rep(tStat, each = nMod), nrow = nMod),
                  two.sided = matrix(rep(tStat, each = nMod), nrow = nMod))
  pVals <- numeric(nMod)
  for(i in 1:nMod){
    pVals[i] <-   1 - pmvt(lower[,i], upper[,i], df = df,
                           corr = corMat, algorithm = ctrl)
  }
  pVals
}
