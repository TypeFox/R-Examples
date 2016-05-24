#' @rdname ComparisonTable
#' @export cattable
#' 

cattable <- function(data, vars, byVar, fisher=NULL, fisher.arg=NULL,
                     cmh=NULL, row.score=NULL, col.score=NULL,
                     mcnemar=NULL, correct=NULL,
                     none=NULL,
                     odds=NULL, row.p=TRUE, alpha=0.05, minl=5){
                 
  if (missing(byVar)){
    byVar <- "PlAcE_hOlDeR_fOr_CaTcOnTtAbLe"
    data[, byVar] <- factor("")
  }
  
  withWarnings <- function(expr) {
  	myWarnings <- NULL
	  wHandler <- function(w) {
  	    myWarnings <<- c(myWarnings, list(paste(w$message, " - ", w$call[1])))
	      invokeRestart("muffleWarning")
  	}
	  val <- withCallingHandlers(expr, warning = wHandler)
  	list(value = val, warnings = myWarnings)
  }
  
  if (!all(vars %in% names(data))){
    bad.vars <- c(vars, byVar)[!c(vars, byVar) %in% names(data)]
    bad.vars.msg <- paste("The following variables are not found in 'data':", paste(bad.vars, collapse=", "))
    stop(bad.vars.msg)
  }
  
  all.missing <- sapply(data[, c(vars, byVar)], function(x) all(is.na(x)))
  if (any(all.missing)){
    miss.vars <- c(vars, byVar)[all.missing]
    miss.vars.msg <- paste("The following variables contain only missing values:", paste(miss.vars, collapse=", "))
    stop(miss.vars.msg)
  }
  
  if ("tbl_df" %in% class(data)) data <- as.data.frame(data)

  var.info <- function(v){
    if (!is.factor(data[, v])){
      v.lab <- Hmisc::label(data[, v])
      data[, v] <- factor(data[,v])
      Hmisc::label(data[, v]) <- v.lab
    }
    
    nlev <- nlevels(data[, byVar])
    nlev.v <- nlevels(data[, v])
    
    .name <- c(v, rep(NA, nlev.v))
    .label <- c(if (Hmisc::label(data[,v]) %in% "") v else Hmisc::label(data[, v]), rep(NA, nlev.v))
    .level <- c(NA, levels(data[, v]))
    .total <- c(sum(table(data[, v])), table(data[, v]))
    .count <- rbind(NA, table(data[, v], data[, byVar]))
    .missing <- c(.total[1] - sum(table(data[, v])), rep(NA, nlev.v))
    .missingp <- c(.missing[1] / nrow(data), rep(NA, nlev.v))
    .prop <- rbind(NA, prop.table(table(data[, v], data[, byVar]), 2 - row.p))
    .meanb <- matrix(NA, nrow=nlev.v + 1, ncol=nlevels(data[, byVar]))
    .lowerb <- matrix(NA, nrow=nlev.v + 1, ncol=nlevels(data[, byVar]))
    .upperb <- matrix(NA, nrow=nlev.v + 1, ncol=nlevels(data[, byVar]))
    .mean <- matrix(NA, nrow=nlev.v + 1, ncol=nlevels(data[, byVar]))
    .sd <- matrix(NA, nrow=nlev.v + 1, ncol=nlevels(data[, byVar]))
    .min <- matrix(NA, nrow=nlev.v + 1, ncol=nlevels(data[, byVar]))
    .p25 <- matrix(NA, nrow=nlev.v + 1, ncol=nlevels(data[, byVar]))
    .median <- matrix(NA, nrow=nlev.v + 1, ncol=nlevels(data[, byVar]))
    .p75 <- matrix(NA, nrow=nlev.v + 1, ncol=nlevels(data[, byVar]))
    .max <- matrix(NA, nrow=nlev.v + 1, ncol=nlevels(data[, byVar]))




    
    .odds.scale <- .odds.unit <- rep(NA, nlev.v + 1)
    if (v %in% odds){
      if (nlev == 2 & nlev.v > 1 & !0 %in% rowSums(.count)){
        warn <- withWarnings(stats::glm(data[, byVar] ~ data[, v], family=stats::binomial))
        if (!is.null(warn$warnings)) warning(v, "(glm): ", warn$warnings)
        m <- warn$value
        m$method <- "Logistic Regression"
        warn <- withWarnings(stats::confint(m, level=1 - alpha))
        if (!is.null(warn$warnings)) warning(paste(v, "(confint): ", 
                                                unique(unlist(warn$warnings)), 
                                                collapse="\n"))
        ci <- warn$value
        .odds <- c(NA, 1.0, exp(stats::coef(m)[-1]))
        .odds.lower <- c(NA, NA, exp(ci[-1,1]))
        .odds.upper <- c(NA, NA, exp(ci[-1,2]))
      }
      else{
        warning(paste("Odds Ratio are only calulated when 'byVar' has exactly",
                      "two levels AND 'var' has more than 1 level."))
        .odds <- .odds.lower <- .odds.upper <- .odds.scale <- .odds.unit <- NA
      }
    }
    else{
      .odds <- .odds.lower <- .odds.upper <- .odds.scale <- .odds.unit <- NA
    }
    
    EvalTable <- table(data[, v], data[, byVar])
    nlev.effective <- ncol(EvalTable) - sum(colSums(EvalTable) %in% 0)
    nlev.v.effective <- nrow(EvalTable) - sum(rowSums(EvalTable) %in% 0)
    
    if ((nlev > 1 && nlev.effective %in% 1) || nlev.v.effective %in% 1)
      warning(paste(v, ": 'x' and 'y' must have at least 2 levels.  No comparison will be performed'", sep=""))
    
#*** 4. Hypothesis Test and Subsequent Info
    if (nlev == 1 || nlev.v==1 || nlev.effective==1 || nlev.v.effective == 1 || v %in% none){
      test.obj <- NULL
      .test.method <- NA
      .test.stat <- NA
      .test.mark <- c("N", rep(NA, nlev.v))
      .pvalue <- NA
    }
    else{
      if (v %in% fisher){
        test.obj <- do.call("fisher.test", append(list(x=data[, v], y=data[, byVar], conf.level=1-alpha), fisher.arg))
#           fisher.test(data[, v], data[, byVar],
#                                 conf.level = 1 - alpha)
        .test.method <- c(test.obj$method, rep(NA, nlev.v))
        .test.mark <- c("F", rep(NA, nlev.v))
        .test.stat <- ifelse(is.null(test.obj$estimate), NA, test.obj$estimate)
        .test.stat <- c(.test.stat, rep(NA, nlev.v))
        .pvalue <- c(test.obj$p.value, rep(NA, nlev.v))
      }
      else if (v %in% cmh){
        if (!v %in% row.score){ 
          row.score <- append(row.score, list("equal"))
          names(row.score)[length(row.score)] <- v
        }
        if (!v %in% col.score){
          col.score <- append(col.score, list("equal"))
          names(col.score)[length(col.score)] <- v
        }
        test.obj <- mantel.test(data[, v], data[, byVar],
                                row.scores=row.score[[v]],
                                col.scores=col.score[[v]])
        .test.method <- c(test.obj$method, rep(NA, nlev.v))
        .test.mark <- c("MH", rep(NA, nlev.v))
        .test.stat <- ifelse(is.null(test.obj$correlation), NA, test.obj$correlation)
        .test.stat <- c(.test.stat, rep(NA, nlev.v))
        .pvalue <- c(test.obj$p.value, rep(NA, nlev.v))
      }
      else if(v %in% mcnemar){
        if (v %in% correct) cont <- TRUE else cont <- FALSE  
        test.obj <- stats::mcnemar.test(data[, v], data[, byVar], correct=cont)
        .test.method <- c(test.obj$method, rep(NA, nlev.v))
        .test.mark <- c("M", rep(NA, nlev.v))
        .test.stat <- ifelse(is.null(test.obj$estimate), NA, test.obj$estimate)
        .test.stat <- c(.test.stat, rep(NA, nlev.v))
        .pvalue <- c(test.obj$p.value, rep(NA, nlev.v))
      }
      else if (v %in% odds){
        if (nlev == 2){
          test.obj <- m
          .test.method <- c(NA, NA, rep(test.obj$method, nlev.v - 1))
          .test.mark <- c(NA, NA, rep("L", nlev.v - 1))
          .test.stat <- c(NA, 1.0, exp(stats::coef(test.obj)[-1]))
          .pvalue <- c(NA, NA, stats::coef(summary(test.obj))[-1, 4])
        }
        else{
          .test.method <- NA
          .test.mark <- c("N", rep(NA, nlev.v))
          .test.stat <- rep(NA, nlev.v + 1)
          .pvalue <- rep(NA, nlev.v + 1)
        }
      }
      else{
        warn <- withWarnings(stats::chisq.test(data[, v], data[, byVar]))
        if (!is.null(warn$warnings)) warning(v, ": ", warn$warnings)
        test.obj <- warn$value
        .test.method <- c(test.obj$method, rep(NA, nlev.v))
        .test.mark <- c("C", rep(NA, nlev.v))
        .test.stat <- c(test.obj$statistic, rep(NA, nlev.v))
        .pvalue <- c(test.obj$p.value, rep(NA, nlev.v))
      }
    }
    
    
    if (v %in% odds) .type <- "Logistic"
    else if (v %in% cmh) .type <- "Mantel-Haenszel"
    else if (v %in% fisher) .type <- "Fisher"
    else .type <- "Chi-Square"
    
    

    names.df <- c("count", "prop", "boot", "lowerb", "upperb",
                  "mean", "sd", "min", "p25", "median", "p75", "max")
    names.df <- rep(names.df, rep(nlev, length(names.df)))
    names.df <- paste(names.df, levels(data[, byVar]), sep=".")

    .df <- as.data.frame(cbind(.count, .prop, .meanb, .lowerb, .upperb,
                               .mean, .sd, .min, .p25, .median, .p75, .max))

    .df <- cbind(.name, .label, .level, .total, .df, .missing, .missingp,
                 .odds, .odds.lower, .odds.upper, .odds.scale, .odds.unit,
                 .test.method, .test.mark, .test.stat, .pvalue,
                 is_significant(.pvalue), .type,
                 stringsAsFactors=FALSE)
    names(.df) <- c("name", "label", "level", "total", names.df, "missing", "missing.perc",
                    "odds", "odds.lower", "odds.upper",
                    "odds.scale", "odds.unit",
                    "test", "test.mark", "test.stat",
                    "pvalue", "significant", "type")
                    
    rownames(.df) <- 
       c(v, paste(v, abbreviate(levels(data[, v]), minlength=minl), sep="-"))

    return(.df)
  }

  if (missing(byVar)){
    byVar <- "PlAcE_hOlDeR_fOr_CaTtAbLe"
    data[, byVar] <- factor("")
  }
#   if (!("ccf.df" %in% class(data))) data <- as.ccf.df.data.frame(data)

  if (!is.factor(data[, byVar])) data[, byVar] <- factor(data[, byVar])
  
  #toFactor <- vars[sapply(vars, function(x) !is.factor(data[, x]))]
  #data[, toFactor] <- lapply(data[, toFactor, drop=FALSE], factor)
  
  ctable <- do.call("rbind", lapply(vars, var.info))
  ctable$type <- factor(ctable$type)
  attributes(ctable)$byVar <- data[, byVar]
  Hmisc::label(attributes(ctable)$byVar) <- Hmisc::label(data[, byVar])
  attributes(ctable)$vars <- vars  
  class(ctable) <- c("ctable", "data.frame")
  return(ctable)
  
}
