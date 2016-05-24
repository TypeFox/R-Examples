#' @rdname ComparisonTable
#' @export conttable
#' 

'conttable' <- function(data, vars, byVar,
                 normal = NULL, var.equal = NULL, 
                 median = NULL,  
                 none = NULL,
                 odds = NULL, odds.scale=NULL, odds.unit=NULL,
                 alpha = 0.05, B=1000, seed=NULL){
                
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
    bad.vars <- vars[!vars %in% names(data)]
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

#******************************************************************************
#* Subroutine for information extraction
#* 1. Variable name, label, and level.  Level is for compatibilty with 
#*    cattable.
#* 2. Numeric Summary variables
#* 3. Odds Ratio
#* 4. Hypothesis Test and Subsequent Info
#* 5. 'type' is used to determine what is printed by print and write methods
#* 6. Build output
#******************************************************************************
  var.info <- function(v){

#*** 1. Variable name, label, and level.  
    nlev <- nlevels(data[, byVar])
    m.boot <- do.call("rbind",tapply(data[, v], data[, byVar],
                                     Hmisc::smean.cl.boot, B=B))
    quant <- do.call("rbind", tapply(data[, v], data[, byVar], stats::quantile,
                     probs=seq(0, 1, by=.25), na.rm=TRUE))
                     
#*** 2. Numeric Summary variables
    .name <- v
    .label <- if (Hmisc::label(data[,v]) %in% "") v else Hmisc::label(data[, v])
    .level <- NA
    .count <- colSums(table(data[, v], data[, byVar]))
    .prop <- rep(NA, nlev); names(.prop) <- levels(data[, byVar])
    .total <- sum(.count)
    .missing <- nrow(data) - .total
    .missingp <- .missing / nrow(data)
    .meanb <- m.boot[, 1]
    .lowerb <- m.boot[, 2]
    .upperb <- m.boot[, 3]
    .mean <- tapply(data[, v], data[, byVar], mean, na.rm=TRUE)
    .sd <- tapply(data[, v], data[, byVar], stats::sd, na.rm=TRUE)
    .min <- quant[, 1]
    .p25 <- quant[, 2]
    .median <- quant[, 3]
    .p75 <- quant[, 4]
    .max <- quant[, 5]
    
#*** 3. Odds Ratio
    if (v %in% odds & !0 %in% .count){
      if (nlev == 2){
        .odds.scale <- if (v %in% names(odds.scale)) odds.scale[[v]] else 1
        .odds.unit <- if (v %in% names(odds.unit))  odds.unit[[v]] else "units"
        m <- stats::glm(data[, byVar] ~ data[, v], family=stats::binomial)
        ci <- stats::confint(m, level=1 - alpha)
        .odds <- exp(stats::coef(m)[2] * .odds.scale)
        .odds.lower <- exp(ci[2,1] * .odds.scale)
        .odds.upper <- exp(ci[2,2] * .odds.scale)
      }
      else{
        warning(paste("Odds Ratio are only calulated when 'byVar' has exactly",
                      "two levels"))
        .odds <- .odds.lower <- .odds.upper <- .odds.scale <- .odds.unit <- NA
      }
    }
    else{
      .odds <- .odds.lower <- .odds.upper <- .odds.scale <- .odds.unit <- NA
    }


    EvalTable <- table(data[, v], data[, byVar])
    nlev.effective <- ncol(EvalTable) - sum(colSums(EvalTable) %in% 0)
    if ((nlev > 1 && nlev.effective %in% 1))
      warning(paste(v, ": grouping factor must have at least two levels.  No comparison is performed", sep=""))
    
#*** 4. Hypothesis Test and Subsequent Info
    if (nlev == 1 || nlev.effective %in% 1 || v %in% none){
      test.obj <- NULL
      test.obj$method <- NA
      test.obj$statistic <- NA
      test.obj$p.value <- NA
      .test.mark <- "N"
    }
    else if (nlev == 2){
      if (v %in% normal){
        v.eq <- v %in% var.equal
        test.obj <- stats::t.test(data[, v] ~ data[, byVar], var.equal=v.eq,
                           conf.level=1 - alpha)
        .test.mark <- "T"
      }
      else{
        warn <- withWarnings(stats::wilcox.test(data[, v] ~ data[, byVar]))
        if (!is.null(warn$warnings)) warning(v, ": ", warn$warnings)
        test.obj <- warn$value
        .test.mark <- "W"
      }
    }
    else{
      if (v %in% normal){
        test.obj <- stats::aov(data[, v] ~ data[, byVar])
        test.obj$method <- "Analysis of Variance"
        test.obj$statistic <- stats::anova(test.obj)[1, 4]
        test.obj$p.value <- stats::anova(test.obj)[1, 5]
        .test.mark <- "A"
      }
      else{
        test.obj <- stats::kruskal.test(data[, v] ~ data[, byVar])
        .test.mark <- "K"
      }
    }
    
#*** 5. 'type' is used to determine what is printed by print and write methods
    if (v %in% normal) .type <- "Parametric Mean"
    else if (v %in% median) .type <- "Median"
    else .type <- "Bootstrap Mean"

#*** 6. Build output
    names.df <- c("count", "prop", "boot", "lowerb", "upperb",
                  "mean", "sd", "min", "p25", "median", "p75", "max")
    names.df <- rep(names.df, rep(nlev, length(names.df)))
    names.df <- paste(names.df, levels(data[, byVar]), sep=".")

    .df <- c(.count, .prop, .meanb, .lowerb, .upperb, .mean,
             .sd,    .min,  .p25,   .median, .p75,    .max)
    .df <- as.data.frame(t(.df))
    .df <- cbind(.name, .label, .level, .total, .df, .missing, .missingp,
                 .odds, .odds.lower, .odds.upper, .odds.scale, .odds.unit,
                 test.obj$method, .test.mark, test.obj$statistic, 
                 test.obj$p.value,
                 is_significant(test.obj$p.value), .type,
                 stringsAsFactors=FALSE)
    names(.df) <- c("name", "label", "level", "total", names.df, "missing", "missing.perc",
                    "odds", "odds.lower", "odds.upper", "odds.scale",
                    "odds.unit", "test",
                    "test.mark", "test.stat", "pvalue", "significant", "type")
    rownames(.df) <- v
    return(.df)
  }
  
#*****************************************************************************
#* require Hmisc (for smean.cl.boot)
#* Set seed, if provided (recommended for reproducibility)
#* Convert byVar to a factor, if necessary
#* Send variables through var.info subroutine
#* Change class and return
#*****************************************************************************
  if (missing(byVar)){
    byVar <- "PlAcE_hOlDeR_fOr_CoNtTabLe"
    data[, byVar] <- factor("")
  }
  if (!is.null(seed)) set.seed(seed)
  if (!is.factor(data[, byVar])) data[, byVar] <- factor(data[, byVar])
#   if (!("ccf.df" %in% class(data))) data <- as.ccf.df.data.frame(data)
  ctable <- do.call("rbind", lapply(vars, var.info))
  ctable$type <- factor(ctable$type)
  class(ctable) <- c("ctable", "data.frame")
  attributes(ctable)$byVar <- data[, byVar]
  Hmisc::label(attributes(ctable)$byVar) <- Hmisc::label(data[, byVar])
  attributes(ctable)$vars <- vars  
  return(ctable)
}


