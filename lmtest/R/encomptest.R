## See Greene (2003), Section 8.3.2, p.154
## original reference: Davidson-MacKinnon 1993, p.384-7 ex Greene

encomptest <- function(formula1, formula2, data = list(), vcov. = NULL, ...)
{
  ## merge two models (if possible) and compute
  ## response y and regressor matrices X and Z
  ## 1. If formulas are available: build big model first
  if(inherits(formula2, "formula") && inherits(formula2, "formula")) {
    formula <- formula2
    if(length(formula) > 2) formula[[2]] <- NULL
    formula[[2]] <- as.call(list(as.name("+"), as.name("."), formula[[2]]))
    formula <- update(formula1, formula)
    mf <- model.frame(formula, data = data, ...)

    ## cannot compute directly on mf because of (potential)
    ## transformations in the variables, try matching via row.names.
    ## formula <- formula1
    ## formula[[3]] <- as.name(".") 
    if(is.data.frame(data)) data <- data[row.names(mf),]
       
    fm1 <- lm(formula1, data = data)
    fm2 <- lm(formula2, data = data)
    fm <- lm(formula, data = data)
  } else {
    fm1 <- if(inherits(formula1, "formula")) lm(formula1, data = data, ...) else formula1
    if(inherits(formula2, "formula")) {
      fm2 <- update(fm1, formula2)
    } else {
      fm2 <- formula2
      if(deparse(terms(fm1)[[2]]) != deparse(terms(fm2)[[2]]))
        stop("models have different response variables")
      if(length(fm1$residuals) != length(fm2$residuals))
        stop("models were not all fitted to the same size of dataset")
    }  
    
    ## compute encompassing model
    fm <- attr(terms(fm2), "term.labels")
    fm <- fm[!(fm %in% attr(terms(fm1), "term.labels"))]
    fm <- paste(fm, collapse = " + ")
    fm <- as.formula(paste(". ~ . +", fm))
    fm <- update(fm1, fm)
  }
  m1 <- paste("Model 1:", paste(deparse(formula(fm1)), collapse = "\n"))
  m2 <- paste("Model 2:", paste(deparse(formula(fm2)), collapse = "\n"))
  m  <- paste("Model E:", paste(deparse(formula(fm )), collapse = "\n"))

  ## check vcov.
  if(!is.null(vcov.) && !is.function(vcov.)) stop("`vcov.' needs to be a function")

  rval1 <- as.data.frame(waldtest(fm, fm1, vcov = vcov., ...))
  rval2 <- as.data.frame(waldtest(fm, fm2, vcov = vcov., ...))
  rval <- rbind(rval1[2,], rval2[2,])
  rval[,1] <- c(rval1[1,1], rval2[1,1])
  rownames(rval) <- c("M1 vs. ME", "M2 vs. ME")

  ## put results together
  title <- "Encompassing test\n"
  topnote <- paste(c(m1, m2, m), collapse="\n")
  rval <- structure(rval, heading = c(title, topnote), class = c("anova", "data.frame"))
  return(rval)
}
