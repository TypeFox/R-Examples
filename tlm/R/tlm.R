tlm <-
function(y, x, z = "none", family = gaussian, data, ypow = 1, xpow = 1, ...)
 {
  datal <- deparse(substitute(data))
  ################ family gaussian, binomial or poisson:
  if (is.character(family))
   {
   	responses <- c("gaussian", "binomial", "poisson")
    familymatch <- charmatch(family, responses, nomatch = 0)
    if (familymatch == 0)
     stop("'family' not allowed")
    family <- get(responses[familymatch], mode = "function", envir = parent.frame())
   }
  if (is.function(family)) family <- family()
  
  ################# 'ypow' control (in addition, family binomial or poisson --> ypow = 1):
  if (!inherits(ypow, c("numeric", "integer")) || length(ypow) != 1L)
    stop("'ypow' must be a number")
    
  if ((ypow != 1) && (family$family != "gaussian"))
   {
   	warning("assuming 'ypow' equal to 1 because the response is not gaussian")
    ypow <- 1
   }

  ################# 'xpow' control:
  if (!inherits(xpow, c("numeric", "integer")) || length(xpow) != 1L)
  	stop("'xpow' must be a number")
  	
  ################ y existence:
  yl <- deparse(substitute(y))
  if (!yl %in% names(data))
    stop(paste("variable '", yl, "' not found", sep = ""))
  
  ################ x existence:
  xl0 <- deparse(substitute(x))
  if (!xl0 %in% names(data))
   	stop(paste("variable '", xl0, "' not found", sep = ""))
  
 ################# x to factor
  xaux <- data[, xl0]
  if (inherits(xaux, c("character", "logical")) || (inherits(xaux, c("numeric", "integer")) && length(unique(xaux[!is.na(xaux)])) == 2))
   data[, xl0] <- as.factor(data[, xl0])

  if ((xpow != 1) && inherits(data[, xl0], "factor"))
   {
   	warning("assuming 'xpow' equal to 1 because x is categorical")
    xpow <- 1
   }

 ################# y to factor
  yaux <- data[, yl]
  if (inherits(yaux, c("character", "logical")) || (inherits(yaux, c("numeric", "integer")) && length(unique(yaux[!is.na(yaux)])) == 2))
   data[, yl] <- as.factor(data[, yl])
   
  if (inherits(data[, yl], "factor") && length(levels(data[, yl])) > 2)
   	stop("categorical response variable with more than 2 levels is not allowed")
  
  if (inherits(data[, yl], "factor") && family$family != "binomial")
   	stop("categorical response variable is only allowed under 'binomial' family")

  if (ypow != 1 && inherits(data[, yl], "factor"))
   {
   	warning("assuming 'ypow' equal to 1 because y is categorical")
    xpow <- 1
   }
  
  ################# bijectivity control
  
  xbij <- T
  if (!inherits(xaux, "factor") && xpow != 1) xbij <- bijectivityback(xt = xaux, power = xpow)
  if (!xbij)
   	stop("non bijective 'x': the provided value for 'xpow' is incompatible with the provided values of 'x'")
  
  ybij <- T
  if (!inherits(yaux, "factor") && ypow != 1) ybij <- bijectivityback(xt = yaux, power = ypow)
  if (!ybij)
   	stop("non bijective 'y': the provided value for 'ypow' is incompatible with the provided values of 'y'")
   	
  
  ################ power labels for plot and print:
  ypowlabel <- deparse(substitute(ypow))
  xpowlabel <- deparse(substitute(xpow))

  ################ formula building:
  # if (!is.factor(data[, xl0]) && length(unique(data[, xl0][!is.na(data[, xl0])])) == 2)
   # xl <- paste("as.factor(", xl0, ")", sep = "") else xl <- xl0
  xl <- xl0
  yxl <- paste(yl, xl, sep = " ~ ")
  zl <- deparse(substitute(z), width.cutoff = 500L)
  if (zl == "\"none\"")
   formula <- yxl else
   {
  	if (zl == "." & ncol(data) > 2)
  	 zl <- paste(names(data)[!names(data) %in% c(yl, xl)], collapse = " + ")
  	formula <- paste(yxl, zl, sep = " + ")
   }
  ################ fitting model in the transformed space:
  #datal <- deparse(substitute(data))
  mod <- NULL
  if (family$family == "gaussian")
    modtext <- paste("mod <- lm(formula = ", formula, ", data = ", datal, ", ...)", sep = "") else
    modtext <- paste("mod <- glm(formula = ", formula, ", family = ", family$family, ", data = ", datal, ", ...)", sep = "")
  eval(parse(text = modtext))
  ################ results:
  mf <- model.frame(mod)
  mt <- attr(mf, "terms")
  Xclass <- attr(mt, "dataClasses")[2]
  if (Xclass == "factor" && xpow != 1)
   {
   	warning("assuming 'xpow' = 1 because the explanatory variable '", xl0, "' is a factor", sep = "")
   	xpow <- 1
   }
  res <- list(model = mod, ypow = ypow, xpow = xpow)
  attr(res, "ypowlabel") <- ypowlabel
  attr(res, "xpowlabel") <- xpowlabel
  class(res) <- "tlm"
  res
 }
