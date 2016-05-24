### likelihood ratio test for models of formal class glimML
if(!isGeneric("anova"))
  setGeneric(name = "anova", def = function(object, ...) standardGeneric("anova"))

setMethod(f = "anova", signature(object = "glimML"), definition = function(object, ...){

#### copied from anova method for lmer (packages Matrix / lme4)

  mCall <- match.call(expand.dots = TRUE)
  dots <- list(...)
  modp <- sapply(dots, inherits, "glimML")
  mods <- c(list(object), dots[modp])
  names(mods) <- sapply(as.list(mCall)[c(FALSE, TRUE, modp)], as.character)

####

  method <- mods[[1]]@method
  
# management of methods "BB" and "NB"
  is.meth <- unlist(lapply(mods, function(x){
    meth <- x@method
    meth == method
    }))
  mods <- mods[is.meth]
  n <- length(mods)
  if(n < 2)
    stop("At least 2 valid models are needed.")
  dfr <- data.frame("logL" = rep(NA, n), k = rep(NA, n), AIC = rep(NA, n), AICc = rep(NA, n), BIC = rep(NA, n),
                    "Resid. dev." = rep(NA, n), "Resid. Df" = rep(NA, n), Test = rep(NA, n), Deviance = rep(NA, n), 
                    Df = rep(NA, n), "P(> Chi2)" = rep(NA, n), check.names = FALSE)
  mod <- vector(mode = "character", length = length(mods))

# utility function to remove white spaces at the beginning and at the end of character strings
  tr <- function(string) gsub("^[[:space:]]+|[[:space:]]+$", "", string)

  for(i in seq(n)){
#    fm <- get(nam[i])
    fm <- mods[[i]]
# model name and formula
# paste + tr are necessary when there are interaction terms in the formula
    fb <- paste(tr(deparse(fm@formula)), collapse = " ")
    ft <- deparse(fm@random)
    mod[i] <- paste(names(mods)[i], ": ", "fixed = ", fb, "; random = ", ft, sep = "")
# anova table
    k <- fm@nbpar; akic <- AIC(fm)@istats
    dfres <- df.residual(fm)
    dfr$logL[i] <- fm@logL
    dfr$k[i]    <- k
    dfr$AIC[i]  <- akic[,2]
    dfr$AICc[i] <- akic[,3]
    dfr$BIC[i]  <- -2 * fm@logL + log(k + dfres) * k
    dfr[i, 6]   <- deviance(fm)
    dfr[i, 7]   <- dfres
    }
  dfr$Test <- c(NA, paste(names(mods)[-nrow(dfr)], names(mods)[-1], sep = "-"))
  dfr$Deviance  <- c(NA, 2 * diff(dfr[,1]))
  dfr$Df <- c(NA, diff(dfr[,2]))
  for(i in 2:n)
    dfr[i , 11] <- 1 - pchisq(abs(dfr$Deviance[i]), df = abs(dfr$Df[i]))
  rownames(dfr) <- names(mods)
  type <- switch(method, BB = "beta-binomial", NB = "negative-binomial")
  new(Class = "anova.glimML", models = mod, anova.table = dfr, type = type)
  })

setMethod("show", signature = "anova.glimML", definition = function(object){
  mod <- object@models
  dfr <- object@anova.table
  type <- object@type
  nam <- rownames(dfr)
  cat("Analysis of Deviance Table (", type, " models)\n\n", sep = "")
  sapply(mod, function(x) cat(x, "\n"))
  cat("\n")
  List <- lapply(dfr, function(x) ifelse(is.na(x), "", format(x, digits = 4)))
  dfr <- as.data.frame(t(do.call("rbind", List)))
  rownames(dfr) <- nam
  print(dfr)
  invisible(dfr)
  })
