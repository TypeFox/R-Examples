anova.aodml <- function(object, ...) {

  #### copied from anova method for lmer (packages Matrix / lme4)

  mCall <- match.call(expand.dots = TRUE)
  dots <- list(...)
  modp <- sapply(dots, inherits, "aodml")
  mods <- c(list(object), dots[modp])
  names(mods) <- sapply(as.list(mCall)[c(FALSE, TRUE, modp)], as.character)

	####
	
	n <- length(mods)
  if(n < 2)
    stop("At least 2 valid models are needed.")
  
	dfr <- data.frame(
    df.resid = rep(NA, n), df.model = rep(NA, n), "2 x Log-lik" = rep(NA, n),
  	AIC = rep(NA, n), AICc = rep(NA, n), BIC = rep(NA, n),
  	Test = rep(NA, n), df = rep(NA, n), "LR stat." = rep(NA, n), 
  	"P(>LR stat.)" = rep(NA, n), check.names = FALSE
    )
	
	mod <- vector(mode = "character", length = length(mods))
	
	# utility function to remove white spaces at the beginning and at the end of character strings
  tr <- function(string) gsub("^[[:space:]]+|[[:space:]]+$", "", string)

  for(i in seq(n)){
  	#    fm <- get(nam[i])
    fm <- mods[[i]]
  	
    # model name and formula
  	# paste + tr are necessary when there are interaction terms in the formula
    fb <- paste(tr(deparse(fm$formula)), collapse = " ")
    ft <- deparse(fm$phi.formula)
    mod[i] <- paste(names(mods)[i], ": ", "Mu = ", fb, " ; Phi = ", ft, " ($phi = ", round(fm$phi, 4), ")", sep = "")
  	
    # anova table
    akic <- AIC(fm)
    df.resid <- df.residual(fm)
    df.model <- fm$df.model
    
    dfr$df.model[i] <- df.model
    dfr$df.resid[i] <- df.resid 
    dfr$AIC[i]  <- akic[, 3]
    dfr$AICc[i] <- akic[, 4]
    dfr$BIC[i]  <- -2 * fm$logL + log(df.model + df.resid) * df.model

    dfr[i, 3] <- 2 * fm$logL
  }
  
	dfr$Test <- c(NA, paste(names(mods)[-nrow(dfr)], names(mods)[-1], sep = " vs "))
  dfr$df <- c(NA, abs(diff(dfr$df.model)))
  dfr[, 9]  <- c(NA, abs(diff(dfr[, 3])))
  
	for(i in 2:n)
    dfr[i , 10] <- 1 - pchisq(dfr[i, 9], df = dfr$df[i])
  rownames(dfr) <- names(mods)
	
	structure(list(models = mod, anova.table = dfr), class = "anova.aodml")

}
