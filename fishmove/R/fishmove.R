# Package Code 'fishmove'
# 
# Author: Johannes Radinger
###############################################################################


#Fishmove main default
fishmove <- function (species = NA, L = NA, AR = NA, SO = 6, T = 30, interval = "confidence", 
          rep = 50, seed = NA, ...) 
{
  if (missing(species) && missing(L)) 
    stop("No fish length or species provided")
  if (missing(species) && missing(AR)) 
    stop("No aspect ratio or species provided")
  if (!missing(species) && !is.element(species, speciesfishmove$SPECIES)) 
    stop("Fish species is not in speciesfishmove")
  if (!missing(species) && is.element(species, speciesfishmove$SPECIES) && 
        missing(L)) {
    L.tmp <- speciesfishmove$LENGTH[speciesfishmove$SPECIES == 
                                      species]
  }
  if (!missing(species) && is.element(species, speciesfishmove$SPECIES) && 
        missing(AR)) {
    AR.tmp <- speciesfishmove$ASPECT.RATIO[speciesfishmove$SPECIES == 
                                             species]
  }
  if (!missing(species) && is.element(species, speciesfishmove$SPECIES) && 
        !missing(L)) {
    warning("new.length will be used, not fish length from species-data", 
            call. = FALSE)
    L.tmp <- L
  }
  if (!missing(species) && is.element(species, speciesfishmove$SPECIES) && 
        !missing(AR)) {
    warning("new.aspect.ratio will be used, not aspect ratio from species", 
            call. = FALSE)
    AR.tmp <- AR
  }
  if (missing(species) && !missing(AR) && !missing(L)) {
    L.tmp <- L
    AR.tmp <- AR
  }
  if (min(datafishmove$LENGTH) > min(L.tmp) | max(L.tmp) > 
        max(datafishmove$LENGTH)) 
    warning("Fish length is outside the range of original regression", 
            call. = FALSE)
  if (min(datafishmove$ASPECT.RATIO) > min(AR.tmp) | max(AR.tmp) > 
        max(datafishmove$ASPECT.RATIO)) 
    warning("Aspect Ratio is outside the range of original regression", 
            call. = FALSE)
  if (min(datafishmove$STREAM.ORDER) > min(SO) | max(SO) > 
        max(datafishmove$STREAM.ORDER)) 
    warning("Stream order is outside the range of original regression", 
            call. = FALSE)
  if (min(datafishmove$TIME) > min(T) | max(T) > max(datafishmove$TIME)) 
    warning("Time is outside the range of original regression", 
            call. = FALSE)
  if (!is.numeric(rep)) 
    stop("rep is not numeric")
  if (rep > 5000) 
    stop("rep to large")
  if (!is.numeric(L.tmp)) 
    stop("L is not numeric")
  if (!is.numeric(AR.tmp)) 
    stop("AR is not numeric")
  if (!is.numeric(SO)) 
    stop("SO is not numeric")
  if (!is.numeric(T)) 
    stop("T is not numeric")
  if (interval != "prediction" && interval != "confidence") 
    stop("Interval must either be 'confidence' or 'prediction'")
  L <- L.tmp
  AR <- AR.tmp
  coefs = function(model) {
    c(coef(model), # statistical parameters
      summary(model.fishmove.stat)$coefficients[ ,2], # SE statistical paramters
      summary(model)$coefficients[-1, 4], # p-values
      pf(summary(model)$fstatistic[1], summary(model)$fstatistic[2], summary(model)$fstatistic[3], lower.tail = FALSE), #overall p-value
      summary(model)$r.squared # R-squared
    )
  }
  reg <- array(NA, dim = c(rep, 16, 2), dimnames = list(NULL, 
                                                        NULL, c("sigma_stat", "sigma_mob")))
  pred <- array(NA, dim = c(rep, 3, 2, length(L), length(AR), 
                            length(SO), length(T)), dimnames = list(NULL, NULL, c("sigma_stat", 
                                                                                  "sigma_mob"), paste("L", L, sep = "="), paste("AR", AR, 
                                                                                                                                sep = "="), paste("SO", SO, sep = "="), paste("T", T, 
                                                                                                                                                                              sep = "=")))
  if (!missing(seed)) {
    set.seed(seed)
  }
  seed_vect <- sample(1:10000, rep, replace = TRUE)
  for (i in 1:rep) {
    REP = NULL
    subsample_seed <- seed_vect[i]
    set.seed(subsample_seed)
    subsample <- ddply(datafishmove, .(REP), function(x) {
      x[sample(nrow(x), 1), ]
    })
    model.fishmove.stat <- lm(log(SIGMA_STAT) ~ log(LENGTH) + 
                                ASPECT.RATIO + sqrt(STREAM.ORDER) + log(TIME), data = subsample)
    model.fishmove.mob <- lm(log(SIGMA_MOB) ~ log(LENGTH) + 
                               ASPECT.RATIO + sqrt(STREAM.ORDER) + log(TIME), data = subsample)
    reg[i, , "sigma_stat"] <- coefs(model.fishmove.stat)
    reg[i, , "sigma_mob"] <- coefs(model.fishmove.mob)
    for (j in 1:length(L)) {
      for (k in 1:length(AR)) {
        for (l in 1:length(SO)) {
          for (m in 1:length(T)) {
            newdata <- data.frame(LENGTH = L[j], ASPECT.RATIO = AR[k], 
                                  STREAM.ORDER = SO[l], TIME = T[m])
            pred[i, , "sigma_stat", j, k, l, m] <- exp(predict.lm(model.fishmove.stat, 
                                                                  newdata = newdata, interval = interval))
            pred[i, , "sigma_mob", j, k, l, m] <- exp(predict.lm(model.fishmove.mob, 
                                                                 newdata = newdata, interval = interval))
          }
        }
      }
    }
  }
  coef.fishmove <- apply(reg, c(2, 3), mean, na.rm = TRUE)
  rownames(coef.fishmove) <- c("Intercept", "log(LENGTH)", "ASPECT.RATIO", "sqrt(STREAM.ORDER)", "log(TIME)",
                               "SE Intercept", "SE log(LENGTH)", "SE ASPECT.RATIO", "SE sqrt(STREAM.ORDER)", "SE log(TIME)",
                               "p-value log(LENGTH)", "p-value ASPECT.RATIO", "p-value sqrt(STREAM.ORDER)", "p-value log(TIME)",
                               "overall p value", "R.squared")
  pred.fishmove <- apply(pred, c(2, 3, 4, 5, 6, 7), mean, na.rm = TRUE)
  rownames(pred.fishmove) <- c("fit", "lwr", "upr")
  out <- list(coef.fishmove = coef.fishmove, pred.fishmove = pred.fishmove)
  class(out) <- "fishmove"
  return(out)
}

# Print method for fishmove 
print.fishmove=function(x,...){
	cat("Predicted movement for selected parameters:\n")
	print(x$pred.fishmove)
}



# Summary method
summary.fishmove <- function(object,...){
	cat("Summary\n")
	cat("Regression coefficients:\n")
	print(object$coef.fishmove)
	cat("Predicted movement for selected parameters:\n")
	print(object$pred.fishmove)
}

