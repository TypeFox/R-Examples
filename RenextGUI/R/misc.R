##*********************************************************************
## Functions in this file are service function. Their arg are
## not GUI handler, and no reference to a widget must be found
##
##*********************************************************************



##=====================================================================
## Guess name from a candidate name and existing ones
## Note that the order of 'existing' has importance  
##=====================================================================

guessNames <- function(cand, existing = NULL) {
  x <- make.names(c(existing, cand), unique = TRUE)
  x[length(existing) + 1]
}

##=====================================================================
## round quantiles and confidence limitss using a 
## suitable number of digits
##=====================================================================

roundPred <- function(pred, dig.quant = NA) {

  cn <- colnames(pred)
  ## find the narrowest confint
  pct.L <- grep("L.[0-9]", cn, value = TRUE)
  pct.U <- grep("U.[0-9]", cn, value = TRUE)
  
  if (is.na(dig.quant) || length(dig.quant) == 0) {
    
    disp <- FALSE
    if (length(pct.L)) {
      i <- which.min(substring(pct.L, first = 3))
      vn1 <- pct.L[i]
      disp <- TRUE
    } else vn1 <- "quant"
    if (length(pct.U)) {
      i <- which.min(substring(pct.U, first = 3))
      vn2 <- pct.U[i]
      disp <- TRUE
    } else vn2 <- "quant"
    
    prec <-  mean(pred[ ,"quant"])/100
    if (disp)  prec <- pmin(prec, min(pred[ ,vn2] - pred[, vn1]))
    dig.quant <- -floor(log(prec, base = 10))
    if (dig.quant < 0) dig.quant <-  0 
  } 

  ## roudn all selected variables
  pred.mod <- pred
  ind <- c("quant", pct.L, pct.U)
  pred.mod[ , ind] <- round(pred.mod[ , ind], digits = dig.quant) 
  pred.mod
  
}

##======================================================================
## Coeff Summary
##======================================================================

makeSummary <- function(results) {

  est <- results$estimate
  s <- sqrt(diag(results$cov))
  t <- est/s
  prec <- s/10
  ind <- is.na(prec)
  if (any(ind)) prec[ind] <- est[ind]/1000
  dig.est <- -floor(log(prec, base = 10))
  dig.est[dig.est < 0] <-  0
  est <- round(est, digits = dig.est)
  s <- round(s, digits = dig.est)
  t <- round(t, digits = 3)
  data.frame(param. = I(names(est)),
             estim. = est,
             sd = s,
             t = t)

}

##======================================================================
## Find a reasonnable range of threshold from data values
##
##======================================================================

guessThreshold <- function(x) {
  
  rx <- range(x, na.rm = TRUE)
  qx <- quantile(x[!is.na(x)], prob = 0.95)
  
  step <- diff(pretty(rx))[1] / 50
  from <- floor(rx[1]/step)*step
  to <- ceiling(rx[2]/step)*step
  
  list(min = rx[1], max = rx[2], U = qx[1],
       from = from, to = to, by = step)
}

##==============================================================================
## Read demo datasets from Renext
##
##==============================================================================

makeDemos <- function( ) {
  ## initialisation to avoid a 'NOTE' when checking the packag
  Brest <- NULL;
  rm(Brest)         
  Garonne <- NULL
  rm(Garonne)     
  block <- 1
  Flow <- 1
  
  ## data(Brest)
  BrestData <- Brest$OTdata
  colnames(BrestData) <- c("date", "x")
  
  Brest.demo <- list(source = "demo",
                     csvPar = list(nCol = 2L, header = TRUE,
                       sep = ";", skip = 0L, dec = "."),
                     hasDate =  TRUE,
                     colDate = 1L,
                     dateFormat =  "%Y-%m-%d",
                     colNum = 2L,
                     main = Brest$info$shortLab,
                     info = Brest$info,
                     path = file.path(system.file("Rendata", package = "Renext"),
                       "Brest.csv"),
                     data = BrestData,
                     varName = Brest$info$varName,
                     effDuration = Brest$OTinfo$effDuration,
                     thresholds = guessThreshold(Brest$OTdata[ , 2]),
                     MAX = list(flag = FALSE),
                     rMax1 = "",
                     rMax1Dur = "",
                     rMax2 = "",
                     rMax2Dur = "")
  
  ## data(Garonne)
  
  GaronneData <- Garonne$OTdata
  colnames(GaronneData) <- c("date", "x")

  Gdata <- subset(Garonne$MAXdata, subset = (block == 1), select = Flow)$Flow
  
  ## keep block 1 only (in case other blocks would be added)
  rMax1 <-
    paste(subset(Garonne$MAXdata, subset = (block == 1), select = Flow)$Flow,
          collapse = ";")
  
  rMax1Dur <- Garonne$MAXinfo$duration[1]
  
  Garonne.demo <- list(source = "demo",
                       csvPar = list(nCol = 2L, header = TRUE,
                         sep = ";", skip = 0L, dec = "."),
                       hasDate =  TRUE,
                       colDate = 1L,
                       dateFormat =  "%Y-%m-%d",
                       colNum = 2L,
                       main = Garonne$info$shortLab,
                       info = Garonne$info,
                       path = file.path(system.file("Rendata", package = "Renext"),
                         "Garonne.csv"),
                       data = GaronneData,
                       varName = Garonne$info$varName,
                       effDuration = Garonne$OTinfo$effDuration,
                       thresholds = guessThreshold(Garonne$OTdata[ , 2]),
                       MAX = list(flag = TRUE,
                         block = factor(rep(1, length(Gdata)), levels = 1),
                         effDuration = rMax1Dur, 
                         r = length(Gdata),
                         data = Gdata),
                       rMax1 = rMax1,
                       rMax1Dur = rMax1Dur,
                       rMax2 = "",
                       rMax2Dur = "")

  list(Brest.demo = Brest.demo,
       Garonne.demo = Garonne.demo)
  
}

##======================================================================
## make a report with tables
##======================================================================

checkNumList <- function(text, min = -Inf, max = Inf) {

  nums <- as.numeric(unlist(strsplit(x = text, split = ";")))

  if ( any(is.na(nums)) ) stop("non-numeric value")
  if ( any(nums < min) ) stop("numeric value less than possible min")
  if ( any(nums > max) ) stop("numeric value greater than possible max")
  
  nums

}

##======================================================================
## abbreviate
##======================================================================

abbrNumList <- function(text) {

  if (length(text) && nchar(text)) {
    nums <- as.numeric(unlist(strsplit(x = text, split = ";")))
    if (length(nums) > 1) {
      rn <- range(nums)
      text <- sprintf("%d values, min = %s max = %s", length(nums),
                      format(rn[1]), format(rn[2]))
    } else {
      text <- sprintf("1 value: %s", format(nums))
    }
  } else {
    text <- "no value"
  }
  
  text

}

