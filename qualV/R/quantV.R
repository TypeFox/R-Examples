N <- function (o, p) {
  if (length(o) == 0 || length(p) == 0) stop ("vector of length zero")
  if (length(o) != length(p) && length(p) != 1) stop ("incompatible dimensions")
  length(o)
}

MeasureType <- function (dis, bigdis,
                         type = c("dissimilarity", "normalized",
                           "similarity", "reference", "formula")) {
  switch (EXPR = match.arg(type),
          dissimilarity = dis,
          normalized = dis / bigdis,
          similarity = 1 - dis / bigdis,
          reference = bigdis,
          formula = paste(deparse(substitute(dis)), "/", deparse(substitute(bigdis)))
         )
}

RelativeMeasureType <- function (dis, bigdis,
                                 type = c("dissimilarity", "normalized",
                                   "similarity", "reference", "formula")) {
  switch (EXPR = match.arg(type),
          dissimilarity = dis,
          normalized = dis / bigdis,
          similarity = 1 - log(dis / bigdis),
          reference = bigdis,
          formula = paste(deparse(substitute(dis)), "/", deparse(substitute(bigdis)))
         )
}

MAE  <- function (o, p, type = "dissimilarity") {
  MeasureType(
    (sum(abs(o - p))) / N(o, p),
    mean(abs(o - median(o))), type
  )
}

MAPE <- function (o, p, type = "dissimilarity") {
  m <- median(o)
  MeasureType(
    100 * (sum(abs((o - p) / o))) / N(o, p),
    100 * (sum(abs((o - m) / m))) / N(o, p),
    type
  )
}

MSE  <- function (o, p, type = "dissimilarity") {
  MeasureType((sum ((o - p)^2)) / N(o, p), var(o), type)
}
  
RMSE <- function (o, p, type = "dissimilarity") {
  MeasureType(sqrt(MSE(o, p)), sd(o), type)
}

CMAE <- function (o, p, type = "dissimilarity") {
  MeasureType(
    sum(abs(o - p - median(o - p))) / N(o, p),
    mean(abs(o - median(o))), type
  )
}

CMSE <- function (o, p, type = "dissimilarity") {
  MeasureType(sum((o - p - mean(o - p))^2) / (N(o, p) - 1), var(o), type)
}

RCMSE <- function (o, p, type = "dissimilarity") {
  MeasureType(sqrt(CMSE(o, p)), sd(o), type)
}

SMAE <- function(o, p, type = "dissimilarity") {
  N(o, p)
  r <- resid(lm(o ~ p))
  erg <- sum(abs(r)) / (length(r))
  MeasureType(erg, mean(abs(o - median(o))), type)  
}

SMSE <- function (o, p, type = "dissimilarity") {
  N(o, p)
  r <- resid(lm(o ~ p))
  erg <- sum(r^2) / (length(r) - 2)
  MeasureType(erg, var(o), type)
}

RSMSE <- function (o, p, type = "dissimilarity") {
  N(o, p)
  r <- resid(lm(o ~ p))
  erg <- sqrt(sum(r^2) / (length(r) - 2))
  MeasureType(erg, sd(o), type)
}

MALE <- function (o, p, type = "dissimilarity") {
  MeasureType(
    sum(abs(log(o / p))) / N(o, p),
    mean(abs(o - median(o))), type
  )
}
    
MAGE <- function (o, p, type = "dissimilarity") {
  RelativeMeasureType(
    exp(MALE(o, p)),
    exp(mean(abs(log(o) - median(log(o))))),
    type
  )
}

MSLE <- function (o, p, type = "dissimilarity") {
  MeasureType( sum(log(o / p)^2) / N(o, p), var(log(o)), type)
}


RMSLE <- function (o, p, type = "dissimilarity") {
  MeasureType(sqrt(MSLE(o, p)), sd(log(o)), type)
}

RMSGE <- function (o, p, type = "dissimilarity") {
  RelativeMeasureType(exp(RMSLE(o, p)), exp(sd(log(o))), type)
}

SMALE <- function (o, p, type = "dissimilarity") {
  MeasureType(
    sum(abs(log(o / p) - median(log(o / p)))) / N(o, p),
    mean(abs(log(o) - median(log(o)))), type
  )
}

SMAGE <- function (o, p, type = "dissimilarity") {
  RelativeMeasureType(
    exp(SMALE(o, p)),
    exp(mean(abs(log(o) - median(log(o))))),
    type
  )
}

SMSLE <- function (o, p, type = "dissimilarity") {
  MeasureType(
    sum((log(o / p) - mean(log(o / p)))^2) / (N(o, p) - 1),
    var(log(o)), type
  )
}

RSMSLE <- function (o, p, type = "dissimilarity") {
  MeasureType(sqrt(SMSLE(o, p)), sd(log(o)), type)
}

RSMSGE <- function (o, p, type = "dissimilarity") {
  RelativeMeasureType(exp(RSMSLE(o, p)), exp(sd(log(o))), type)
}

rankportion <- function (x) (rank(x) - 1) / (length(x) - 1)

MAOE <- function (o, p, type = "dissimilarity") {
  MeasureType(
    sum(abs(rankportion(o) - rankportion(p))) / N(o, p),
    0.25, type
  )
}

MSOE <- function (o, p, type = "dissimilarity") {
  MeasureType(
    sum((rankportion(o) - rankportion(p))^2) / N(o, p),
    0.25, type
  )
}

RMSOE <- function (o, p, type = "dissimilarity") {
  MeasureType(sqrt(MSOE(o, p)), 0.5, type)
}

generalME <- function (o, p, ignore = c("raw", "centered", "scaled", "ordered"),
                       geometry = c("real", "logarithmic", "geometric", "ordinal"),
                       measure = c("mad", "var", "sd"),
                       type = c("dissimilarity", "normalized", "similarity",
                         "reference", "formula", "name", "function"), method = NULL) {

  geometry <- match.arg(geometry)
  type     <- match.arg(type)
  measure  <- match.arg(measure)
  ignore   <- match.arg(ignore)
  if (is.null(method)) {
    ERR <- function (o, p, ...) {
      stop ("The requested deviance measusure does not exist")
    }
    method <- switch (geometry,
                      real = switch (ignore,
                        raw = list(mad = "MAE", var = "MSE", sd = "RMSE"),
                        centered = list(mad = "CMAE", var = "CMSE", sd = "RCMSE"),
                        scaled = list(mad = "SMAE", var = "SMSE", sd = "RSMSE"),
                        ordered = list(mad = "MAOE", var = "MSOE", sd = "RMSOE")
                        ),
                      logarithmic = switch(ignore,
                        raw = list(mad = "MALE", var = "MSLE", sd = "RMSLE"),
                        centered = list(mad = "SMALE", var = "SMSLE", sd = "RSMSLE"),
                        scaled = list(mad = "SMALE", var = "SMSLE", sd = "RSMSLE"),
                        ordered = list(mad = "MAOE", var = "MSOE", sd = "RMSOE")
                        ),
                      geometric = switch(ignore,
                        raw = list(mad = "MAGE", var = "RMSGE", sd = "RMSGE"),
                        centered = list(mad = "SMAGE", var = "RSMSGE", sd = "RSMSGE"),
                        scaled = list(mad = "SMAGE", var = "RSMSGE", sd = "RSMSGE"),
                        ordered = list(mad = "MAOE", var = "MSOE", sd = "RMSOE")
                        ),
                      ordinal = list(mad = "MAOE", var = "MSOE", sd = "RMSOE")
                     )[[measure]]
  }
  if (type == "name") {
    if (is.character(method))
      return (method)
    else {
      for (x in objects(pos = "package:qualV"))
        if(identical(method, get(x, pos = "package:qualV")))
          return (x)
      stop(paste("method ", deparse(method), "not found", col = " "))  
    }
  }
  if (is.character(method))
    method <- get(method, envir = environment(MSE), mode = "function")
  if (type == "function")
    return (method)
  method(o, p, type)
}
