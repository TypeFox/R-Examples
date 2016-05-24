################################################################################
# corRStruct method functions
################################################################################

Initialize.corRStruct <- function(object, data, ...)
{
   form <- formula(object)
   ## obtaining the groups information, if any
   if (!is.null(getGroupsFormula(form))) {
      attr(object, "groups") <- getGroups(object, form, data = data)
      attr(object, "Dim") <- Dim(object, attr(object, "groups"))
   } else { # no groups
      attr(object, "Dim") <- Dim(object, as.factor(rep(1, nrow(data))))
   }
   ## obtaining the covariate(s)
   attr(object, "covariate") <- getCovariate(object, data = data)

   object
}


Dim.corRStruct <- function(object, groups, ...)
{
   if (missing(groups)) return(attr(object, "Dim"))
   ugrp <- unique(groups)
   groups <- factor(groups, levels = ugrp)
   len <- table(groups)

   list(N = length(groups),
        M = length(len),
        maxLen = max(len),
        sumLenSq = sum(len^2),
        len = len,
        start = match(ugrp, groups) - 1)
}


print.corRStruct <- function(x, ...)
{
   aux <- coef(x)
   if (length(aux) > 0) {
      cat("Correlation structure of class", class(x)[1], "representing\n")
      print(invisible(aux), ...)
   } else {
      cat("Uninitialized correlation structure of class", class(x)[1], "\n")
   }
}


################################################################################
# corRSpatial method functions
################################################################################

Initialize.corRSpatial <- function(object, data, ...)
{
   if (!is.null(attr(object, "covariate"))) { # already initialized
      return(object)
   }
   object <- Initialize.corRStruct(object, data)

   val <- as.vector(object)
   if (length(val) == 0) {
      val <- attr(getCovariate(object), "minD") * 0.9
   } else if (!all(inbounds(val, attr(object, "bounds")))) {
      stop()
   }
   attributes(val) <- attributes(object)

   val
}


Dim.corRSpatial <- function(object, groups, ...)
{
   if (missing(groups)) return(attr(object, "Dim"))
   val <- Dim.corRStruct(object, groups)
   val[["start"]] <-
      c(0, cumsum(val[["len"]] * (val[["len"]] - 1)/2)[-val[["M"]]])
   ## will use third component of Dim list for spClass
   names(val)[3] <- "spClass"
   val[[3]] <- match(class(object)[1], c("corRExp", "corRExpwr", "corRGaus",
                     "corRGneit", "corRLin", "corRMatern", "corRCauchy",
                     "corRSpher"), 0)

   val
}


getCovariate.corRSpatial <- function(object, form = formula(object), data)
{
   covar <- attr(object, "covariate")

   if (is.null(covar)) { # need to calculate it
      if (missing(data)) {
         stop("Need data to calculate covariate")
      }
      covForm <- terms(getCovariateFormula(form))
      attr(covForm, "intercept") <- 0
      if (length(all.vars(covForm)) > 0) { # covariate present
         covar <- model.matrix(covForm,
                     model.frame(covForm, data, drop.unused.levels = TRUE))
      } else {
         covar <- as.matrix(1:nrow(data))
      }

      if (is.null(getGroupsFormula(form))) { # no groups
         attr(covar, "assign") <- NULL
         attr(covar, "contrasts") <- NULL
         attr(covar, "dist") <- as.vector(dist2(covar,
                                   method = attr(object, "metric"),
                                   r = attr(object, "radius")))
         attr(covar, "minD") <- min(attr(covar, "dist"))
      } else { # by groups
         grps <- getGroups(object, data = data)
         covar <- lapply(split(as.data.frame(covar), grps),
            function(el, metric, radius) {
               el <- as.matrix(el)
               attr(el, "dist") <- as.vector(dist2(el, metric, r = radius))
               el
            }, metric = attr(object, "metric"), radius = attr(object, "radius"))
         attr(covar, "minD") <- min(unlist(lapply(covar, attr, which = "dist")))
      }

      if (attr(covar, "minD") == 0) {
         stop("Cannot have zero distances in \"corRSpatial\"")
      }
   }

   covar
}


corMatrix.corRSpatial <- function(object, covariate = getCovariate(object),
   corr = TRUE, ...)
{
   if (data.class(covariate) == "list") {
      dist <- unlist(lapply(covariate, attr, which = "dist"))
      len <- unlist(lapply(covariate, nrow))
   } else {
      dist <- attr(covariate, "dist")
      len <- nrow(covariate)
      names(len) <- 1
   }

   par <- coef(object)
   val <- switch(class(object)[1],
      corRExp    = cor.exp(dist, par[1]),
      corRExpwr  = cor.exp(dist, par[1], par[2]),
      corRGaus   = cor.exp(dist, par[1], 2),
      corRGneit  = cor.gneiting(dist, par[1]),
      corRLin    = cor.lin(dist, par[1]),
      corRMatern = cor.matern(dist, par[1], par[2]),
      corRCauchy = cor.cauchy(dist, par[1]),
      corRSpher  = cor.spher(dist, par[1]),
      corRWave   = cor.wave(dist, par[1])
   )

   val <- split(val, rep(names(len), len * (len - 1) / 2))
   lD <- NULL
   for(i in names(val)) {
      x <- matrix(0, len[i], len[i])
      x[lower.tri(x)] <- val[[i]]
      if (corr) {
         val[[i]] <- x + t(x)
         diag(val[[i]]) <- 1
      } else {
         diag(x) <- 1
         l <- chol(t(x))
         val[[i]] <- t(backsolve(l, diag(len[i])))
         lD <- c(lD, diag(l))
      }
   }
   if (length(len) == 1) val <- val[[1]]
   if (!is.null(lD)) attr(val, "logDet") <- -1 * sum(log(lD))

   val
}


corFactor.corRSpatial <- function(object, ...)
{
   val <- corMatrix(object, corr = FALSE, ...)
   lD <- attr(val, "logDet")
   if (is.list(val)) val <- unlist(val)
   else val <- as.vector(val)
   names(val) <- NULL
   attr(val, "logDet") <- lD

   val
}


coef.corRSpatial <- function(object, ...)
{
   val <- as.vector(object)
   if (length(val) == 0) {
      return(val)
   }
   names(val) <- rownames(attr(object, "bounds"))

   val
}


"coef<-.corRSpatial" <- function(object, ..., value)
{
   if (!all(inbounds(value, attr(object, "bounds")))) stop()
   object[] <- value

   object
}


################################################################################
# corRExp - exponential spatial correlation structure
################################################################################

corRExp <- function(value = numeric(0), form = ~ 1,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956)
{
   attr(value, "formula") <- form
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "bounds") <- matrix(c(0, Inf, 1), ncol=3,
      dimnames = list("range", c("lower", "upper", "type")))
   class(value) <- c("corRExp", "corRSpatial", "corRStruct")

   value
}


################################################################################
# corRExpwr - Powered exponential spatial correlation structure
################################################################################

corRExpwr <- function(value = numeric(0), form = ~ 1,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956)
{
   attr(value, "formula") <- form
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "bounds") <- matrix(c(0, 0, Inf, 2, 1, 3), ncol=3,
      dimnames = list(c("range", "shape"), c("lower", "upper", "type")))
   class(value) <- c("corRExpwr", "corRSpatial", "corRStruct")

   value
}


Initialize.corRExpwr <- function(object, data, ...)
{
   if (!is.null(attr(object, "covariate"))) { # already initialized
      return(object)
   }
   object <- Initialize.corRStruct(object, data)

   val <- as.vector(object)
   if (length(val) == 0) {
      val <- c(attr(getCovariate(object), "minD") * 0.9, 1)
   } else if (!all(inbounds(val, attr(object, "bounds")))) {
      stop()
   }
   attributes(val) <- attributes(object)

   val
}


################################################################################
#corRGaus - Gaussian spatial correlation structure
################################################################################

corRGaus <- function(value = numeric(0), form = ~ 1,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956)
{
   attr(value, "formula") <- form
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "bounds") <- matrix(c(0, Inf, 1), ncol=3,
      dimnames = list("range", c("lower", "upper", "type")))
   class(value) <- c("corRGaus", "corRSpatial", "corRStruct")
   value
}


################################################################################
#corRGneit - Gneiting spatial correlation structure
################################################################################

corRGneit <- function(value = numeric(0), form = ~ 1,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956)
{
   attr(value, "formula") <- form
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "bounds") <- matrix(c(0, Inf, 1), ncol=3,
      dimnames = list("range", c("lower", "upper", "type")))
   class(value) <- c("corRGneit", "corRSpatial", "corRStruct")
   value
}


################################################################################
# corRLin - Linear spatial correlation structure
################################################################################

corRLin <- function(value = numeric(0), form = ~ 1,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956)
{
   attr(value, "formula") <- form
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "bounds") <- matrix(c(0, Inf, 1), ncol=3,
      dimnames = list("range", c("lower", "upper", "type")))
   class(value) <- c("corRLin", "corRSpatial", "corRStruct")
   value
}


################################################################################
# corRMatern - Matern spatial correlation structure
################################################################################

corRMatern <- function(value = numeric(0), form = ~ 1,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956)
{
   attr(value, "formula") <- form
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "bounds") <- matrix(c(0, 0, Inf, 2, 1, 3), ncol=3,
      dimnames = list(c("range", "scale"), c("lower", "upper", "type")))
   class(value) <- c("corRMatern", "corRSpatial", "corRStruct")
   value
}


Initialize.corRMatern <- function(object, data, ...)
{
   if (!is.null(attr(object, "covariate"))) { # already initialized
      return(object)
   }
   object <- Initialize.corRStruct(object, data)

   val <- as.vector(object)
   if (length(val) == 0) {
      val <- c(attr(getCovariate(object), "minD") * 0.9, 0.5)
   } else if (!all(inbounds(val, attr(object, "bounds")))) {
      stop()
   }
   attributes(val) <- attributes(object)

   val
}


################################################################################
# corRCauchy - Cauchy spatial correlation structure
################################################################################

corRCauchy <- function(value = numeric(0), form = ~ 1,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956)
{
   attr(value, "formula") <- form
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "bounds") <- matrix(c(0, Inf, 1), ncol=3,
      dimnames = list("range", c("lower", "upper", "type")))
   class(value) <- c("corRCauchy", "corRSpatial", "corRStruct")
   value
}


################################################################################
# corRSpher - spherical spatial correlation structure
################################################################################

corRSpher <- function(value = numeric(0), form = ~ 1,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956)
{
   attr(value, "formula") <- form
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "bounds") <- matrix(c(0, Inf, 1), ncol=3,
      dimnames = list("range", c("lower", "upper", "type")))
   class(value) <- c("corRSpher", "corRSpatial", "corRStruct")
   value
}


################################################################################
# corRWave - sine wave spatial correlation structure
################################################################################

corRWave <- function(value = numeric(0), form = ~ 1,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956)
{
   attr(value, "formula") <- form
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "bounds") <- matrix(c(0, Inf, 1), ncol=3,
      dimnames = list("range", c("lower", "upper", "type")))
   class(value) <- c("corRWave", "corRSpatial", "corRStruct")
   value
}


################################################################################
# corRSpatioTemporal method functions
################################################################################

Initialize.corRSpatioTemporal <- function(object, data, ...)
{
   if (!is.null(attr(object, "covariate"))) { # already initialized
      return(object)
   }
   object <- Initialize.corRStruct(object, data)

   val <- as.vector(object)
   if (length(val) == 0) {
      val <- attr(getCovariate(object), "minD") * 0.9
      val[val == 0] <- 1
   } else if (!all(inbounds(val, attr(object, "bounds")))) {
      stop()
   }
   attributes(val) <- attributes(object)

   val
}


Dim.corRSpatioTemporal <- function(object, groups, ...)
{
   if (missing(groups)) return(attr(object, "Dim"))
   val <- Dim.corRStruct(object, groups)
   val[["start"]] <-
      c(0, cumsum(val[["len"]] * (val[["len"]] - 1)/2)[-val[["M"]]])
   ## will use third component of Dim list for spClass
   names(val)[3] <- "spClass"
   val[[3]] <- match(class(object)[1], c("corRExp2", "corRExpwr2"), 0)

   val
}


getCovariate.corRSpatioTemporal <- function(object, form = formula(object), data)
{
   covar <- attr(object, "covariate")

   if (is.null(covar)) { # need to calculate it
      if (missing(data)) {
         stop("Need data to calculate covariate")
      }
      covForm <- terms(getCovariateFormula(form))
      attr(covForm, "intercept") <- 0
      tcovar <- length(all.vars(covForm))
      if (tcovar >= 2) { # covariates present
         covar <- model.matrix(covForm,
                     model.frame(covForm, data, drop.unused.levels = TRUE))
      } else if (tcovar == 1) {
         covar <- model.matrix(covForm,
                     model.frame(covForm, data, drop.unused.levels = TRUE))
         covar <- cbind(covar, 1:nrow(data))
         tcovar <- 2
      } else {
         covar <- cbind(1:nrow(data), 1:nrow(data))
         tcovar <- 2
      }

      if (nrow(covar) > nrow(unique(covar))) {
         stop("Cannot have zero distances in \"corRSpatioTemporal\"")
      }

      if (is.null(getGroupsFormula(form))) { # no groups
         attr(covar, "assign") <- NULL
         attr(covar, "contrasts") <- NULL
         x <- as.vector(dist2(covar[, -tcovar], method = attr(object, "metric"),
                                                r = attr(object, "radius")))
         attr(covar, "dist") <- x
         minD <- ifelse(any(x > 0), min(x[x > 0]), 0)

         idx <- lower.tri(matrix(0, nrow(covar), nrow(covar)))
         x <- abs(covar[col(idx)[idx], tcovar] - covar[row(idx)[idx], tcovar])
         attr(covar, "period") <- x
         minD <- c(minD, ifelse(any(x > 0), min(x[x > 0]), 0))

      } else { # by groups
         grps <- getGroups(object, data = data)
         covar <- lapply(split(as.data.frame(covar), grps),
            function(el, metric, radius) {
               el <- as.matrix(el)
               attr(el, "dist") <- as.vector(dist2(el[, -tcovar], metric,
                                                   r = radius))
               idx <- lower.tri(matrix(0, nrow(el), nrow(el)))
               attr(el, "period") <- abs(el[col(idx)[idx], tcovar] -
                                         el[row(idx)[idx], tcovar])
               el
            }, metric = attr(object, "metric"), radius = attr(object, "radius"))
         x <- unlist(lapply(covar, attr, which = "dist"))
         minD <- ifelse(any(x > 0), min(x[x > 0]), 0)
         x <- unlist(lapply(covar, attr, which = "period"))
         minD <- c(minD, ifelse(any(x > 0), min(x[x > 0]), 0))
      }
      attr(covar, "minD") <- minD
   }

   covar
}


corMatrix.corRSpatioTemporal <- function(object, covariate = getCovariate(object),
   corr = TRUE, ...)
{
   if (data.class(covariate) == "list") {
      dist <- unlist(lapply(covariate, attr, which = "dist"))
      period <- unlist(lapply(covariate, attr, which = "period"))
      len <- unlist(lapply(covariate, nrow))
   } else {
      dist <- attr(covariate, "dist")
      period <- attr(covariate, "period")
      len <- nrow(covariate)
      names(len) <- 1
   }

   par <- coef(object)
   val <- switch(class(object)[1],
      corRExp2   = cor.exp2(dist, period, par[1], 1, par[2], 1, par[3]),
      corRExpwr2 = cor.exp2(dist, period, par[1], par[2], par[3], par[4], par[5])
   )

   val <- split(val, rep(names(len), len * (len - 1) / 2))
   lD <- NULL
   for(i in names(val)) {
      x <- matrix(0, len[i], len[i])
      x[lower.tri(x)] <- val[[i]]
      if (corr) {
         val[[i]] <- x + t(x)
         diag(val[[i]]) <- 1
      } else {
         diag(x) <- 1
         l <- chol(t(x))
         val[[i]] <- t(backsolve(l, diag(len[i])))
         lD <- c(lD, diag(l))
      }
   }
   if (length(len) == 1) val <- val[[1]]
   if (!is.null(lD)) attr(val, "logDet") <- -1 * sum(log(lD))

   val
}


corFactor.corRSpatioTemporal <- function(object, ...)
{
   val <- corMatrix(object, corr = FALSE, ...)
   lD <- attr(val, "logDet")
   if (is.list(val)) val <- unlist(val)
   else val <- as.vector(val)
   names(val) <- NULL
   attr(val, "logDet") <- lD

   val
}


coef.corRSpatioTemporal <- function(object, ...)
{
   val <- as.vector(object)
   if (length(val) == 0) {
      return(val)
   }
   names(val) <- rownames(attr(object, "bounds"))

   val
}


"coef<-.corRSpatioTemporal" <- function(object, ..., value)
{
   if (!all(inbounds(value, attr(object, "bounds")))) stop()
   object[] <- value

   object
}


################################################################################
# corRExp2 - Exponential spatio-temporal correlation structure
################################################################################

corRExp2 <- function(value = numeric(0), form = ~ 1,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956)
{
   attr(value, "formula") <- form
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "bounds") <- matrix(c(0, 0, 0, Inf, Inf, Inf, 1, 1, 2), ncol=3,
      dimnames = list(c("spatial range", "temporal range", "interaction"),
                      c("lower", "upper", "type")))
   class(value) <- c("corRExp2", "corRSpatioTemporal", "corRStruct")

   value
}


Initialize.corRExp2 <- function(object, data, ...)
{
   if (!is.null(attr(object, "covariate"))) { # already initialized
      return(object)
   }
   object <- Initialize.corRStruct(object, data)

   val <- as.vector(object)
   if (length(val) == 0) {
      val <- attr(getCovariate(object), "minD") * 0.9
      val[val == 0] <- 1
      val <- c(val, 0)
   } else if (!all(inbounds(val, attr(object, "bounds")))) {
      stop()
   }
   attributes(val) <- attributes(object)

   val
}


################################################################################
# corRExpwr2 - Powered exponential spatio-temporal correlation structure
################################################################################

corRExpwr2 <- function(value = numeric(0), form = ~ 1,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956)
{
   attr(value, "formula") <- form
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "bounds") <- matrix(c(0, 0, 0, 0, 0, Inf, 2, Inf, 2, Inf,
                                     1, 3, 1, 3, 2), ncol=3,
      dimnames = list(c("spatial range", "spatial shape", "temporal range",
                        "temporal shape", "interaction"),
                      c("lower", "upper", "type")))
   class(value) <- c("corRExpwr2", "corRSpatioTemporal", "corRStruct")

   value
}


Initialize.corRExpwr2 <- function(object, data, ...)
{
   if (!is.null(attr(object, "covariate"))) { # already initialized
      return(object)
   }
   object <- Initialize.corRStruct(object, data)

   val <- as.vector(object)
   if (length(val) == 0) {
      val <- attr(getCovariate(object), "minD") * 0.9
      val[val == 0] <- 1
      val <- c(val[1], 1, val[2], 1, 0)
   } else if (!all(inbounds(val, attr(object, "bounds")))) {
      stop()
   }
   attributes(val) <- attributes(object)

   val
}


################################################################################
# corRExpwrDt - Joint powered exponential spatial and integrated exponential
#               temporal correlation structure
################################################################################

corRExpwr2Dt <- function(value = numeric(0), form = ~ 1,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956)
{
   attr(value, "formula") <- form
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "bounds") <- matrix(c(0, 0, 0, 0, Inf, 2, Inf, Inf,
                                     1, 3, 1, 2), ncol=3,
      dimnames = list(c("spatial range", "spatial shape", "temporal range",
                        "interaction"), c("lower", "upper", "type")))
   class(value) <- c("corRExpwr2Dt", "corRSpatioTemporal", "corRStruct")

   value
}


Initialize.corRExpwr2Dt <- function(object, data, ...)
{
   if (!is.null(attr(object, "covariate"))) { # already initialized
      return(object)
   }
   object <- Initialize.corRStruct(object, data)

   val <- as.vector(object)
   if (length(val) == 0) {
      val <- attr(getCovariate(object), "minD") * 0.9
      val[val == 0] <- 1
      val <- c(val[1], 1, val[2], 0)
   } else if (!all(inbounds(val, attr(object, "bounds")))) {
      stop()
   }
   attributes(val) <- attributes(object)

   val
}


Dim.corRExpwr2Dt <- function(object, groups, ...)
{
   if (missing(groups)) return(attr(object, "Dim"))
   val <- Dim.corRStruct(object, groups)
   val[["start"]] <-
      c(0, cumsum(val[["len"]] * (val[["len"]] - 1)/2)[-val[["M"]]])
   ## will use third component of Dim list for spClass
   names(val)[3] <- "spClass"
   val[[3]] <- match(class(object)[1], c("corRExpwr2Dt"), 0)

   val
}


getCovariate.corRExpwr2Dt <- function(object, form = formula(object), data)
{
   covar <- attr(object, "covariate")

   if (is.null(covar)) { # need to calculate it
      if (missing(data)) {
         stop("Need data to calculate covariate")
      }
      covForm <- terms(getCovariateFormula(form))
      attr(covForm, "intercept") <- 0
      tcovar <- length(all.vars(covForm)) + c(-1, 0)
      if (tcovar[1] >= 2) { # covariates present
         covar <- model.matrix(covForm,
                     model.frame(covForm, data, drop.unused.levels = TRUE))
      } else if (tcovar[1] == 1) {
         covar <- model.matrix(covForm,
                     model.frame(covForm, data, drop.unused.levels = TRUE))
         covar <- cbind(1:nrow(data), covar)
         tcovar <- tcovar + 1
      } else {
         covar <- matrix(1:nrow(data), nrow(data), 3)
         tcovar <- c(2, 3)
      }

      if (nrow(covar) > nrow(unique(covar))) {
         stop("Cannot have duplicate sites in \"corRExpwr2Dt\"")
      } else if (any(covar[,tcovar[1]] > covar[,tcovar[2]])) {
         stop("Temporal limits must be ascending in \"corRExpwr2Dt\"")
      }

      if (is.null(getGroupsFormula(form))) { # no groups
         attr(covar, "assign") <- NULL
         attr(covar, "contrasts") <- NULL
         x <- as.vector(dist2(covar[, -tcovar], method = attr(object, "metric"),
                                                r = attr(object, "radius")))
         attr(covar, "dist") <- x
         minD <- ifelse(any(x > 0), min(x[x > 0]), 0)

         idx <- lower.tri(matrix(0, nrow(covar), nrow(covar)))
         t1 <- covar[col(idx)[idx], tcovar]
         t2 <- covar[row(idx)[idx], tcovar]
         attr(covar, "t1") <- t1
         attr(covar, "t2") <- t2
         x <- abs((t2 - t1) %*% c(0.5, 0.5))
         minD <- c(minD, ifelse(any(x > 0), min(x[x > 0]), 0))

      } else { # by group
         grps <- getGroups(object, data = data)
         covar <- lapply(split(as.data.frame(covar), grps),
            function(el, metric, radius) {
               el <- as.matrix(el)
               attr(el, "dist") <- as.vector(dist2(el[, -tcovar], metric,
                                                   r = radius))
               idx <- lower.tri(matrix(0, nrow(el), nrow(el)))
               attr(el, "t1") <- el[col(idx)[idx], tcovar]
               attr(el, "t2") <- el[row(idx)[idx], tcovar]
               el
            }, metric = attr(object, "metric"), radius = attr(object, "radius"))
         x <- unlist(lapply(covar, attr, which = "dist"))
         minD <- ifelse(any(x > 0), min(x[x > 0]), 0)
         x <- rapply(covar,
            function(x) abs((attr(x, "t2") - attr(x, "t1")) %*% c(0.5, 0.5)))
         minD <- c(minD, ifelse(any(x > 0), min(x[x > 0]), 0))
      }
      attr(covar, "minD") <- minD
   }

   covar
}


corMatrix.corRExpwr2Dt <- function(object, covariate = getCovariate(object),
   corr = TRUE, ...)
{
   if (data.class(covariate) == "list") covar <- covariate
   else covar <- list(covariate)

   par <- coef(object)

   val <- list()
   lD <- NULL
   for(i in seq(covar)) {
      r <- cor.exp2dt(attr(covar[[i]], "dist"),
                      attr(covar[[i]], "t1"), attr(covar[[i]], "t2"),
                      par[1], par[2], par[3], par[4])
      x <- matrix(0, nrow(covar[[i]]), nrow(covar[[i]]))
      x[lower.tri(x)] <- r
      idx <- ncol(covar[[i]]) + c(-1, 0)
      if (corr) {
         val[[i]] <- x + t(x)
         diag(val[[i]]) <- cor.exp2dt(0, covar[[i]][,idx], covar[[i]][,idx],
                                      par[1], par[2], par[3], par[4])
      } else {
         diag(x) <- cor.exp2dt(0, covar[[i]][,idx], covar[[i]][,idx],
                               par[1], par[2], par[3], par[4])
         l <- chol(t(x))
         val[[i]] <- t(backsolve(l, diag(length(covar[[i]]))))
         lD <- c(lD, diag(l))
      }
   }
   if (length(val) == 1) val <- val[[1]]
   if (!is.null(lD)) attr(val, "logDet") <- -1 * sum(log(lD))

   val
}


corFactor.corRExpwr2Dt <- function(object, ...)
{
   val <- corMatrix(object, corr = FALSE, ...)
   lD <- attr(val, "logDet")
   if (is.list(val)) val <- unlist(val)
   else val <- as.vector(val)
   names(val) <- NULL
   attr(val, "logDet") <- lD

   val
}


################################################################################
# Distance and correlation functions
################################################################################

dist2 <- function(x, method = c("euclidean", "maximum", "manhattan", "canberra",
   "binary", "minkowski", "haversine"), diag = FALSE, upper = FALSE,
   p = 2, r = 3956)
{
   METHOD <- match.arg(method)
   switch(METHOD,
      haversine = {
              m <- matrix(NA, nrow(x), nrow(x))
              idx <- lower.tri(m)
              m[idx] <- haversine(x[col(m)[idx],1:2], x[row(m)[idx],1:2], r)
              d <- as.dist(m, diag = diag, upper = upper)
            },
            {
              f <- get("dist", envir = as.environment("package:stats"))
              d <- f(x, method = METHOD, diag = diag, upper = upper, p = p)
            }
   )

   d
}

## Anisotropic transformation of coordinates
anisotropic <- function(x, par, system = c("cartesian", "polar", "spherical"))
{
   X <- as.matrix(x)
   system <- match.arg(system)

   d <- ncol(x)
   if (d == 2) {
      if (length(par) != 2)
         stop("parameter vector must consist of two elements",
              " - an anisotropy angle and ratio")
      if ((r <- par[2]) < 1) stop("anisotropy ratios must be >= 1")
      S <- diag(c(1, 1 / r))
      alpha <- par[1]
      R <- matrix(c(cos(alpha), -sin(alpha),
                    sin(alpha),  cos(alpha)), 2, 2)
      switch(system,
         cartesian = {
            Y <- X %*% S %*% R
         },
         polar = {
            theta <- X[,2]
            Z <- X[,1] * cbind(cos(theta), sin(theta)) %*% S %*% R
            x <- Z[,1]
            y <- Z[,2]
            theta <- if (x >= 0 && y >= 0) atan(abs(y / x))
                     else if (x < 0 && y >= 0) pi - atan(abs(y / x))
                     else if (x < 0 && y < 0) pi + atan(abs(y / x))
                     else 2 * pi - atan(abs(y / x))
            Y <- cbind(sqrt(x^2 + y^2), theta)
         },
         spherical = {
         }
      )
   } else if (d == 3) {
      if (length(par) != 5)
         stop("parameter vector must consist of five elements",
              " - three anisotropy angles and two ratios")
      if (any((r <- par[4:5]) < 1)) stop("anisotropy ratios must be >= 1")
      S <- diag(c(1, 1 / r))
      alpha <- par[1]
      beta <- par[2]
      theta <- par[3]
      R1 <- matrix(c(cos(alpha), -sin(alpha), 0,
                     sin(alpha),  cos(alpha), 0,
                               0,          0, 1), 3, 3)
      R2 <- matrix(c( cos(beta), 0, sin(beta),
                              0, 1,          0,
                     -sin(beta), 0, cos(beta)), 3, 3)
      R3 <- matrix(c(1,          0,           0,
                     0, cos(theta), -sin(theta),
                     0, sin(theta),  cos(theta)), 3, 3)
      R <- R1 %*% R2 %*% R3
      switch(system,
         cartesian = {
            Y <- X %*% S %*% R
         },
         polar = {
            stop("polar coordinates must be of two dimensions")
         },
         spherical = {
            phi <- X[,2]
            theta <- X[,3]
            Z <- X[,1] * cbind(sin(phi) * cos(theta), sin(phi) * sin(theta),
                               cos(theta)) %*% S %*% R
            x <- Z[,1]
            y <- Z[,2]
            z <- Z[,3]
            Y <- cbind(sqrt(x^2 + y^2 + z^2), atan(y / x),
                       atan(sqrt(x^2 + y^2) / z))
         }
      )
   } else {
      stop("anisotropy supported only for 2-D and 3-D coordinate systems")
   }

   Y
}

## Great circle distance
## Coordinates specified as cbind(longitude, latidute)
haversine <- function(x, y, r = 3956)
{
   if(is.vector(x)) x <- matrix(x, 1, 2)
   if(is.vector(y)) y <- matrix(y, 1, 2)

   rad <- pi / 180
   z <- sin((y - x) * (rad / 2))^2
   a <- z[,2] + cos(rad * x[,2]) * cos(rad * y[,2]) * z[,1]
   (2 * r) * atan2(sqrt(a), sqrt(1 - a))
}

## Powered exponential correlation function
cor.exp <- function(x, range = 1, p = 1)
{
   if (range <= 0 || p <= 0)
      stop("Exponential correlation parameter must be > 0")

   if (p == 1) exp(x / (-1 * range))
   else exp(-1 * (x / range)^p)
}

## Gneiting correlation function
cor.gneiting <- function(x, range = 1)
{
   if (range <= 0)
      stop("Gneiting correlation parameter must be > 0")

   range <- range / 0.3008965026325734   # sqrt(2) * 10 / 47
   r <- (x < range)
   x0 <- x[r] / range
   r[r] <- (1 + 8 * x0 + 25 * x0^2 + 32 * x0^3) * (1 - x0)^8
   r
}

## Linear correlation function
cor.lin <- function(x, range = 1)
{
   if (range <= 0)
      stop("Linear correlation parameter must be > 0")

   r <- (x < range)
   r[r] <- 1 - x[r] / range
   r
}

## Matern correlation function
cor.matern <- function(x, range = 1, scale = 1)
{
   if(range <= 0 || scale <= 0)
      stop("Matern correlation parameters must be > 0")

   idx <- (x > 0)
   r <- as.double(!idx)
   x0 <- x[idx] / range
   r[idx] <- x0^scale * besselK(x0, scale) / (2^(scale - 1) * gamma(scale))
   r
}

## Cauchy correlation function
cor.cauchy <- function(x, range = 1)
{
   if (range <= 0)
      stop("Cauchy correlation parameter must be > 0")

   1 / (1 + (x / range)^2)
}

## Sperical correlation function
cor.spher <- function(x, range = 1)
{
   if (range <= 0)
      stop("Spherical correlation parameter must be > 0")

   r <- (x < range)
   x0 <- x[r] / range
   r[r] <- 1 - 1.5 * x0 + 0.5 * x0^3
   r
}

## Sine wave correlation function
cor.wave <- function(x, range = 1)
{
   if (range <= 0)
      stop("Sine wave correlation parameter must be > 0")

   x0 <- (x / range)
   sin(x0) / x0
}

## Non-separable exponential spatio-temporal correlation function
cor.exp2 <- function(x, t, x.range = 1, x.p = 1, t.range = 1, t.p = 1, lambda = 0)
{
   if (t.range <= 0 || x.range <= 0 || x.p <= 0 || lambda < 0)
      stop("Exponential correlation parameters must be > 0")

   x0 <- if (x.p == 1) x / (-1 * x.range)
         else -1 * (x / x.range)^x.p
   t0 <- if (t.p == 1) t / (-1 * t.range)
         else -1 * (t / t.range)^t.p

   exp(x0 - lambda * x0 * t0 + t0) 
}

## Non-separable temporally integrated exponential spatial correlation function
cor.exp2dt <- function(x, t1, t2, x.range = 1, x.p = 1, t.range = 1, lambda = 0)
{
   if (t.range <= 0 || x.range <= 0 || x.p <= 0 || lambda < 0)
      stop("Exponential correlation parameters must be > 0")

   if (is.vector(t1)) t1 <- matrix(t1, 1, 2)
   if (is.vector(t2)) t2 <- matrix(t2, 1, 2)

   x0 <- if (x.p == 1) x / (-1 * x.range)
         else -1 * (x / x.range)^x.p

   overlap <- pmin(t1[,2], t2[,2]) - pmax(t1[,1], t2[,1])
   overlap[overlap < 0] <- 0
   norm <- (t1 %*% c(-1, 1)) * (t2 %*% c(-1, 1))

   if (lambda == 0) theta <- t.range
   else theta <- t.range / (1 - lambda * x0)

   val <- (theta^2 * exp(abs(t1[, c(1,1,2,2)] - t2[, c(2,1,2,1)]) /
             (-1 * theta)) %*% c(1, -1, -1, 1) + 2 * theta * overlap) / norm
   exp(x0) * as.vector(val)
}
