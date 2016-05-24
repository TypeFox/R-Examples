predict.bicreg <-
function( object, newdata, quantiles = c(.1,.5,.9), ...){
# CF August 2011 - January 2012

  cdfBMAnormal <-
function (x, WEIGHTS, MEAN, SD, offset = 0)
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
  sum(WEIGHTS*pnorm(x, mean = MEAN, sd = SD)) - offset
}

quantBMAnormal <-
function(alpha, WEIGHTS, MEAN, SD)
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
  lower <- min(MEAN-6*SD)
  upper <- max(MEAN+6*SD)

  if (cdfBMAnormal(lower, WEIGHTS, MEAN, SD, 0) > alpha) return(NA)
  if (cdfBMAnormal(upper, WEIGHTS, MEAN, SD, 0) < alpha) return(NA)

  uniroot(cdfBMAnormal, lower = lower, upper = upper,
          WEIGHTS=WEIGHTS, MEAN=MEAN, SD=SD, offset = alpha)$root
}

   dropcols <- function(x, y, wt, maxCols = 31) {
# CF: copied from bicreg (undocumented)
        x1.ldf <- data.frame(x, y = y)
        temp.wt <- wt
        lm.out <- lm(y ~ ., data = x1.ldf, weights = temp.wt)
        form.vars <- all.vars(formula(lm.out))[-1]
        any.dropped <- FALSE
        dropped.which <- NULL
        while (length(lm.out$coefficients) > maxCol) {
            any.dropped <- TRUE
            droplm <- drop1(lm.out, test = "none")
            dropped <- row.names(droplm)[which.min(droplm$RSS[-1]) +
                1]
            dropped.index <- match(dropped, form.vars)
            form.vars <- form.vars[-dropped.index]
            formla <- formula(paste("y", "~", paste(form.vars,
                collapse = " + "), sep = " "))
            lm.out <- lm(formla, data = x1.ldf, weights = temp.wt)
            dropped.which <- c(dropped.which, dropped)
        }
        new.var.names <- names(lm.out$coefficients)
        return(list(mm = model.matrix(lm.out)[, -1, drop = FALSE],
            any.dropped = any.dropped, dropped = dropped.which,
            var.names = new.var.names))
    }
  

  get.names <- function(x,maxCol,drop.factor.levels) {
# CF: written from portion of bicreg code
    x <- data.frame(x)

# dummy y, wt
    y <- rep(0,nrow(x))
    wt <- rep(1,nrow(x))

    if (is.null(dimnames(x)))
        dimnames(x) <- list(NULL, paste("X", 1:ncol(x), sep = ""))
#   y <- as.numeric(y)
    options(contrasts = c("contr.treatment", "contr.treatment"))
    xnames <- dimnames(x)[[2]]
    x2 <- na.omit(data.frame(x))
    used <- match(row.names(data.frame(x)), row.names(x2))
    omitted <- seq(nrow(x))[is.na(used)]
    if (length(omitted) > 0) {
      stop("NAs in newdata")
#       wt <- wt[-omitted]
        x <- x2
#       y <- y[-omitted]
        warning(paste("There were ", length(omitted), "records deleted due to NA'\
s"))
    }
    if (drop.factor.levels) {
        cdf <- cbind.data.frame(y = y, x)
        mm <- model.matrix(formula(cdf), data = cdf)[, -1, drop = FALSE]
        x <- mm
    }
    xx <- dropcols(x, y, wt, maxCol)
    xnames <- xx$var.names[-1]
    x <- xx$mm
    reduced <- xx$any.dropped
    dropped <- NULL
    if (reduced)
        dropped <- xx$dropped

    xnames
  }

  factor.names <- function(x) {
# CF: copied from bic.glm.data.frame (undocumented)
        out <- list()
        for (i in 1:ncol(x)) if (is.factor(x[, i])) 
            out[[i]] <- levels(x[, i])
        else out <- c(out, list(NULL))
        attributes(out)$names <- names(x)
        return(out)
    }

  
  create.assign <- function(x) {
# CF: copied from bic.glm.data.frame (undocumented)
        asgn <- list()
        asgn[[1]] <- 1
        cnt <- 2
        for (i in 1:ncol(x)) {
            if (!is.factor(x[, i])) 
                size <- 1
            else size <- length(levels(x[, i])) - 1
            asgn[[i + 1]] <- cnt:(cnt + size - 1)
            cnt <- cnt + size
        }
        names(asgn) <- c("(Intercept)", attributes(x)$names)
        return(asgn)
    }

    newdata <- as.data.frame(newdata[,object$input.names])
    nObs <- nrow(newdata)
    callList <- as.list(object$call)
    callFunc <- callList[[1]]
    defaults <- as.list(args(deparse(substitute(callFunc))))

    drop.factor.levels <- callList$drop.factor.levels
    if (is.null(drop.factor.levels)) {
      if (as.character(callFunc) ==  "bic.glm.matrix") {
        drop.factor.levels <- FALSE      
      }
      else {
        drop.factor.levels <- defaults$drop.factor.levels
      }
    }
    else drop.factor.levels <- eval(drop.factor.levels)

    maxCol <- callList$maxCol
    maxCol <- if (is.null(maxCol)) defaults$maxCol else eval(maxCol)

    newnam <- get.names( newdata, maxCol, drop.factor.levels)
# object$namesx should be a subset of newnam
# caution: newdata may not have the same column order as the original
    mvars <- match( object$namesx, newnam, nomatch = 0)

    if (any(mvars == 0)) stop("newdata is missing variables")

    if (!all(1:length(mvars) == sort(mvars))) 
      stop("newdata has extra variables")

    fn <- factor.names(newdata)

    new.assign <- create.assign(newdata)
    fac.levels <- unlist(lapply(new.assign, length)[-1])
    fac.ident <- sapply(newdata,is.factor)
 
    nModels <- length(object$postprob)

    nVars <- length(fac.levels) # some may not be used

    nModels <- length(object$postprob)
    linpred <- matrix(NA,nModels,nObs)

    dimnames(object$which) <- list(rownames(object$which),
                                     object$xnames)

    dimnames(object$mle) <- list(rownames(object$mle),
                                 c("(Intercept)",object$namesx))


      nVars <- length(fac.ident) 
      for (k in 1:nModels) {
         linpred[k,] <- object$mle[k,1]
         IN <- object$mle[k,] != 0
         j <- 1
         for (i in 1:nVars) {
            if (fac.ident[i]) {
               if (IN[j+1]) {
                 coefs <- object$mle[k,j+(1:as.numeric(fac.levels[i]))]
                 if (!all(coefs == 0)) {
                   fac.vals <- newdata[, i]
                   m <- match(fac.vals,fn[[i]],nomatch=NA)
                   if (any(is.na(m))) stop("NA")
                   linpred[k,] <- linpred[k,] + as.vector(c(0,coefs)[m])
                 }
               }
               j <- j + as.numeric(fac.levels[i])
            }
          else {
            j <- j + 1
            linpred[k,] <- linpred[k,]+ object$mle[k,j]*newdata[,i]
          }
        }
      }

# print(linpred) # predictions for the individual models

  predmean <- apply(object$postprob * linpred, 2, sum)

  objVARterm <- sum(object$postprob * object$residvar)
 
  predSD <- sqrt(objVARterm +
                 apply(object$postprob * (predmean - linpred)^2, 2, sum))

  predInt <- matrix( NA, nrow(newdata), length(quantiles))
  rownames(predInt) <- names(predmean) <- names(predSD) <- rownames(newdata)
  colnames(predInt) <- quantiles

  for (i in 1:nrow(newdata)) {
    predInt[i,] <- sapply( quantiles, quantBMAnormal, WEIGHTS=object$postprob,
                        MEAN=linpred[,i], SD=predSD[i]) 
  }

  list(mean  = predmean, sd = predSD, quantiles = predInt)
}

