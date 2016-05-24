  predict.sitar <- function(object, newdata=getData(object), level=1, ...,
                            deriv=0, abc=ranef(object),
                            xfun=function(x) x, yfun=xfun) {
# create x and id variables in newdata
    oc <- object$call.sitar
    if (is.null(newdata$x)) newdata$x <- eval(oc$x, newdata)
    x <- newdata$x
    if (is.null(xoffset <- object$xoffset)) {
      xoffset <- mean(getCovariate(object))
      warning('xoffset set to mean(x) - best to refit model')
    }
    newdata$x <- newdata$x - xoffset
# create id in newdata
    if (any(level == 1)) newdata$id <- eval(oc$id, newdata)
    else newdata$id <- rep.int(getGroups(object)[1], nrow(newdata))
    id <- newdata$id <- factor(newdata$id)
# check abc
    if (abcset <- !is.data.frame(abc)) {
      abc <- data.frame(t(abc))
      abc[, letters[1:3][!letters[1:3] %in% names(abc)]] <- 0 # fill with zeros
      level <- 0
      id <- rep.int(1, length(x))
    }
    abc[, letters[1:3][!letters[1:3] %in% names(ranef(object))]] <- 0 # zeros if not in model
    abc <- abc[id, ]
# check if old-style object lacking fitnlme
    if(!'fitnlme' %in% names(object)) {
      warning('fitnlme missing - best to refit model')
      object <- update(object, control=nlmeControl(maxIter=0, pnlsMaxIter=0, msMaxIter=0))
    }
# attach object for fitnlme
    on.exit(detach(object))
    eval(parse(text='attach(object)'))
# identify covariates needed in newdata, omitting fixed effects and x
    argnames <- names(formals(fitnlme))
    argnames <- argnames[!argnames %in% names(fixef(object))][-1]
    if (length(argnames) > 0) {
# check if newdata subsetted (from plot)
      if (is.null(subset <- attr(newdata, 'subset'))) {
# identify covariates in newdata other than x and id
    		covnames <- names(newdata)
    		covnames <- covnames[!covnames %in% c('x', 'id')]
# set to 0 covariates not in newdata
    		newdata[, argnames[!argnames %in% covnames]] <- 0
# centre each needed covariate in newdata
    		covnames <- covnames[covnames %in% argnames]
    		if (length(covnames) > 0) {
    		  gd <- getData(object)
    		  for (i in covnames) {
# continuous variable
    		    if (i %in% argnames) newdata[[i]] <- newdata[[i]] - mean(gd[[i]])
    	      else {
# factor as instrumental variable(s)
    	        lev <- levels(gd[[i]])
    	        for (j in 2:length(lev)) {
    	          k <- paste0(i, lev[j])
    	          newdata[[k]] <- as.numeric(newdata[[i]] == lev[j]) - mean(gd[[i]] == lev[j])
    	        }
    	      }
    		  }
    		}
      }
# newdata subsetted (in plot)
      else {
        gd <- update(object, returndata=TRUE)[subset, argnames]
        argnames <- unlist(lapply(gd, mean))
        newdata <- data.frame(newdata, t(argnames))
      }
    }
# set class to nlme
    class(object) <- class(object)[-1]
# simple prediction
    if (deriv == 0 && !abcset) {
      pred <- yfun(predict(object=object, newdata=newdata, level=level, ...))
      return(pred)
    }
# complex prediction
    else { # deriv == 1 || abcset
# mean distance curve
      pred <- predict(object=object, newdata=newdata, level=0, ...)
# DISTANCE
      if (deriv == 0) { # abcset
# level 1 prediction based on x changed to reflect individual b and c
        pred <- spline(list(x=x, y=pred), method='natural',
                      xout=xyadj(x=x, id=id, object=object, abc=abc)$x)$y
# add individual a to prediction
        if (!is.null(abc$a)) pred <- yfun(pred + abc$a)
      }
# VELOCITY
      else { # deriv == 1
# mean velocity curve on back-transformed axes
        vel0 <- predict(makess(x, pred, xfun=xfun, yfun=yfun), xfun(x), deriv=1)
        if (any(level == 0) && !abcset) pred0 <- pred <- vel0$y
        if (any(level == 1) || abcset) {
# level 1 prediction based on x changed to reflect individual b and c
          pred <- spline(vel0, method='natural',
                         xout=xfun(xyadj(x=x, id=id, object=object, abc=abc)$x))$y
# multiply velocity by individual c (inexact when xfun != I)
          if (!is.null(abc$c)) pred <- pred * exp(abc$c)
        }
      }
# return data frame if level 0:1
      if (length(level) > 1) return(data.frame(id=id, predict.fixed=pred0, predict.id=pred))
# add names or split by id if level 1
      if (level == 1) {
        asList <- ifelse(is.null(asList <- list(...)$asList), FALSE, asList)
        if (asList) pred <- split(pred, id) else names(pred) <- id
      }
      attr(pred, 'label') <- 'Predicted values'
      return(pred)
    }
  }

  getData.sitar <- function(object) {
    object$call <- object$call.sitar
    class(object) <- 'lme'
    getData(object)
  }

  getVarCov.sitar <- function(obj, ...) {
    class(obj) <- 'lme'
    getVarCov(obj)
  }

  getCovariate.sitar <- function(object, ...)
  {
    eval(object$call.sitar$x, getData(object))
  }
