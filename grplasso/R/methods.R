print.grplasso <- function(x, digits = 4, ...)
{
  ## Purpose: Print an object of class "grplasso"
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## x:      Object of class "grplasso"
  ## digits: nr of digits to print
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 30 Mar 2006, 18:01

  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")

  cat("* Nr. of (penalized) predictor groups:",
      length(unique(na.omit(x$index))), "\n")
  cat("* Nr. of predictors (dummy variables):", length(x$index), "\n")
  cat("  whereof", sum(is.na(x$index)), "are not penalized", "\n")
  cat("* Nr. of observations:", length(x$y), "\n")
  cat("* Penalty parameter lambda:\n")
  cat("    Number of grid points: ", length(x$lambda), "\n")
  cat("    Min value: ", format(min(x$lambda), digits = digits), "\n")
  cat("    Max value: ", format(max(x$lambda), digits = digits))

  cat("\n\n")
  
  invisible(x)
}

plot.grplasso <- function(x, type = "coefficients", col = NULL, ...)
{
  ## Purpose: Plots the solution path of a "grplasso" object. The x-axis
  ##          is the penalty parameter lambda, the y-axis can be
  ##          coefficients or the l2-norm of the coefficient groups.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## x:    a grplasso object
  ## type: what should be on the y-axis? Coefficients (dummy
  ##       variables)?
  ## col:  a vector indicating the color of the different solution
  ##       paths. The length should equal the number of coefficients. 
  ## ...:  other parameters to be passed to the plotting functions.
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  7 Apr 2006, 08:34

  if(length(x$lambda) == 1)
    stop("Plot function not available for a one dimensional lambda")
  
  type <- match.arg(type)
  
  if(is.null(col))
    col <- 1:length(unique(x$index))

  dict.pen <- unique(na.omit(x$index))
  nr.npen <- sum(is.na(x$index))
  dict.pen.ord <- ((nr.npen) + 1 : length(dict.pen))

  if(type == "coefficients"){
    coef.col <- numeric()
    index.ord <- numeric()
    index.ord[is.na(x$index)] <- 1:nr.npen

    for(j in 1:length(dict.pen))
      index.ord[x$index == sort(dict.pen)[j]] <- dict.pen.ord[j]
        
    coef.col <- col[index.ord]

    matplot(x$lambda, t(coef(x)), type = "l",
            xlab = "Lambda", ylab = "Coefficients", col = coef.col,
            main = "Coefficient paths", xlim = c(max(x$lambda), min(x$lambda)),
            ...)
    axis(4, at = coef(x)[, ncol(coef(x))], labels = rownames(coef(x)))
  }
}

predict.grplasso <- function(object, newdata, 
                             type = c("link", "response"),
                             na.action = na.pass, ...)
{
  ## Purpose: Obtains predictions from a "grplasso" object.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## object:  a "grplasso" object
  ## newdata: data.frame or design matrix of observations at which
  ##          predictions are to be made.
  ## type: the type of prediction. type = "link" is on the
  ##       scale of linear predictors, whereas type = "response" is on
  ##       the scale of the response variable, i.e. type = "response"
  ##       applies the inverse link function on the linear predictors.
  ## ...:  other options to be passed to the predict function.
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  7 Apr 2006, 09:04

  type <- match.arg(type)
  
  na.act <- object$na.action

  ## If no new data is available, use the information in the fit object
  if(missing(newdata) || is.null(newdata)){
    pred <- switch(type,
                   link = object$linear.predictors,
                   response = fitted(object))
    if(!is.null(na.act))
      pred <- napredict(na.act, pred)

    if(dim(pred)[2] == 1)
      pred <- pred[,1,drop = TRUE]

    attr(pred, "lambda") <- object$lambda

    return(pred)
  }
  if(!is.null(tt <- object$terms)){ ## if we have a terms object in the fit
    newdata <- as.data.frame(newdata)
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action,
                     xlev = object$xlevels)
    offset <- attr(tt, "offset")

    if(!is.null((cl <- attr(Terms, "dataClasses"))))
      .checkMFClasses(cl, m)
    x <- model.matrix(Terms, m, contrasts = object$contrasts)
    pred <- x %*% coef(object)

    ## new code for offset handling
    offset <- rep(0, nrow(x))
    if (!is.null(off.num <- attr(tt, "offset"))) 
      for (i in off.num) offset <- offset +
        eval(attr(tt, "variables")[[i + 1]], newdata)
    if (!is.null(object$call$offset)) 
      offset <- offset + eval(object$call$offset, newdata)

    pred <- pred + offset
    ## end new code for offset handling

    ## old code
##-     if(!is.null(offset)){
##-       offset <- eval(attr(tt, "variables")[[offset]], newdata)
##-       pred <- pred + offset
##-     }
    ## end old code
  }else{ ## if the object comes from grplasso.default
    x <- as.matrix(newdata)
    pred <- x %*% coef(object)
    if(any(object$offset != 0))
      warning("Possible offset not considered!")
  }
    
  pred <- switch(type,
                 link = pred,
                 response = object$model@invlink(pred))
                 ##apply(pred, 2, object$model@invlink))

  if(!is.null(na.action))
    pred <- napredict(na.action, pred)

  if(dim(pred)[2] == 1)
    pred <- pred[,1,drop = TRUE]

  attr(pred, "lambda") <- object$lambda
  pred
}

fitted.grplasso <- function(object, ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 26 Jun 2006, 12:11

  out <- object$fitted
  attr(out, "lambda") <- object$lambda
  out
}

"[.grplasso" <- function(x, i){
  
  ## First get dimensions of the original object x

  nrlambda <- length(x$lambda)

  if(missing(i))
    i <- 1:nrlambda

  ## Error checking
  ## ...

  ## Subset the object
  fit.red <- x
  
  fit.red$coefficients <- coef(x)[,i,drop = FALSE]
  if(length(fit.red$coefficients) == 0)
    stop("Not allowed to remove everything!")

  fit.red$lambda       <- x$lambda[i]
  fit.red$ngradient    <- x$ngradient[,i,drop = FALSE]
  fit.red$nloglik      <- x$nloglik[i]
  fit.red$fitted       <- fitted(x)[,i]
  fit.red$linear.predictors <- x$linear.predictors[,i,drop = FALSE]
  fit.red$fn.val       <- x$fn.val[i]
  
  fit.red
}


