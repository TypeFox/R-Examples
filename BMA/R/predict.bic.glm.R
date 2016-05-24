predict.bic.glm <-
function(object,newdata,...){
# CF October 2011

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

get.names <- function(x,maxCol) {
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
	warning(paste("There were ", length(omitted), "records deleted due to N\
A'\
s"))
    }
#    if (drop.factor.levels) {
#       cdf <- cbind.data.frame(y = y, x)
#        mm <- model.matrix(formula(cdf), data = cdf)[, -1, drop = FALSE]
#        x <- mm
#   }
    xx <- dropcols(x, y, wt, maxCol)
    xnames <- xx$var.names[-1]
    x <- xx$mm
    reduced <- xx$any.dropped
    dropped <- NULL
    if (reduced)
        dropped <- xx$dropped

    xnames
  }

linkinvBinom <- function(x) {
## written by CF

  dimx <- dim(x)
  x <- as.vector(x)  

  fine <- is.finite(x)

  pos <- x >= 0

  small <- -x < log(.Machine$double.eps)

  index <- fine & pos & small
  x[index] <- 1

  index <- fine & pos & !small
  xindex   <- exp(-x[index])
  x[index] <- 1/(1+xindex)

  small <- x < log(.Machine$double.xmin)

  index <- fine & !pos & small
  x[index] <- 0

  index <- fine & !pos & !small
  xindex   <- exp(x[index])
  x[index] <- xindex/(1+xindex)

  if (!is.null(dimx)) array(x,dimx) else x
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
    callList <- as.list(object$call)
    callFunc <- callList[[1]]
    defaults <- as.list(args(deparse(substitute(callFunc))))

    factor.type <- callList$factor.type
    factor.type <- if (is.null(factor.type)) defaults$factor.type else eval(factor.type)
    if (as.character(callFunc) ==  "bic.glm.matrix") factor.type <- FALSE      

    maxCol <- callList$maxCol
    maxCol <- if (is.null(maxCol)) defaults$maxCol else eval(maxCol)

    newdata <- data.frame(newdata)

    if (!is.null(object$formula)) {
       newdata <- model.matrix(object$formula, data = newdata)[,-1]
    }
    else {
       y <- rnorm(nrow(newdata))
       newdata <- model.matrix(formula(cbind.data.frame(y = y,
                         newdata)), data = newdata)[,-1]
    } 

    nam <- colnames(object$mle)[-1]

    mvars <- match( nam, colnames(newdata), nomatch = 0)

    if (any(mvars == 0)) stop("newdata is missing variables")

    newdata <- newdata[,nam, drop = FALSE]

    nObs <- nrow(newdata)

    newnam <- get.names( newdata, maxCol)

    fn <- factor.names(newdata)

    new.assign <- create.assign(newdata)

    fac.levels <- unlist(lapply(new.assign, length)[-1])
    fac.ident <- sapply(newdata,is.factor)

    nModels <- length(object$postprob)
    linpred <- matrix(NA,nModels,nObs)

    dimnames(object$which) <- list(dimnames(object$which)[[1]],
                                     object$xnames)

    colnames(object$mle) <- c("(Intercept)",newnam)

    if (any(fac.ident)) {
      nVars <- length(fac.ident) 
      for (k in 1:nModels) {
         linpred[k,] <- object$mle[k,1]
         IN <- object$mle[k,] != 0
         j <- 1
         for (i in 1:nVars) {
            if (fac.ident[i]) {
               L <- (1:as.numeric(fac.levels[i]))
               if (any(IN[L])) {
# works for factor.type = FALSE if all levels are in some model
                 coefs <- object$mle[k,j+L]
                 if (!all(coefs == 0))  {
                   fac.vals <- newdata[, i]
                   m <- match(fac.vals,fn[[i]],nomatch=NA)
                   if (any(is.na(m))) stop("NA")
                   linpred[k,] <- linpred[k,] + as.vector(c(0,coefs)[m])
                 }
               }
               j <- j + length(L)
            }
          else {
            j <- j + 1
            linpred[k,] <- linpred[k,]+ object$mle[k,j]*newdata[,i]
          }
        }
      }
     }
   else {
      linpred <- tcrossprod(object$mle, cbind(1, as.matrix(newdata)))
   }


  rhs <- as.function(object$linkinv)(linpred)

# rhs =  predictions for individual models
# if (show) print(rhs[,1:min(9,ncol(rhs))])

  pred <- apply(object$postprob * rhs, 2, sum)
  names(pred) <- rownames(newdata)

  pred
}

