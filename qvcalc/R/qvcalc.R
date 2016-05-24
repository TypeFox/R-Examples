qvcalc <- function(object, factorname = NULL, coef.indices = NULL,
                   labels = NULL, dispersion = NULL,
                   estimates = NULL, modelcall = NULL)
{
  coef.indices.saved <- coef.indices
  if (!is.matrix(object)) {
      model <- object
      ## special case of an unstructured Bradley-Terry model
      ## more standard lm, glm, etc. objects
      if (is.null(factorname) && is.null(coef.indices)) {
          stop("arguments \"factorname\" and \"coef.indices\" are both NULL")
      }
      if (is.null(coef.indices)) {   ## try to use factorname
          term.index <- which(attr(terms(model),"term.labels") == factorname)
          modelmat <- model.matrix(model)
          has.coef <- colnames(modelmat) %in% names(coef(model))
          coef.indices <- which(attr(modelmat,"assign")[has.coef] == term.index)
          if (length(model$xlevels[[factorname]]) == length(coef.indices)){
      ## factor has no constraint applied, eg if no intercept in model
              contmat <- diag(length(coef.indices))}
          else {
              contmat <- eval(call(model$contrasts[[factorname]],
                                   model$xlevels[[factorname]]))}
#          rownames(contmat) <- model$xlevels[[factorname]]  ## not needed?
          if (is.null(estimates))
              estimates <- contmat %*% coef(model)[coef.indices]
          covmat <- vcov(model, dispersion = dispersion)
          covmat <- covmat[coef.indices, coef.indices, drop = FALSE]
          covmat <- contmat %*% covmat %*% t(contmat)
      }
      else {
          k <- length(coef.indices)
          refPos <- numeric(0)
          if (0 %in% coef.indices) { ## there's a reference level to include
              refPos <- which(coef.indices == 0)
              coef.indices <- coef.indices[-refPos]
          }
          covmat <- vcov(model, dispersion = dispersion)
          covmat <- covmat[coef.indices, coef.indices, drop = FALSE]
          if (is.null(estimates)) estimates <- coef(model)[coef.indices]
          if (length(refPos) == 1) {
              if (length(estimates) != k) estimates <- c(0, estimates)
              covmat <- cbind(0, rbind(0, covmat))
              names(estimates)[1] <- rownames(covmat)[1] <-
                  colnames(covmat)[1] <- "(reference)"
              if (refPos != 1) {
                  if (refPos == k){
                      perm <- c(2:k, 1)
                  } else {
                      perm <- c(2:refPos, 1, (refPos + 1):k)
                  }
                  estimates <- estimates[perm]
                  covmat <- covmat[perm, perm, drop = FALSE]
              }
          }
      }
      return(qvcalc(covmat,
                    factorname = factorname,
                    coef.indices = coef.indices.saved,
                    labels = labels,
                    dispersion = dispersion,
                    estimates = estimates,
                    modelcall = model$call)
             )
      }
  else {  ##  the basic QV calculation, on a covariance matrix
      if (inherits(object, "BTabilities")) {  ## catch this special case first
          vc <- vcov(object)
          cf <- coef(object)
          if (is.null(factorname)) factorname <- attr(object, "factorname")
          if (is.null(modelcall)) modelcall <- attr(object, "modelcall")
          return(qvcalc(vc,
                        factorname = factorname,
                        labels = labels,
                        dispersion = dispersion,
                        estimates = cf,
                        modelcall = modelcall))
      }
      covmat <- object
      if (!is.null(labels)) rownames(covmat) <- colnames(covmat) <- labels
      n <- dim(covmat)[1]
      if (n <= 2) stop(
           "qvcalc works only for factors with 3 or more levels")
      simple.contrasts <- function(n, levelnames = 1:n){
          result <- list()
          for (i in 1:(n-1)){
	          for (j in (i+1):n){
		          result[[paste(levelnames[i],
                                        levelnames[j],
                                        sep = ",")]] <- c(i, j)}}
          result
      }
      qvdesign <- function(n){
          nrows <- choose(n, 2)
          m <- matrix(0, nrows, n)
          indices <- simple.contrasts(n)
          for (i in 1:nrows){
	          m[i, indices[[i]][1]] <- 1
	          m[i, indices[[i]][2]] <- 1}
          m
      }
      level <- qvdesign(n)
      contrast.variance <- function(contrast, covmat){
          if (!(is.matrix(covmat) &&
              (dim(covmat)[1] == dim(covmat)[2])))
          stop("covmat must be a square matrix")
          n <- dim(covmat)[1]
         if (length(contrast) == n && sum(contrast) == 0)
             ## arbitrary contrast vector
             return(as.vector(contrast %*% covmat %*% contrast))
         if (length(contrast) == 2 && all(contrast %in% 1:n)){
             ## simple contrast specified as an index pair
	         i <- contrast[1]
	         j <- contrast[2]
	         return(covmat[i,i] + covmat[j,j] - 2*covmat[i,j])}
         else stop("invalid contrast")
      }
      simple.contrast.variances <- function(n, covmat){
          if (!is.null(rownames(covmat)))
              levelnames <- rownames(covmat)
          else levelnames <- 1:n
          sapply(simple.contrasts(n, levelnames),
	          function(contrast){contrast.variance(contrast, covmat)})
      }
      response <- simple.contrast.variances(n,covmat)
      if (any(response <= 0)) {
          stop("not all contrasts have positive variance")
      } else response <- log(response)
      expLinear <- structure(list(
        	           family = "expLinear",
        	           link = "exp",
	    	           linkfun = function(mu) exp(mu),
        	           linkinv = function(eta) log(eta),
       		           variance = function(mu) rep(1, length(mu)),
        	           dev.resids = function(y, mu, wt) wt *
                           ((y - mu)^2),
        	           aic = function(y, n, mu, wt, dev) sum(wt) *
                           (log(dev/sum(wt) * 2 * pi) + 1) + 2,
        	           mu.eta = function (eta) 1/eta,
                           initialize = expression({
                               n <- rep(1, nobs)
                               mustart <- y}),
                       validmu = function(mu) TRUE),
                   class = "family")
      model <- glm(response ~ 0 + level, family = expLinear)
      qv <- coef(model)
      NAs <- rep(NA, length(qv))
      if (!is.null(rownames(covmat))) names(qv) <- rownames(covmat)
      frame <- data.frame(estimate = NAs,
                          SE = sqrt(diag(covmat)),
                          quasiSE = sqrt(qv),
                          quasiVar = qv,
                          row.names = names(qv))
      if (!is.null(estimates)) frame$estimate <- estimates
      relerrs <-  sqrt(exp(- residuals(model))) - 1
      ##  The above formula was corrected in v0.8-9; it
      ##  previously said 1 - sqrt(exp(residuals(model)), which is
      ##  not what should be expected for "relative error" here.
      ##  This corrected version agrees with the Biometrika paper.
      ##  Thanks to Shaun Killingbeck for spotting this error in the
      ##  previous version.
      names(relerrs) <- names(response)
      return(structure(list(covmat = covmat,
                            qvframe = frame,
 	      	            dispersion = dispersion,
                            relerrs = relerrs,
                            factorname = factorname,
                            coef.indices = coef.indices,
                            modelcall = modelcall),
 	            class="qv"))}
  }

worstErrors <- function(qv.object)
{
    reducedForm <- function(covmat, qvmat){
        nlevels <- dim(covmat)[1]
 	firstRow <- covmat[1, ]
 	ones <- rep(1, nlevels)
 	J <- outer(ones, ones)
 	notzero <- 2:nlevels
 	r.covmat <- covmat + (firstRow[1]*J) -
 		             outer(firstRow, ones) -
 	   		     outer(ones, firstRow)
 	r.covmat <- r.covmat[notzero, notzero]
 	qv1 <- qvmat[1, 1]
 	r.qvmat <- (qvmat + qv1*J)[notzero, notzero]
 	list(r.covmat, r.qvmat)}
    covmat <- qv.object$covmat
    qvmat <- diag(qv.object$qvframe$quasiVar)
    r.form <- reducedForm(covmat, qvmat)
    r.covmat <- r.form[[1]]
    r.qvmat <- r.form[[2]]
    inverse.sqrt <- solve(chol(r.covmat))
    evalues <- eigen(t(inverse.sqrt) %*% r.qvmat %*% inverse.sqrt,
 	               symmetric=TRUE)$values
    sqrt(c(min(evalues), max(evalues))) - 1
}

indentPrint <- function(object, indent = 4, ...){
    zz <- ""
    tc <- textConnection("zz", "w", local = TRUE)
    sink(tc)
    try(print(object, ...))
    sink()
    close(tc)
    indent <- paste(rep(" ", indent), sep = "", collapse = "")
    cat(paste(indent, zz, sep = ""), sep = "\n")}

print.qv <- function(x, ...){
    print(x$qvframe)
}

summary.qv <- function(object, ...)
{
    if (!is.null(object$modelcall))
        cat("Model call: ",
            deparse(object$modelcall), "\n")
    if (!is.null(object$dispersion))
        cat("Dispersion: ", object$dispersion, "\n")
    if (!is.null(object$factorname))
        cat("Factor name: ",object$factorname,"\n")
    indentPrint(object$qvframe,...)
    if (!is.null(object$relerrs)){
        minErrSimple <- round(100*min(object$relerrs), 1)
        maxErrSimple <- round(100*max(object$relerrs), 1)
        errors<-worstErrors(object)
        minErrOverall<-round(100*errors[1], 1)
        maxErrOverall<-round(100*errors[2], 1)
        cat("Worst relative errors in SEs of simple contrasts (%): ",
            minErrSimple, maxErrSimple, "\n")
        cat("Worst relative errors over *all* contrasts (%): ",
            minErrOverall, maxErrOverall, "\n")
    }
}

plot.qv <- function(x,
                    intervalWidth = 2,
                    ylab = "estimate",
                    xlab = x$factorname,
                    ylim = NULL,
                    main = "Intervals based on quasi standard errors",
                    levelNames = NULL,
                    ...) {
    frame <- x$qvframe
    if (!is.null(levelNames)) {
        if (nrow(frame) != length(levelNames)) stop(
                "levelNames is not a vector of the right length"
                )
        row.names(frame) <- levelNames
    }
    if (is.null(frame$quasiSE))
        stop("Cannot plot, because there are no quasi standard errors")
    if (is.na(frame$estimate[1]))
        stop("No parameter estimates to plot")
    if (any(is.nan(frame$quasiSE)))
        stop(paste("No comparison intervals available,\n",
                   "since one of the quasi variances is negative.",
                   "  See ?qvcalc for more.",
                   sep = ""))
    faclevels <- factor(row.names(frame), levels = row.names(frame))
    xvalues <- seq(along = faclevels)
    est <- frame$estimate
    se <- frame$quasiSE
    tops <- est + (intervalWidth * se)
    tails <- est - (intervalWidth * se)
    range <- max(tops) - min(tails)
    if (is.null(ylim)) ylim <- c(min(tails) - range/10, max(tops) + range/10)
    if (is.null(xlab)) xlab <- "factor level"
    plot(faclevels, frame$estimate, border = "transparent", ylim = ylim,
         xlab = xlab, ylab = ylab,
         main = main, ...)
    points(frame$estimate, ...)
    segments(xvalues, tails, xvalues, tops)
    invisible(x)
}

