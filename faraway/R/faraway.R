# Functions for use with Faraway's book
#
"maxadjr" <-
  function (l,best=3)
# Display the best (3) models from a leaps() object
{
  i <- rev(order(l$a))
  nopreds <- max(l$size)-1
  labels <- apply(l$which,1,function(x) paste(as.character((1:nopreds)[x]),collapse=","))
  m <- round(l$a[i[1:best]],3)
  names(m) <- labels[i[1:best]]
#  m <- cbind(round(l$a[i[1:best]],3),labels[i[1:best]])
#  dimnames(m) <- list(NULL,c("Adj R^2","Model"))
  m
}
"qqnorml" <-
  function(y,main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles",
    ylab = "Sample Quantiles",...)
# labeled Q-Q plot
  {
    n <- length(y)
    u <- qnorm((1:n)/(n+1))
    i <- order(y)
    plot(u,y[i],xlab=xlab,ylab=ylab,main=main,type="n")
    text(u,y[i],as.character(1:n)[i])
  }
"Cpplot" <-
  function (cp)
# Construct a Cp plot
{
  p <- max(cp$size)
  i <- (cp$Cp < (p+1.5))
  plot(cp$size[i],cp$Cp[i],xlab="p",ylab="Cp",type="n")
  labels <- apply(cp$which,1,function(x) paste(as.character((1:(p-1))[x]),collapse=""))
  text(cp$size[i],cp$Cp[i],labels[i])
  abline(0,1)
}
vif <- function(object)
UseMethod("vif")

vif.default <- function(object) {
  if(!is.data.frame(object) & !is.matrix(object)) stop("Not matrix or data frame")
  if(is.data.frame(object)) object <- as.matrix(object)
  ncols <- dim(object)[2]
  v <- numeric(ncols)
  names(v) <- dimnames(object)[[2]]
  for(i in 1:ncols) v[i] <- 1/(1-summary(lm(object[,i]~object[,-i]))$r.squared)
  v
}

# function from Bill Venables post on R-digest
vif.lm <- function(object) {
  V <- summary(object)$cov.unscaled
  Vi <- crossprod(model.matrix(object))
        nam <- names(coef(object))
  if(k <- match("(Intercept)", nam, nomatch = FALSE)) {
                v1 <- diag(V)[-k]
                v2 <- (diag(Vi)[-k] - Vi[k, -k]^2/Vi[k,k])
                nam <- nam[-k]
        } else {
                v1 <- diag(V)
                v2 <- diag(Vi)
                warning("No intercept term detected.  Results may surprise.")
        }
        structure(v1*v2, names = nam)
}

prplot <- function(g,i)
{
# Partial residuals plot for predictor i
  xl <- attributes(g$terms)$term.labels[i]
  yl <- paste("beta*",xl,"+res",sep="")
  x <- model.matrix(g)[,i+1]
  plot(x,g$coeff[i+1]*x+g$res,xlab=xl,ylab=yl)
  abline(0,g$coeff[i+1])
  invisible()
}
"halfnorm" <-
function (x, nlab = 2, labs = as.character(1:length(x)), ylab = "Sorted Data",
            ...)
{
  x <- abs(x)
  labord <- order(x)
  x <- sort(x)
  i <- order(x)
  n <- length(x)
  ui <- qnorm((n + 1:n)/(2 * n + 1))
  plot(ui, x[i], xlab = "Half-normal quantiles", ylab = ylab, ylim=c(0,max(x)),
       type = "n", ...)
  if(nlab < n)
    points(ui[1:(n - nlab)], x[i][1:(n - nlab)])
  text(ui[(n - nlab + 1):n], x[i][(n - nlab + 1):n], labs[labord][(n -
                                                              nlab + 1):n])
}
# logit and inverse logit
logit <- function(x){
  if(any(omit <- (is.na(x) | x <=0 | x >= 1))){
    lv <- x
    lv[omit] <- NA
    if(any(!omit))
      lv[!omit] <- Recall(x[!omit])
    return(lv)
  }
  log(x/(1-x))
}
ilogit <- function(x){
  if(any(omit <- is.na(x))){
    lv <- x
    lv[omit] <- NA
    if(any(!omit))
      lv[!omit] <- Recall(x[!omit])
    return(lv)
  }
  exp(x)/(1 + exp(x))
}

# Essential regression summary (idea from Gelman and Hill)

if (!isGeneric("sumary")) {
    setGeneric("sumary",
               function(object, ...)
               standardGeneric("sumary"))
}

setMethod("sumary", signature(object = "lm"),
    function(object)
    {
  digits <- options()$digits
  summ <- summary (object)
  sigma.hat <- summ$sigma
  r.squared <- summ$r.squared
  coef <- summ$coef[,,drop=FALSE]
  n <- summ$df[1] + summ$df[2]
  p <- summ$df[1]
  if (nsingular <- summ$df[3] - summ$df[1]) cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", sep = "")
  printCoefmat(coef,signif.stars=FALSE)
  cat("\n")
  cat (paste ("n = ", n, ", p = ", p,
    ", Residual SE = ", format(round(sigma.hat, digits-2),nsmall=digits-2),
    ", R-Squared = ", format(round(r.squared, 2)), "\n", sep=""))
  invisible(summ)
    }
)

setMethod("sumary", signature(object = "glm"),
    function(object, dispersion=NULL)
    {
        digits <- options()$digits
        summ <- summary(object, dispersion = dispersion)
        n <- summ$df[1] + summ$df[2]
        p <- summ$df[1]
        coef <- summ$coef[,,drop=FALSE]
        printCoefmat(coef,signif.stars=FALSE)
        cat("\n")
        if (summ$dispersion != 1) {
            cat(paste0("Dispersion parameter = ", fround(summ$dispersion,digits-2),"\n"))
        }
        cat(paste0("n = ", n, " p = ", p,"\n"))
        cat(paste0("Deviance = ",fround(summ$deviance,digits-2),
                   " Null Deviance = ", fround(summ$null.deviance,digits-2),
                   " (Difference = ", fround(summ$null.deviance-summ$deviance,digits-2), ")"),"\n")
    return(invisible(summ))
  }
)

# This one taken from the arm package. Edited out the model printing.

setMethod("sumary", signature(object = "merMod"),
    function(object, digits=2, detail=FALSE)
    {
    out <- NULL
    out$call <- object@call
    #print (out$call)
    #object <- summary(object)
    #summ <- summary(object)
    fcoef <- fixef(object)
    #coefs <- attr(summ, "coefs")
    #useScale <- attr (VarCorr (object), "sc")
    useScale <- getME(object, "devcomp")$dims["useSc"]
    corF <- vcov(object)@factors$correlation
    coefs <- cbind(fcoef, corF@sd)
    if (length (fcoef) > 0){
      if (!useScale) {
        coefs <- coefs[, 1:2, drop = FALSE]
        out$z.value <- coefs[, 1]/coefs[, 2]
        out$p.value <- 2 * pnorm(abs(out$z.value), lower.tail = FALSE)
        coefs <- cbind(coefs, `z value` = out$z.value, `Pr(>|z|)` = out$p.value)
      }
      else {
        out$t.value <- coefs[, 1]/coefs[, 2]
        coefs <- cbind(coefs, `t value` = out$t.value)
      }
    dimnames(coefs)[[2]][1:2] <- c("coef.est", "coef.se")
      cat("Fixed Effects:\n")
      if(detail){
        pfround (coefs, digits)
      }
      else{
        pfround(coefs[,1:2], digits)
      }
    }
    out$coef <- coefs[,"coef.est"]
    out$se <- coefs[,"coef.se"]
    cat("\nRandom Effects:\n")
    vc <- as.matrix.VarCorr (VarCorr (object), useScale=useScale, digits)
    print (vc[,c(1:2,4:ncol(vc))], quote=FALSE)
    out$ngrps <- lapply(object@flist, function(x) length(levels(x)))
    is_REML <- isREML(object)
    llik <- logLik(object, REML=is_REML)
    out$AIC <- AIC(llik)
    out$deviance <- deviance(refitML(object))     # Dbar
    out$n <- getME(object, "devcomp")$dims["n"]
    Dhat <- -2*(llik) # Dhat
    pD <- out$deviance - Dhat              # pD
    out$DIC <- out$deviance + pD               # DIC=Dbar+pD=Dhat+2pD
    cat("---\n")
    cat(sprintf("number of obs: %d, groups: ", out$n))
    cat(paste(paste(names(out$ngrps), out$ngrps, sep = ", "), collapse = "; "))
    cat(sprintf("\nAIC = %g, DIC = ", round(out$AIC,1)))
    cat(round(out$DIC, 1))
    cat("\ndeviance =", fround (out$deviance, 1), "\n")
    if (useScale < 0){
      out$sigma.hat <- .Call("mer_sigma", object, FALSE, PACKAGE = "lme4")
      cat("overdispersion parameter =", fround (out$sigma.hat, 1), "\n")
    }
    return(invisible(out))
    }
)

# Taken from arm package

fround <- function (x, digits) {
    format (round (x, digits), nsmall=digits)
}
pfround <- function (x, digits) {
    print (fround (x, digits), quote=FALSE)
}


# Taken from arm package

as.matrix.VarCorr <- function (varc, useScale, digits){
# VarCorr function for lmer objects, altered as follows:
#   1.  specify rounding
#   2.  print statement at end is removed
#   3.  reMat is returned
#   4.  last line kept in reMat even when there's no error term
                  sc <- attr(varc, "sc")[[1]]
                  if(is.na(sc)) sc <- 1
#                  recorr <- lapply(varc, function(el) el@factors$correlation)
                  recorr <- lapply(varc, function(el) attr(el, "correlation"))
                  #reStdDev <- c(lapply(recorr, slot, "sd"), list(Residual = sc))
                  reStdDev <- c(lapply(varc, function(el) attr(el, "stddev")), list(Residual = sc))
                  reLens <- unlist(c(lapply(reStdDev, length)))
                  reMat <- array('', c(sum(reLens), 4),
                                 list(rep('', sum(reLens)),
                                      c("Groups", "Name", "Variance", "Std.Dev.")))
                  reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
                  reMat[,2] <- c(unlist(lapply(reStdDev, names)), "")
#                  reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
#                  reMat[,4] <- format(unlist(reStdDev), digits = digits)
                  reMat[,3] <- fround(unlist(reStdDev)^2, digits)
                  reMat[,4] <- fround(unlist(reStdDev), digits)
                  if (any(reLens > 1)) {
                      maxlen <- max(reLens)
                      corr <-
                          do.call("rbind",
                                  lapply(recorr,
                                         function(x, maxlen) {
                                             x <- as(x, "matrix")
#                                             cc <- format(round(x, 3), nsmall = 3)
                                             cc <- fround (x, digits)
                                             cc[!lower.tri(cc)] <- ""
                                             nr <- dim(cc)[1]
                                             if (nr >= maxlen) return(cc)
                                             cbind(cc, matrix("", nr, maxlen-nr))
                                         }, maxlen))
                      colnames(corr) <- c("Corr", rep("", maxlen - 1))
                      reMat <- cbind(reMat, rbind(corr, rep("", ncol(corr))))
                  }
#                  if (!useScale) reMat <- reMat[-nrow(reMat),]
          if (useScale<0) reMat[nrow(reMat),] <- c ("No residual sd", rep("",ncol(reMat)-1))
          return (reMat)
      }
