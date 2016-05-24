#' The adaptive fence procedure
#'
#' This function implements the adaptive fence procedure to
#' first find the optimal cstar value and then finds the
#' corresponding best model as described in Jiang et. al.
#' (2009) with some practical modifications.
#'
#' The initial stepwise procedure performs forward stepwise model
#' selection using the AIC and backward stepwise model selection
#' using BIC.  In general the backwise selection via the more
#' conservative BIC will tend to select a smaller model than that
#' of the forward selection AIC approach.  The size of these two
#' models is found, and we go two dimensions smaller and larger
#' to estimate a sensible range of \code{c} values over which to
#' perform a parametric bootstrap.
#'
#' This procedure can take some time.  It is recommended that you start
#' with a relatively small number of bootstrap samples (\code{B})
#' and grid of boundary values (\code{n.c}) and increase both as
#' required.
#'
#' If you use \code{initial.stepwise=TRUE} then in general you will
#' need a smaller grid of boundary values than if you select
#' \code{initial.stepwise=FALSE}.
#' It can be useful to check \code{initial.stepwise=FALSE} with a
#' small number of bootstrap replications over a sparse grid to ensure
#' that the \code{initial.stepwise=TRUE} has landed you in a reasonable
#' region.
#'
#' The \code{best.only=FALSE} option when plotting the results of the
#' adaptive fence is a modification to the adaptive fence procedure
#' which considers all models at a particular size that pass the fence
#' hurdle when calculating the p* values.  In particular,
#' for each value of c and at each bootstrap replication,
#' if a candidate model is found that passes the fence, then we look to see
#' if there are any other models of the same size that also pass the fence.
#' If no other models of the same size pass the fence, then that model is
#' allocated a weight of 1.  If there are two models that pass the fence, then
#' the best model is allocated a weight of 1/2.  If three models pass the fence,
#' the best model gets a weight of 1/3, and so on. After \code{B} bootstrap
#' replications, we aggregate the weights by summing over the various models.
#' The p* value is the maximum aggregated weight divided by the number of bootstrap
#' replications.
#' This correction penalises the probability associated with the best model if
#' there were other models of the same size that also passed the fence hurdle.
#' The rationale being that if a model has no redundant variables
#' then it will be the only model at that size that passes the fence over a
#' range of values of c.
#' The result is more pronounced peaks which can help to determine
#' the location of the correct peak and identify the optimal c*.
#'
#' @param mf a fitted 'full' model, the result of a call
#'   to lm or glm (and in the future lme or lmer).
#' @param B number of bootstrap replications at each fence
#'   boundary value
#' @param n.c number of boundary values to be considered
#' @param initial.stepwise logical.  Performs an initial stepwise
#'   procedure to look for the range of model sizes where attention
#'   should be focussed. See details for implementation.
#' @param force.in the names of variables that should be forced
#'   into all estimated models
#' @param n.cores number of cores to be used when parallel
#'   processing the bootstrap
#' @param nvmax size of the largest model that can still be
#'   considered as a viable candidate.  Included for performance
#'   reasons but if it is an active constraint it could lead to
#'   missleading results.
#' @param c.max manually specify the upper boundary limit.
#'   Only applies when \code{initial.stepwise=FALSE}.
#' @param screen logical, whether or not to perform an initial
#'   screen for outliers.  Highly experimental, use at own risk.
#'   Default = FALSE.
#' @param ... further arguments (currently unused)
#' @references Jiang J., Nguyen T., Sunil Rao J. (2009),
#'   A simplified adaptive fence procedure, Statistics &
#'   Probability Letters, 79(5):625-629.  doi: 10.1016/j.spl.2008.10.014
#'
#'   Jiang J., Sunil Rao J., Gu Z, Nguyen T. (2008),
#'   Fence methods for mixed model selection, Annals of Statistics,
#'   36(4):1669-1692. doi: 10.1214/07-AOS517
#' @export
#' @import foreach
#' @import parallel
#' @family fence
#' @examples
#' n = 100
#' set.seed(11)
#' e = rnorm(n)
#' x1 = rnorm(n)
#' x2 = rnorm(n)
#' x3 = x1^2
#' x4 = x2^2
#' x5 = x1*x2
#' x6 = rep(c("A","B"),n=50)
#' y = 1 + x1 + x2 + e
#' dat = data.frame(y,x1,x2,x3,x4,x5,x6)
#' lm1 = lm(y~.,data=dat)
#' \dontrun{
#' af1 = af(lm1, n.cores=4, initial.stepwise=TRUE)
#' summary(af1)
#' plot(af1)
#' }

af = function(mf,
              B=60, n.c=20,
              initial.stepwise=FALSE,
              force.in=NULL,
              n.cores, nvmax, c.max, screen=FALSE,  ...){
  method="ML"
  af.call = match.call()
  if(!missing(c.max) & initial.stepwise==TRUE) {
    initial.stepwise=FALSE
    warning("When c.max is specified, initial.stepwise=FALSE")
  }
  if(!any(class(mf)=="lm")){
    warning("Adaptive fence currently only implemented for lm and glm")
    model.type = "lm"
  }
  if(any(class(mf)=="glm")==TRUE){
    family=stats::family(mf)
    if(!is.null(force.in)){
      warning("force.in is not implemented for glms")
    }
    model.type="glm"
  } else if(class(mf)=="lm"){
    model.type="lm"
  }
  m = mextract(mf,screen=screen)
  fixed = m$fixed
  yname = m$yname
  family = m$family
  Xy = m$X
  kf = m$k
  n = m$n
  Xy$initial.weights = m$wts
  initial.weights = m$wts
  null.ff = stats::as.formula(paste(yname,"~1"))
  if(model.type=="glm"){
    m0 = stats::glm(null.ff, data = Xy, family=family, weights = initial.weights)
    mfstar = stats::glm(fixed, data = Xy, family=family, weights = initial.weights)
  } else {
    m0 = stats::lm(null.ff, data = Xy, weights = initial.weights)
    mfstar = stats::lm(fixed, data = Xy, weights = initial.weights)
  }
  Qm0 = Qm(m0, method=method)
  Qmfstar = Qm(mfstar, method=method)
  if(!is.null(force.in)){
    small.ff = stats::as.formula(paste(yname,"~",paste(force.in,collapse="+")))
  } else {
    small.ff = null.ff
  }
  if (initial.stepwise) {
    if(model.type=="glm"){
      small.mod = stats::glm(small.ff, data = Xy, family=family, weights = initial.weights)
    } else {
      small.mod = stats::lm(small.ff, data = Xy, weights = initial.weights)
    }
    # backwards and forwards model selection using
    # BIC (conservative) and AIC (less conservative)
    bwds.BIC = stats::step(mfstar,
                           scope = list(lower=small.ff, upper=fixed),
                           direction="backward", k=log(n), trace=0)
    fwds.BIC = stats::step(small.mod,
                           scope = list(lower=small.ff, upper=fixed),
                           direction="forward", k=log(n), trace=0)
    bwds.AIC = stats::step(mfstar,
                           scope = list(lower=small.ff, upper=fixed),
                           direction="backward", k=2, trace=0)
    fwds.AIC = stats::step(small.mod,
                           scope = list(lower = small.ff, upper = fixed),
                           direction = "forward", k = 2, trace = 0)
    k.vals = c(length(bwds.BIC$coef),length(fwds.BIC$coef),
               length(bwds.AIC$coef),length(fwds.AIC$coef))
    k.min = max(min(k.vals) - 2, 1)
    k.max = min(max(k.vals) + 2, kf)
    k.range = list(k.min=k.min,k.max=k.max)
    Q.range = qrange(k.range=k.range, data=Xy,
                     yname=yname, fixed=fixed,
                     method=method, force.in=force.in,
                     model.type = model.type,
                     family = family)
    c.max = (Q.range$Q.max - Qmfstar)*1.1
    c.min = max(Q.range$Q.min - Qmfstar,0)*0.9
    c.range = seq(c.min,c.max,length.out=n.c)
    if(missing(nvmax)) nvmax = k.max
  } else {
    k.range = list(k.min=1,k.max=kf)
    if(missing(c.max)) c.max = (Qm0-Qmfstar)*1.1
    c.min = 0.1
    c.range = seq(c.min,c.max,length.out=n.c)
    if(missing(nvmax)) nvmax = kf
  }
  
  if(missing(n.cores)) n.cores = max(detectCores()-1,1)
  cl.af = makeCluster(n.cores)
  doParallel::registerDoParallel(cl.af)
  j=NULL # avoid global variable NOTE in R CMD check
  p.star.all = foreach(j = 1:n.c, .combine=rbind, .packages = c("mplot")) %dopar% {
    fence.mod = list()
    fence.rank = list()
    ystar = stats::simulate(object=mfstar,nsim=B)
    initial.weights <<- m$wts
    if(model.type=="glm"){
      for(i in 1:B){
        Xy[yname] = ystar[,i]
        mfstarB = do.call("glm",list(fixed,data=Xy,family=family, weights = initial.weights))
        fms = glmfence(mfstarB, cstar=c.range[j],
                       trace=FALSE, nvmax=nvmax, adaptive=TRUE)
        fence.mod = c(fence.mod,fms)
        fence.rank = c(fence.rank,1:length(fms))
      }
    } else {
      for(i in 1:B){
        Xy[yname] = ystar[,i]
        mfstarB = do.call("lm",list(fixed, data=Xy, weights = initial.weights))
        fms = lmfence(mfstarB, cstar=c.range[j], trace=FALSE,
                      nvmax = nvmax, force.in=force.in, adaptive=TRUE)
        fence.mod = c(fence.mod,fms)
        fence.rank = c(fence.rank,1:length(fms))
      }
    }
    process.fn(fence.mod,fence.rank)
  }
  stopCluster(cl.af)
  
  # Another function that processes results within af function
  #
  # This function is used by the af function to process
  # the results when iterating over different boundary values
  
  pstar.fn = function(input,type){
    if(type=="bo"){
      p.star = input[,1:2]
    } else if(type=="all"){
      p.star = input[,3:4]
    }
    pstarmods = sort(table(p.star[,2]),decreasing=TRUE)
    n.pstarmods = length(pstarmods)
    p.star = data.frame(pstar = as.numeric(p.star[,1]),
                        model = as.character(p.star[,2]),
                        modelident = match(p.star[,2], names(pstarmods)))
    redundent.vars = grepl("REDUNDANT.VARIABLE",as.character(p.star$model))
    c.range = c.range[!redundent.vars]
    p.star = p.star[!redundent.vars,]
    p.star$model = droplevels(p.star$model)
    # if want runs of (near) maximums
    max.p = max(p.star[,1]) #- 2/B
    tf = p.star[,1] >= max.p
    a = rle(tf)
    pos = which(a$values==TRUE)
    # what if there was a sequence of falses of the same length?
    # keep only the position of these runs where we had true values
    pos = pos[a$values[pos]==TRUE]
    mid = NA
    for(i in 1:length(pos)){
      # find the midpoint
      if(pos[i]==1){
        mid[i] = sum(a$lengths[1:pos[i]]+1)/2
      } else {
        mid[i] = (sum(a$length[1:pos[i]]) + sum(a$lengths[1:(pos[i]-1)])+1)/2
      }
    }
    mid = floor(mid)
    c.star = min(c.range[mid])
    if(model.type=="glm"){
      afmod = glmfence(mf, cstar=c.star,
                       trace=FALSE, nvmax=nvmax)[[1]]
    } else {
      afmod = lmfence(mf, cstar=c.star, trace=FALSE,
                      nvmax=nvmax, force.in=force.in)[[1]]
    }
    p.star[,1] = as.numeric(as.character(p.star[,1]))
    return(list(p.star=p.star,
                c.range=c.range,
                c.star=c.star,
                model = afmod))
  }
  # set up the output object class
  afout = list()
  afout$bestOnly = pstar.fn(p.star.all,type="bo")
  afout$all = pstar.fn(p.star.all,type="all")
  afout$call = af.call
  afout$screen = screen
  if(initial.stepwise){
    afout$initial.stepwise = list(
      fwds.AIC = stats::as.formula(fwds.AIC),
      fwds.BIC = stats::as.formula(fwds.BIC),
      bwds.AIC = stats::as.formula(bwds.AIC),
      bwds.BIC = stats::as.formula(bwds.BIC))
  } else {afout$initial.stepwise=NULL}
  afout$k.range = k.range
  class(afout) = "af"
  return(afout)
}

#' Summary method for an af object
#'
#' Provides comprehensive  output of the bootstrap results of an
#' af object.
#'
#' @param object \code{af} object, the result of \code{\link{af}}
#' @param best.only logical determining whether the output used the
#'   standard fence approach of only considering the best models
#'   that pass the fence (\code{TRUE}) or if it should take into
#'   account all models that pass the fence at each boundary
#'   value (\code{FALSE}).
#' @param ... further arguments (currently unused)
#' @export
# S3 method for class 'af'
summary.af = function (object,best.only=TRUE,...) {
  if(best.only){
    xsub = object$bestOnly
  } else {
    xsub = object$all
  }
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Adaptive fence model (c*=")
  cat(round(xsub$c.star,1))
  cat("):\n")
  cat(deparse(xsub$model))
  cat("\n\n")
  cat("Model sizes considered: ")
  cat(object$k.range$k.min)
  cat(" to ")
  cat(object$k.range$k.max)
  cat(" (including intercept).")
  cat("\n\n")
  if(!is.null(object$initial.stepwise)){
    cat("Stepwise procedures:\n")
    cat("Forwards AIC: ")
    cat(deparse(object$initial.stepwise$fwds.AIC))
    cat("\n")
    cat("Backwards AIC: ")
    cat(deparse(object$initial.stepwise$bwds.AIC))
    cat("\n")
    cat("Forwards BIC: ")
    cat(deparse(object$initial.stepwise$fwds.BIC))
    cat("\n")
    cat("Backwards BIC: ")
    cat(deparse(object$initial.stepwise$bwds.BIC))
    cat("\n\n")
  }
  invisible(object)
}


#' Plot diagnostics for an af object
#'
#' Summary plot of the bootstrap results of an
#' af object.
#'
#' @param x \code{af} object, the result of \code{\link{af}}
#' @param classic logical.  If \code{classic=TRUE} a
#'   base graphics plot is provided instead of a googleVis plot.
#'   Default is \code{classic=FALSE}.
#' @param best.only logical determining whether the output used the
#'   standard fence approach of only considering the best models
#'   that pass the fence (\code{TRUE}) or if it should take into
#'   account all models that pass the fence at each boundary
#'   value (\code{FALSE}).
#' @param pch plotting character, i.e., symbol to use
#' @param tag Default NULL. Name tag of the objects to be extracted 
#' from a gvis (googleVis) object. 
#' 
#' The default tag for is NULL, which will 
#' result in R opening a browser window.  Setting \code{tag='chart'} 
#' or setting \code{options(gvis.plot.tag='chart')} is useful when 
#' googleVis is used in scripts, like knitr or rmarkdown. 
#' 
#' @param shiny Default FALSE. Set to TRUE when using in a shiny interface.
#' 
#' @param width Width of the googleVis chart canvas area, in pixels.
#'   Default: 800.
#' @param height Height of the googleVis chart canvas area, in pixels.
#'   Default: 400.
#' @param chartWidth googleVis chart area width.
#'   A simple number is a value in pixels;
#'   a string containing a number followed by \code{\%} is a percentage.
#'   Default: \code{"60\%"}
#' @param chartHeight googleVis chart area height.
#'   A simple number is a value in pixels;
#'   a string containing a number followed by \code{\%} is a percentage.
#'   Default: \code{"80\%"}
#' @param fontSize font size used in googleVis chart.  Default: 12.
#' @param left space at left of chart (pixels?).  Default: "50".
#' @param top space at top of chart (pixels?).  Default: "30".
#' @param options If you want to specify the full set of googleVis
#'   options.
#' @param backgroundColor The background colour for the main area
#'   of the chart. A simple HTML color string,
#'   for example: 'red' or '#00cc00'.  Default: 'transparent'
#' @param legend.position legend position, e.g. \code{"topleft"}
#'   or  \code{"bottomright"}
#' @param ... further arguments (currently unused)
#' @export
# S3 method for class 'af'
plot.af = function(x, pch, classic = FALSE,
                   tag = NULL, shiny = FALSE,
                   best.only = TRUE,
                   width = 800, height = 400, fontSize = 12,
                   left = 50, top = 30, chartWidth = "60%",
                   chartHeight = "80%",
                   backgroundColor = 'transparent',
                   legend.position = "topleft",
                   options=NULL, ...){
  if (best.only) {
    x = x$bestOnly
  } else {
    x = x$all
  }
  if(classic){
    if(missing(pch)) pch=19
    graphics::par(mar=c(3.4,3.4,2.1,0.1),mgp=c(2.0, 0.75, 0))
    graphics::plot(x$p.star[,1]~x$c.range,
                   ylim=c(0,1), pch=pch,
                   col=x$p.star[,3],
                   ylab = "p*", xlab = "c")
    graphics::legend(legend.position, legend=unique(x$p.star[,2]),
                     pch=pch,col=unique(x$p.star[,3]),bty="n")
    graphics::axis(side=3, at=x$c.star,
                   labels=paste("c*=", round(x$c.star,1),sep=""))
  } else {
    dat <- matrix(NA, nrow = nrow(x$p.star),
                  ncol = nlevels(x$p.star$model) + 1)
    for(i in 1:nlevels(x$p.star$model)){
      lvl <- levels(x$p.star$model)[i]
      ind <- which(x$p.star$model == lvl)
      dat[ind, c(1, i+1)] <- x$p.star$pstar[ind]
    }
    plot.dat = data.frame(c.range = as.numeric(x$c.range),
                          dat[,-1])
    colnames(plot.dat) = c("c.range",levels(x$p.star$model))
    plot.dat = round(plot.dat,2)
    # FOR FUN ON A RAINY DAY
    # INCLUDE ANNOTATION USING `ROLES'
    # SEE HERE: http://cran.r-project.org/web/packages/
    # googleVis/vignettes/Using_Roles_via_googleVis.html
    gvis.title = paste("Adaptive fence: c*=",round(x$c.star,1),sep="")
    namefunc <- function(v1) {
      deparse(substitute(v1))
    }
    chartArea = paste("{left:",left,
                      ",top:",top,
                      ",width:'",chartWidth,
                      "',height:'",chartHeight,"'}",sep="")
    if(is.null(options)){
      options=list(title=gvis.title,
                   fontSize=fontSize,
                   vAxis="{title:'p*',minValue:0,maxValue:1,
                           ticks: [0.0,0.2,0.4,0.6,0.8,1.0]}",
                   hAxis="{title:'c'}",
                   axisTitlesPosition="out",
                   chartArea=chartArea,
                   backgroundColor=backgroundColor,
                   width=width, height=height)
    }
    fplot = googleVis::gvisScatterChart(data=plot.dat,options=options)
    if(shiny){
      return(fplot)
    } else {
      return(graphics::plot(fplot, tag = tag))
    }
  }
}


#' Print method for an af object
#'
#' Prints basic output of the bootstrap results of an
#' af object.
#'
#' @param x an \code{af} object, the result of \code{\link{af}}
#' @param best.only logical determining whether the output used the
#'   standard fence approach of only considering the best models
#'   that pass the fence (\code{TRUE}) or if it should take into
#'   account all models that pass the fence at each boundary
#'   value (\code{FALSE}).
#' @param ... further arguments (currently unused)
#' @export
# S3 print method for class 'af'
print.af = function (x, best.only=TRUE, ...) {
  if(best.only){
    x = x$bestOnly
  } else {
    x = x$all
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Adaptive fence model (c*=")
  cat(round(x$c.star,1))
  cat("):\n")
  cat(deparse(x$model))
  cat("\n\n")
  invisible(x)
}
