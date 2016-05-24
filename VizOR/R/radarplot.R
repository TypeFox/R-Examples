##' A multidimensional 'Table 1 at-a-glance'
##'
##' Plots a specialized radar (aka, \sQuote{star}) plot to compare summary
##' statistics on a handful of covariates, across a multi-dimensional set
##' of comparison groups, as in the standard \sQuote{Table 1} of experimental
##' studies, where the groups correspond to treatment assignments.
##'
##' The covariates appear on the radii, or 'spokes' of the radarplot,
##' while the several (up to 3) dimensions defining the comparison
##' groups appear as a color-coded overlay, plus up to two trellis
##' dimensions.  In generic applications, the client code employs
##' \code{Hmisc::summary.formula(method='cross')} to create the
##' multi-dimensional array for plotting, and then rescales the array
##' elements by a linear or affine transformation to the interval
##' [0,1].  Ordered factors may be treated as quasi-continuous
##' variables according to their internal integer representation;
##' logical covariates (and also, conveniently, 2-valued factors)
##' appear as proportions; and unordered categorical variables with n
##' levels may be transformed to n-1 dummy variables, which may then
##' be plotted as proportions.  An optional \code{strength} parameter
##' may be specified, providing a strength-of-association measure
##' connecting the spoke variables to the experimental outcome.  This
##' is reflected graphically by fanning the spokes out into sectors
##' subtending angles proportional to the strength of associations.
##' Used in this way, the radar plot draws attention to those spoke
##' variables with greatest potential to confound the outcome, and so
##' becomes a useful exploratory tool.  To provide convenient support
##' for this specialized usage, \code{radarplot} optionally accepts a
##' regression model in which the outcome of interest is the
##' regressand, and the spoke variables are regressors.  From this
##' model, strength-of-association measures are automatically derived,
##' appropriately to the class of the regression model.  Invoked in
##' this way, \code{radarplot} is also able to free client code from
##' the burden of generating and rescaling a \code{summary.formula.cross}
##' object. At this time, models of class \code{lrm} and \code{cph} are
##' supported. See the example code below.
##'
##' @param x Typically, a formula of the form \code{S~x}, \code{S~x|t}
##' or \code{S~x|t*u}, where \code{S} is the (suitably normalized)
##' summary matrix of a \code{summary.formula.cross} object having factors
##' \{\code{x}[, \code{t}[, \code{u}]]\} as covariates. In the resulting
##' radar plot, the \code{x} factor will correspond to colored polygons
##' overlaid within each panel of a trellis of 0, 1 or 2 dimensions, defined
##' by the (optional) terms \code{t} and \code{u}.  Alternatively, \code{x}
##' may be a fitted model in which the LHS is the experimental outcome,
##' and the RHS includes the spoke variables.
##' @param data Either a \code{summary.formula.cross} object (when
##' \code{x} is a formula), or a \code{data.frame} (when \code{x} is a
##' fitted model).
##' @param datadist An \code{rms:datadist} object describing the
##' underlying dataset. By suitably modifying \code{datadist$limits},
##' greater control may be achieved over the automatic calculation of
##' strength-of-association measures.
##' @param xlim A positioning parameter that really should be automated!
##' @param ylim A positioning parameter that really should be automated!
##' @param rescale How to rescale continuous covariates for plotting on the spokes
##' @param treatment When \code{x} is an outcomes model, identifying
##' the treatment variable (and any variables to stratify on, see
##' below) permits the formula embedded in the model to be used to
##' generate a \code{summary.formula.cross} object
##' @param stratify When \code{x} is an outcomes model, identifying
##' the treatment variable (see above) and any variables to stratify
##' on, permits the formula embedded in the model to be used to
##' generate a \code{summary.formula.cross} object
##' @param strength When the default choice of strength-of-association
##' measure is not suitable, it may be specified explicitly as a named
##' vector with elements corresponding to the spoke variables, in the
##' order they appear in the model formula. The names will be printed
##' at the outer endpoints of the spokes, and will typically be
##' numbers (e.g., odds ratios) rounded to 2-3 significant figures.
##' @param include.na Should the summary include NA values? Passed to
##' the summary.formula call that generates the radarplot data.
##' @param overall Should the summary include an 'ALL' row? Passed to
##' the summary.formula call that generates the radarplot data.
##' @param ... Other parameters to be passed to \code{panel.radarplot}
##' @keywords hplot htest
##'
##' @examples
##' library(rms)
##' df <- upData(mtcars,
##'              cyl=factor(cyl,levels=2*(2:4),labels=paste(2*(2:4),"cyl", sep="-")),
##'              am=factor(am,levels=0:1,labels=c("automatic","manual")),
##'              gear=factor(gear,levels=3:5,labels=paste(3:5,"speed", sep="-")),
##'              labels=c(
##'                mpg="Miles per gallon"
##'                ,cyl="Number of cylinders"
##'                ,disp="Displacement"
##'                ,hp="Gross horsepower"
##'                ,drat="Rear axle ratio"
##'                ,wt="Weight"
##'                ,qsec="1/4 mile time"
##'                ,am="Transmission type"
##'                ,gear="Number of forward gears"
##'                ,carb="Number of carburetors"
##'                ),
##'              units=c(
##'                wt="lb/1000"
##'                ,disp="in^3"
##'                ,qsec="sec"
##'                ),
##'              drop='vs'
##'              )
##' s <- summary(cbind(mpg, disp, hp, drat, wt) ~ cyl + gear + am,
##'              method='cross', overall=TRUE, data=df)
##' dd <- datadist(df)
##' radarplot(S ~ cyl | gear*am, data=s, datadist=dd, rescale="range")
##'
##' ## TODO: Provide example of convenient usage with 'lrm' and 'cph' models.
##' @author David C. Norris
##' @export radarplot
radarplot <- function(x, data, datadist=getOption('datadist'),
                      xlim=c(-1.4, 1.4), ylim=c(-1.2, 1.2), # TODO: automate optimal positioning
                      rescale=c("IQR","range"), # TODO: 1. Permit variable-specific rescaling
                      treatment=NULL,           #       2. Warn if IQR results in negative plots
                      stratify=NULL,
                      strength=NULL,
                      include.na=FALSE, # TODO: Would it *ever* make sense to set this to TRUE? If not, omit!
                      overall=FALSE,
                      ...){
  ## Note that the 'data' parameter may not be required by the x=fit cases below.
  if(is(x, "lrm") || is(x, "cph")) { # Case where strength-of-association model is logistic or Cox PH
    lhs <- attr(x$terms,'variables')
    ## unwrap any rcs-splined variables:
    lhs[[1]] <- as.name('cbind') # convert list to cbind expression
    lhs <- lapply(lhs, function(term) if(length(term)>1 && term[[1]]==as.name('rcs')) term[[2]] else term)
    lhs <- as.call(lhs) # restore to call form after lapply
    lhs[[2]] <- NULL # remove the regression's LHS var
    lhs <- lhs[!as.character(lhs) %in% as.character(substitute(treatment))]
    lhs <- lhs[!as.character(lhs) %in% as.character(substitute(stratify))]
    treatment <- substitute(treatment)
    if(!is.name(treatment))
      treatment <- as.name(treatment)
    stratify <- substitute(stratify)
    ## Select and sort the odds ratios for 'spoke' covariates only
    treatment <- substitute(treatment)
    if(!is.name(treatment))
      treatment <- as.name(treatment)
    stratify <- substitute(stratify)
    if(missing(strength)){ # Unless it was explicitly specified,
      ## we obtain the strengths-of-association as abs(log(odds or hazard ratios))
      strength <- matrix(summary(x)[,'Effect'], nrow=2)[2,]
      names(strength) <- matrix(attr(summary(x),'dimnames')[[1]], nrow=2)[1,]
      ## Eliminate the trailing contrast specifiers to recover the spoke variable names
      names(strength) <- sub(" .*$", "", names(strength))
      ## Select and sort the odds/hazard ratios for 'spoke' covariates only
      strength <- strength[as.character(lhs)[-1]]
      ## Having thus aligned the strengths vector with the spokes,
      ## we now exploit names(strength) as a place to store labels
      ## to appear at the spoke endpoints.
      names(strength) <- format(round(strength, digits=2), digits=2)
      strength <- abs(log(strength)) # sector widths are thus made proportional to abs(log(OR)) or abs(log(HR))
    }
    if(is.null(stratify)){
      env <- list(treatment=as.name(treatment))
      rhs <- substitute(treatment, env)
      x <- as.formula(substitute(S ~ treatment, env))
    } else if (length(stratify)==1){
      env <- list(treatment=as.name(treatment),
                  stratify=as.name(stratify))
      rhs <- substitute(treatment + stratify, env)
      x <- as.formula(substitute(S ~ treatment | stratify, env))
    } else {
      env <- list(treatment=as.name(treatment),
                             stratify.1=as.name(stratify[[2]]),
                             stratify.2=as.name(stratify[[3]]))
      rhs <- substitute(treatment + stratify.1 + stratify.2, env)
      x <- as.formula(substitute(S ~ treatment | stratify.1 * stratify.2, env))
    }
    ##datadist <- datadist(data)
    data <- do.call("summary",
                    list(formula=substitute(lhs ~ rhs),
                         method='cross',
                         include.na=include.na,
                         overall=overall,
                         data=data))
  } # -- end (lrm|cph) case
  else if(is(x, "ols")) { # Case where strength-of-association model is ols
    stop("Radarplot case for 'ols' model not yet implemented") # TODO
  } # -- end 'ols' case
  else { # TODO: Should this branch be folded into the x=fit cases, then disappear?
    ## TODO: If strength is provided, check that it has correct length
    if(!is(data, "data.frame"))
      stop("The 'data' parameter of 'radarplot' must be a data.frame or summary.formula.cross object")
    else if(!is(data, "summary.formula.cross")) { # Replace 'data' with a canonical summary.formula.cross object
      ## RHS of x must be of the form u [| v [* w]]
      if(!is(x,"formula"))
        stop("Argument 'x' to radarplot must be a formula")
      if(length(x[[3]])==1)
        x2 <- x
      else if(length(x[[3]])==3 && as.character(x[[3]][[1]])=="|") {
        if(length(x[[3]][[3]])==1)
          x2 <- substitute(lhs ~ u + v, list(lhs=x[[2]], u=x[[3]][[2]], v=x[[3]][[3]]))
        else if(length(x[[3]][[3]]==3 && as.character(x[[3]][[3]][[1]])=="*"))
          x2 <- substitute(lhs ~ u + v + w, list(lhs=x[[2]], u=x[[3]][[2]], v=x[[3]][[3]][[2]], w=x[[3]][[3]][[3]]))
        else
          stop("The RHS of formula 'x' must be a single variable, or else of the form u | v [* w]")
      }
      data <- eval(substitute(summary(x2, data, method='cross', overall=overall),
                              list(x2=x2)))
      data$S <- t(t(data$S) * (1/sapply(as.data.frame(data$S), max, na.rm=TRUE))) # normalize
      x[[2]] <- as.name("S")
    }
  }
  ## Let the common pathway of all signatures be accomplished via
  ## the is.formula(x) && is(data, "summary.formula.cross") case.
  ## At this point, x should be a formula like S ~ trt (| str (* str2)),
  ## and 'data' should be a summary.formula.cross object.
  ## Furthermore, 'datadist' should be set appropriately.
  if(is.null(datadist)){
      ## TODO: Does the error message below require correction?
    stop("radarplot requires datadist when data is a summary.formula.cross")
  }
  if(is.character(datadist))
    datadist <- get(datadist)
  ## There are three kinds of spoke covariates, each with a natural mapping to the spokes.
  ## For logical variables, the derived proportions should be plotted on the spokes,
  ## mapped linearly to the interval [0, 1].  This requires NO rescaling at all!
  ## For the continuous columns of data$S, we map the spokes to inter-quartile ranges
  ## by default, or else to the ranges themselves.
  ## For bivalued factors, we support the natural interpretation as labelled logicals,
  ## such that e.g., diabetes='Not Diabetic'|'Diabetic' may be treated as if it were
  ## the logical variable diabetic=F/T.
  cvars <- lapply(datadist$limits, is.numeric) # so data$S[,cvars] are the continuous cols
  cvars <- names(cvars)[cvars==TRUE] # now cvars is simply a list of var names
  cvars <- intersect(cvars, colnames(data$S)) # restrict cvars to columns of summary
  cvars <- setdiff(cvars, names(datadist$values))
  if(length(cvars)){
    if(rescale[1]=="IQR"){
      Q1 <- as.numeric(datadist$limits['Low:effect',cvars])
      Q3 <- as.numeric(datadist$limits['High:effect',cvars])
      IQR <- Q3 - Q1
      data$S[,cvars] <- t((t(data$S[,cvars]) - Q1)/IQR)
    } else { # rescale to max-min
      Q0 <- as.numeric(datadist$limits['Low',cvars])
      Q4 <- as.numeric(datadist$limits['High',cvars])
      RANGE <- Q4 - Q0
      data$S[,cvars] <- t((t(data$S[,cvars]) - Q0)/RANGE)
    }
  }
  dvars <- sapply(datadist$values, length) # vector of lengths of discrete variables
  dvars <- dvars[dvars>1] # dvars is now a list of names of discrete variables
  dvars <- dvars[names(dvars) %in% colnames(data$S)] # restrict dvars to columns of summary
  factors <- sapply(datadist$limits, is.factor)
  factors <- names(factors)[factors]
  dvars <- dvars[names(dvars) %in% factors] # further restrict to factors..
  ## TODO: The treatment of factors with 3 or more values as quasi-continuous values
  ##       is conceivably valid only for ordered factors. Unfortunately, I find that
  ##       lrm estimated on a model with ordered RHS factor results in a fit to which
  ##       'summary.rms' cannot be applied. (Generates a 'non-conformable arguments'
  ##       error.) Thus, I cannot feasibly apply this sensible restriction at this time.
  ##  ordfacs <- sapply(datadist$limits, is.ordered)
  ##  ordfacs <- names(ordfacs)[ordfacs]
  ##  dvars <- dvars[dvars==2 | names(dvars) %in% ordfacs] # ..and to ordered factors if >2 levels
  if(length(dvars)){
    data$S[,names(dvars)] <- t(t(data$S[,names(dvars)]) - 1)/dvars
  }

  group <- x[[3]]
  if(length(group)>1)
    group <- group[[2]]
  ## Remove the marginal level for the innermost (group) level
  groupcol <- match(as.character(group), names(data))
  data <- data[data[[groupcol]]!="ALL",]
  ## Convert the grouping column to a factor
  data[,groupcol] <- factor(data[,groupcol])
  ## Eliminate rows for unrepresented factor combinations
  data <- data[data$N>0,]

  ## If any of the radii is negative, warn appropriately:
  if(any(data$S<0)){
      stopifnot(rescale[1]!="range")
      neg.rads <- sapply(as.data.frame(data$S), min, na.rm=TRUE) < 0
      neg.rads <- names(neg.rads)[neg.rads]
      warning(paste("Plotting negative radii for variables: ",
                    paste(neg.rads, collapse=", "),
                    "; consider setting rescale='range'.", sep=""))
  }
  ## Exploit xyplot's standard lattice setup
  eval(substitute(xyplot(x, data=data, groups=group,
                         panel=panel.radarplot, aspect="iso",
                         xlim=xlim, xlab="",
                         ylim=ylim, ylab="",
                         lty=1, lwd=2,
                         scales=list(draw=FALSE),
                         auto.key=list(
                           columns=length(levels(data$group))
                           ,lines=TRUE
                           ,points=FALSE
                           ),
                         radii=colnames(data$S),
                         strength=strength,
                         ...),
                  list(data=data, group=group)))
}
