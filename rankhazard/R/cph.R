rankhazardplot.cph <- function (
  cphobj, data, select = NULL, refpoints = NULL, 
  CI_level = 0.95, x_CI = NULL, draw.confint = FALSE, 
  legendtext = NULL, axistext = NULL, legendlocation = "top",
  axistextposition = -0.1, reftick = TRUE, refline = FALSE, 
  col.refline = 1, lwd.refline = 1, lty.refline = 2, 
  ylab = NULL, ylim = NULL, yticks = NULL, yvalues = NULL, 
  xtext =TRUE, plottype = "hazard", axes = TRUE, na.rm = TRUE,
  draw = TRUE, return = FALSE, col = NULL, lwd = 1, lty = 1, 
  pch = NULL, cex = 1, bg = "transparent", pt.lwd = 1, 
  col.CI = col, lty.CI = lty +1, lwd.CI = lwd, add = FALSE, 
  graphsbefore = 0, args.legend = NULL, ...)					
{
  if (!is.data.frame(data)) 
    stop("Covariate data must be provided as a data frame by the argument 'data'.")

  if (is.null(cphobj$x) && is.null(x_CI))
    stop("To calculate confidence intevals covariate data need to be provided either as argument 'x_CI' or as 'cphobj$x'.")

  if (is.null(x_CI)) x_CI <- as.data.frame(cphobj$x)
  if (!is.data.frame(x_CI)) stop("'x_CI' must be a data frame.")
 
  term_labels <- cphobj$Design$name 

  if (is.null(select))
    select <- 1:length(term_labels)
  if (max(select) > length(term_labels))
    stop("There are fewer covariates in the model than selected.")
    
  labelcheck <- !is.element(term_labels, names(data))
  if (any(labelcheck)) stop (paste("Covariate(s) not in the argument 'data':", paste(paste(term_labels[labelcheck], collapse = ", "), ".", sep ="")))
    
  data_labels <- term_labels						
  x <- data[data_labels]

  if (na.rm) x <- na.omit(x)

  factorlevels <- cphobj$Design$parms
  factorlabs <- names(factorlevels)
  factors <- which(is.element(term_labels, factorlabs))
  nonfactors <- which(!is.element(term_labels, factorlabs))

  refs <- data.frame(setNames(replicate(length(data_labels), numeric(0), simplify = F), data_labels))
  # the defaults for the reference points are calculated for all covariates
  for(i in nonfactors)
  	refs[1, i] <- median(x[, i], na.rm = TRUE)
  j <- 1
  for(i in factors){
    refs[1, i] <- factorlevels[[j]][1]
    refs[i] <- factor(refs[i], levels = factorlevels[[j]])
    j <- j + 1
  }

  if(!is.null(refpoints)){
    
    if(length(refpoints) != length(select)) 
      stop("The length of 'refpoints' must be the same as the number of covariates to be plotted.")
    
    change <- which(!is.na(refpoints))
    
    if (any(is.na(as.numeric(refpoints[is.element(select, intersect(select[change], nonfactors))])))) warning("'refpoints' must be numeric for variables which are not factors.")
    j <- 1
    for (i in factors){
      if (is.element(i, select[change]))
        if (!is.element(refpoints[which(select == i)], factorlevels[[j]])) 
          warning(paste("The given value in 'refpoints' for factor", factorlabs[j], "is not one of the levels of the factor."))
      j <- j + 1
    } 
    # the defaults are replaced by the given refpoints
    if(is.numeric(refpoints)) refs[1, select[change]] <- refpoints[change]
    else{
      refs[1, intersect(factors, select[change])] <- refpoints[is.element(select, intersect(factors, select[change]))]
      refs[1, intersect(nonfactors, select[change])] <- as.numeric(refpoints[is.element(select, intersect(nonfactors, select[change]))])
    }
  }

  predictions <- predict(cphobj, type = "terms", newdata = x)
  pred_refvalues <- predict(cphobj, type = "terms", newdata = refs)

  refvalues <- as.vector(pred_refvalues)
  names(refvalues) <- attr(pred_refvalues, "dimnames")[[2]]

  xp <- as.data.frame(predictions)

  ### Calculating the confidence intervals ###

  refs[factors] <- 0
  refs <- as.vector(as.matrix(refs))
  
  if (CI_level <= 0 || 1 <= CI_level) stop ("'CI_level' must be a number between 0 and 1.")  
  confinterval <- rankhazard_CI.cph(cphobj, x_CI, refs, CI_level)
    
  select_CI <- confinterval$select_CI
  selecttext <- select

  if (draw.confint){
    CI <- confinterval
    select <- which(is.element(select_CI, select))
    selecttext <- select_CI[select]
    if(length(select) == 0)
      stop("Confidence intevals cannot be caluclated for selected covariates.")
  } 
  else CI <- NULL
    
  if (!is.null(legendtext) && is.null(axistext))
    axistext <- legendtext
  if (is.null(legendtext) && !is.null(axistext))
    legendtext <- axistext
  if (is.null(legendtext)){
    legendtext <- attr(cphobj$terms, "term.labels")[selecttext]
    axistext <- term_labels[selecttext]		
  }

  if (draw)
    rankhazardplot.default( 						
      x = x[select], xp = xp[select], refvalues = refvalues[select], 					
      legendtext = legendtext, axistext = axistext,
      na.rm = na.rm, select = select, confinterval = CI, draw.confint = NULL,
      legendlocation = legendlocation, axistextposition = axistextposition, 
      reftick = reftick, refline = refline, col.refline = col.refline, 
      lwd.refline = lwd.refline, lty.refline = lty.refline, ylab = ylab, 
      ylim = ylim, yticks = yticks, yvalues = yvalues, xtext =xtext, plottype = plottype, 
      col = col, lwd = lwd, lty = lty, pch = pch, axes = axes,
      cex = cex, bg = bg, pt.lwd = pt.lwd, 
      col.CI = col.CI, lwd.CI = lwd.CI, lty.CI = lty.CI, add = add, graphsbefore = graphsbefore,
      args.legend = args.legend, ...)

    if (return)
      return(list(x = x, xp = xp, refvalues = refvalues, confinterval = confinterval))
}