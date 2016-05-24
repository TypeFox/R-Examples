	update.sitar <- function (object, ..., evaluate = TRUE)
{
	mcall <- object$call.sitar
	if (is.null(mcall))
		stop("need an object with call.sitar component")
	extras <- as.list(match.call(expand.dots = FALSE)$...)
#	drop null args
	mcall$start <- NULL
	for (i in names(extras))
		if (is.null(extras[[i]]))
			mcall[[i]] <- extras[[i]] <- NULL
#	expand formulae
	if (any(grep('formula', names(extras)))) {
    for (a in letters[1:3]) {
      n <- paste(a, 'formula', sep='.')
      if (!is.null(extras[[n]]) && !is.null(mcall[[n]]))
        extras[[n]] <-  update.formula(mcall[[n]], extras[[n]])
    }
  }
# update args
	mcall[names(extras)] <- extras
#	add start arg if none of these args specified
	if (!sum(pmatch(names(extras), c("x", "y", "id", "fixed", "random", "a.formula", "b.formula",
	                                "c.formula", "start", "returndata"), 0))) {
  	start. <- list(fixed=fixef(object), random=ranef(object))
# update start if any of these args specified
  	if (sum(pmatch(names(extras), c('data', 'subset', 'df', 'knots', 'bounds', 'xoffset', 'bstart'), 0))) {
# get data etc
  		data <- eval(mcall$data)
  		subset <- eval(mcall$subset, data)
  		if (!is.null(subset)) data <- data[subset, ]
  		x <- eval(mcall$x, data)
  		xoffset <- object$xoffset
  		if (is.null(xoffset)) xoffset <- mean(x)
  		x <- x - xoffset
  		df <- object$ns$rank - 1
  		knots <- attr(object$ns$model$ns, 'knots')
  		bounds <- attr(object$ns$model$ns, 'Boundary.knots')
# update random effects
  		if (!is.null(extras$data) || !is.null(extras$subset)) {
  		  id <- factor(eval(mcall$id, data))
  		  levels.obj <- levels(getGroups(object))
  		  if (!identical(levels(id), levels.obj)) {
#	omit random effects for missing levels in id
  		    start.$random <- start.$random[idcheck <- levels.obj %in% levels(id), ]
  		    cat(length(levels.obj) - sum(idcheck), 'subjects omitted\n')
#	add zero random effects for new levels in id
  		    newid <- !levels(id) %in% levels.obj
  		    if (sum(newid) > 0) {
  		      newre <- matrix(0, nrow=sum(newid), ncol=dim(ranef(object))[2],
  		                      dimnames=list(levels(id)[newid], dimnames(ranef(object))[[2]]))
  		      start.$random <- rbind(start.$random, newre)
  		      cat(sum(newid), 'subjects added\n')
  		    }
  		  }
  		}
#	update fixed effects
  		if (length(fixef(object)) > df + 1) fixed.extra <- (df+2):length(fixef(object))
  			else fixed.extra <- NULL
# new arg xoffset
  		if (!is.null(extras$xoffset)) {
  		  xoffset.t <- xoffset
  		  xoffset <- extras$xoffset
  		  xoffset.t <- xoffset - xoffset.t
  		  x <- x - xoffset.t
  		  knots <- knots - xoffset.t
  		  bounds <- bounds - xoffset.t
  		}
# new arg knots
  		if (!is.null(extras$knots)) {
  			knots <- eval(extras$knots) - xoffset
  			df <- length(knots) + 1
  			mcall$df <- NULL
  		}
# new arg df
  		else if (!is.null(extras$df)) {
  			df <- extras$df
  			knots <- quantile(x, (1:(df-1))/df)
  			mcall$knots <- NULL
  		}
# new arg bounds
  		if (!is.null(extras$bounds)) {
  			bounds <- eval(extras$bounds)
  			if (length(bounds) == 1) bounds <- range(x) + abs(bounds) * c(-1,1) * diff(range(x))
  			  else bounds <- bounds - xoffset
  		}
#	get spline start values
  		spline.lm <- lm(predict(object, data, level=0) ~ ns(x, knots=knots, Bound=bounds))
  		start.$fixed <- c(coef(spline.lm)[c(2:(df+1), 1)], start.$fixed[fixed.extra])
# new arg bstart
  		if (!is.null(extras$bstart) && !is.null(start.$fixed['b'])) {
  			bstart <- eval(extras$bstart)
  			if (bstart == 'mean') bstart <- mean(x)
  			  else bstart <- bstart - xoffset
  			start.$fixed['b'] <- bstart
  		}
  	}
#	save start. object
		assign('start.', start., parent.frame())
		mcall[['start']] <- quote(start.)
	}
	if (evaluate)
		eval(mcall, parent.frame())
	else mcall
}
