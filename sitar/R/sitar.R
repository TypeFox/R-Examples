	sitar <- function(x, y, id, data, df, knots, fixed=random, random='a+b+c',
	                  a.formula=~1, b.formula=~1, c.formula=~1, bounds=0.04, start, bstart=xoffset, xoffset='mean',
	                  returndata=FALSE, verbose=FALSE, correlation=NULL, weights=NULL, subset=NULL, method='ML',
	                  na.action=na.fail, control = nlmeControl(returnObject=TRUE))
#
#	fit growth curves of y ~ ns(f(x)) by id
#	df - number of knots (or length of knots + 1)
#	knots - default x quantiles
#	fixed - default random
#	random - random effects, default "a+b+c"
#	a.formula etc -  fixed effect formulae, default ~ 1
#			to omit fixed effects use ~ -1
#	bounds - span of ns, default x-range 4%
#	start - starting values - default estimated
#			requires spline coefficients, any missing zeroes added
#	bstart - starting value for b, default xoffset
#	xoffset - offset for x, default 'mean', alternatives 'apv' or value
# returndata - if TRUE returns nlme data frame, not nlme model
#	verbose etc - arguments passed to nlme

#	returns call, ns, nlme.out if returndata FALSE
#			else fulldata
{
	b.origin <- function(b) {
		if (b == 'mean') return(mean(x))
		if (b == 'apv') {
			spline.lm <- lm(y ~ ns(x, knots=knots, Bound=bounds))
			return(makess(x, fitted(spline.lm))$apv[1])
		}
		if (!is.numeric(b)) stop('must be "mean" or "apv" or numeric')
		b
	}

# get data
	mcall <- match.call()
	data <- eval(mcall$data, parent.frame())
	x <- eval(mcall$x, data)
	y <- eval(mcall$y, data)
	xoffset <- b.origin(xoffset)
	x <- x - xoffset

# get df, knots and bounds
	if (missing(df) & missing(knots)) stop("either df or knots must be specified")
	if (!missing(df) & !missing(knots)) cat("both df and knots specified - df redefined from knots\n")
	if (missing(knots)) knots <- quantile(x, (1:(df-1))/df) else {
	  knots <- knots - xoffset
	  df <- length(knots) + 1
	}
	if (nrow(data) <= df) stop("too few data to fit spline curve")
	if (length(bounds) == 1) bounds <- range(x) + abs(bounds) * c(-1,1) * diff(range(x))
	else if (length(bounds) == 2) bounds <- bounds - xoffset
	else stop("bounds should be length 1 or 2")

#	get spline start values
	spline.lm <- lm(y ~ ns(x, knots=knots, Bound=bounds))

#	get starting values for ss, a and b
	if (missing(start)) start <- coef(spline.lm)[c(2:(df+1), 1)]

#	force fixed effect for a
	fix <- fixed
	if (!grepl('a', fix)) fix <- paste('a', fix, sep='+')

# set up args for fitnlme
	fixed <- ss <- paste('s', 1:df, sep='')
	pars <- c('x', ss)

#	set up model elements for a, b and c
	names(model) <- model <- letters[1:3]
	constant <- mm.formula <- as.formula('~ 1')
	# mm.intercept <- TRUE
	id <- eval(mcall$id, data)
	subset <- eval(mcall$subset, data)
	if (is.null(subset)) subset <- 1:length(x)
	fulldata <- data.frame(x, y, id, subset)
	for (l in model) {
		if (!grepl(l, fix) && !grepl(l, random)) {
			model[l] <- NA
			next
		}
		pars <- c(pars, l)
		formula <- get(paste(l, 'formula', sep='.'))
		if (formula == as.formula('~ -1') || formula == as.formula('~ 1-1') || !grepl(l, fix)) next
		if (formula == constant) mm.intercept <- TRUE
		else {
			if (formula != mm.formula) {
				mm <- model.matrix(formula, data)
				if (dim(mm)[[1]] < length(y))
					stop('Missing values in data')
				mm.formula <- formula
				#	convert spaces to underlines
				colnames(mm) <- gsub(' ', '_', colnames(mm))
				#	convert colons to dots (from interaction)
				colnames(mm) <- gsub(':', '.', colnames(mm), fixed=TRUE)
				mm.intercept <- colnames(mm)[1] == '(Intercept)'
				# omit constant columns and centre others
				cmm <- apply(mm, 2, function(x) sd(x, na.rm=TRUE) > 0)
				mm.names <- colnames(mm)[cmm]
				mm <- as.matrix(mm[, mm.names])
				colnames(mm) <- mm.names
				mm <- scale(mm, scale=FALSE)
			}
			if (exists('mm')) for (i in 1:dim(mm)[[2]]) {
				var <- colnames(mm)[i]
				rc <- paste(l, var, sep='.')
				pars <- c(pars, rc)
				fixed <- c(fixed, rc)
				if (!is.list(start)) start <- c(start, 0)
				model[l] <- paste(model[l], '+', rc, '*', var, sep='')
				if (!var %in% pars) {
					pars <- c(pars, var)
					fulldata <- cbind(fulldata, mm[,i])
					names(fulldata)[dim(fulldata)[2]] <- var
					# assign(var, mm[,i])
				}
			}
		}
		if (mm.intercept) {
			fixed <- c(fixed, l)
			if (!is.list(start) && l != 'a') {
			  if (l == 'b') start <- c(start, b.origin(bstart) - xoffset)
			    else start <- c(start, 0)
			}
		}
	}

# 	if (!is.null(weights)) {
#     if (is.list(weights)) form <- asOneFormula(lapply(weights, function(z) attr(z, 'formula')))
#       else form <- attr(weights, 'formula')
#     if (!is.null(form)) {
#       wt <- as.data.frame(model.frame(form, data))
#       wt.names <- names(wt)[!names(wt) %in% names(fulldata)]
#       wt <- as.data.frame(wt[, wt.names])
#       names(wt) <- wt.names
#       fulldata <- cbind(fulldata, wt)
# 	  }
# 	}

	if (returndata) invisible(fulldata) else {
		pars <- paste(pars, collapse=',')
		fixed <- paste(fixed, collapse='+')
		sscomma <- paste(ss, collapse=',')

#	combine model elements
		nsd <- paste(model['a'], '+')
		nsf <- paste('(x', ifelse(!is.na(model['b']), paste('- (', model['b'], '))'), ')'))
		if (!is.na(model['c'])) nsf <- paste(nsf, '* exp(', model['c'], ')')

#	code to parse
  	fitcode <- c(
	"fitenv <- new.env()",
	"fitenv$fitnlme <- function($pars) {",
	"as.vector( $nsd",
	"(as.matrix(cbind($sscomma)) * as.matrix(ns($nsf,",
	"knots=knots, Boundary.knots=bounds))) %*%",
	"matrix(rep(1,df), ncol=1))",
	"}",
	"on.exit(detach(fitenv))",
	"attach(fitenv)",
	"nlme(y ~ fitnlme($pars),",
	"fixed = $fixed ~ 1,",
	"random = $random ~ 1 | id,",
	"data = fulldata,",
	"start = start, correlation = correlation,",
	"weights = weights, subset = subset, method = method,",
	"na.action = na.action, control = control, verbose = verbose)")

		for (i in c('random', 'pars', 'fixed', 'sscomma', 'nsd', 'nsf'))
  		fitcode <- gsub(paste('$', i, sep=''), get(i), fitcode, fixed=TRUE)

#	print values
		if (verbose) {
			cat('\nconstructed code', fitcode, sep='\n')
			cat('\ndf', df, 'bstart', bstart, 'xoffset', xoffset, '\nknots\n', knots, '\nbounds\n', bounds)
			if (is.list(start)) {
				cat('\nstarting values\n  fixed effects\n', start$fixed)
				if (!is.null(start$random)) {
					cat('\n  random effects\n')
					print(start$random)
				}
			}
			else cat('\nstarting values\n', start, '\n')
		}

#	save fitted model
    nlme.out <- eval(parse(text=fitcode))
#     if (exists('start.')) rm(start., inherits=TRUE)
    nlme.out$fitnlme <- fitenv$fitnlme
		nlme.out$call.sitar <- mcall
		nlme.out$xoffset <- xoffset
		nlme.out$ns <- spline.lm
		if (!'sitar' %in% class(nlme.out)) class(nlme.out) <- c('sitar', class(nlme.out))
		nlme.out
	}
}
