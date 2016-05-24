"yates.effects" <- function(aov.obj, error.term="Within", data=NULL)
{
	currentOptions <- getOption("contrasts")
	options(contrasts=c("contr.helmert", "contr.helmert"))
  multistratum <- inherits(aov.obj, what="aovlist")
	if (multistratum)
	{	Yates.aov <- aov(as.formula(deparse(attr(aov.obj, "terms"))), data=data)
	  if(any(tapply(rep(1, length(Yates.aov[[error.term]][["assign"]])), 
                  Yates.aov[[error.term]][["assign"]], sum) != 1))
       stop("terms are not all single degree of freedom") 
#    neffects <- length(Yates.aov[[error.term]][["coefficients"]])
#	  nterms <-  length(attr(Yates.aov[[error.term]][["terms"]], "term.labels"))
		int <- match("(Intercept)", names(Yates.aov[[error.term]][["coefficients"]]))
    if (!is.na(int))
       estimable[int] <- FALSE
    estimable <- !is.na(Yates.aov[[error.term]][["coefficients"]])
    yeffects <- (Yates.aov[[error.term]][["coefficients"]])[estimable]*2
  	names(yeffects) <- as.character(attr((Yates.aov[[error.term]][["terms"]])
                             [Yates.aov[[error.term]][["assign"]][estimable]], 
                                                          "term.labels"))
	}
	else
	{	Yates.aov <- aov(as.formula(deparse(aov.obj$terms)), data=data)
	  if(any(tapply(rep(1, length(Yates.aov[["assign"]])), 
                  Yates.aov[["assign"]], sum) != 1))
       stop("terms are not all single degree of freedom") 
#    neffects <- length(Yates.aov[["coefficients"]])
#	  nterms <-  length(attr(Yates.aov[["terms"]], "term.labels"))
    estimable <- !is.na(Yates.aov[["coefficients"]])
		int <- match("(Intercept)", names(Yates.aov[["coefficients"]]))
    if (!is.na(int))
       estimable[int] <- FALSE
		yeffects <- Yates.aov[["coefficients"]][estimable]*2
  	names(yeffects) <- as.character(attr((Yates.aov[["terms"]])
                        [Yates.aov[["assign"]][estimable]], "term.labels"))
	}
	(options(contrasts=currentOptions))
	return(yeffects)
}

"qqyeffects" <- function(aov.obj, error.term="Within", data=NULL, pch=16, full=FALSE, ...)
{
#
# check that full is logical
#
	if (!is.logical(full))
	{ stop("full must be a logical.")
	}
#
# get factorial effects and their names; compute probabilities for quantiles and then quantiles
#
#	aov.tables <- model.tables(aov.obj, type="effects")
#	feffects <- as.vector(unlist(aov.tables$tables))
#	names(feffects) <- names(aov.tables$tables)
	feffects <- yates.effects(aov.obj, error.term, data=data)
	quants <- ppoints(feffects)
	if (full == F)
	{
		feffects <- abs(feffects)
		quants <- (0.5+0.5*quants)
	}
	feffects <- sort(feffects)
	quants <- qnorm(quants)
#
# plot
#
	if (full == F)
	{
		xlab="Half-normal quantiles"
	}
	else
	{
		xlab="Normal quantiles"
	}
	labeff <-	identify(qqplot(quants,feffects, ylab="Factorial effects", 
                     xlab=xlab, pch=pch), labels=names(feffects), plot=T, ...)
# next line adds a line through 75th percentile
#	slope<- quantile(feffects[setdiff(1:length(feffects), labeff)], 0.75)/quantile(quants[setdiff(1:length(feffects), labeff)], 0.75)
#
# compute slope for regression line through the origin for unselected points and plot 
#
	feffects.error <- feffects[setdiff(1:length(feffects), labeff)]
	quants.error <- quants[setdiff(1:length(feffects), labeff)]
	slope <- sum(feffects.error*quants.error)/sum(quants.error*quants.error)
	abline(0, slope)
#	The next two lines fit the regression line through the unselected points
#		coeffs <- coef(lm(feffects.error ~ quants.error))
#		abline(coeffs[1], coeffs[2])
#
# print out names of labelled effects (in case illegible on plot)
#
	if (length(labeff) == 0)
	{
		cat("No effects labelled\n")
	}
	else
	{
		cat("Effect(s) labelled:",names(feffects)[labeff],"\n")
	}
# return plotted points invisibly
	invisible(list(x=quants, y=feffects))
}
