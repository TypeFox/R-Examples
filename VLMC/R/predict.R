#### $Id: predict.R,v 1.15 2015/04/16 16:22:33 maechler Exp $
predict.vlmc <-
function(object, newdata,
         type = c("probs", "class","response", "id.node", "depth", "ALL"),
         se.fit = FALSE,
         ## dispersion = NULL,
         ## terms=NULL,
         allow.subset = TRUE, check.alphabet = TRUE,
         ...)
{
    ## Predict probabilities for new series from a fitted "vlmc" object
    ## -----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 10 Apr 2000, 16:36

    ## o type = c("link", "response", "terms"),
    ##		 is just 'aped' from predict.glm ...
    ##   predict.multinom  seems better :
    ##	 	(object, newdata, type = c("class", "probs"), ...)
    ## 	probs := Pr[ X_i | context ]

    if(!is.vlmc(object))
	stop("first argument must be a \"vlmc\" object; see ?vlmc")

    type <- match.arg(type)

    ## FIXME: se.fit & dispersion are NOT YET supported
    if (!missing(se.fit))	.NotYetUsed("se.fit")
    ##if (!missing(dispersion))	.NotYetUsed("dispersion")

    alphabet <- strsplit(object$alpha,NULL)[[1]]

    ivlmc <- object $ vlmc

    ## newdata : must be discrete time series -- just as "dts" in  vlmc()
    if(missing(newdata)) {
        newdata <- object $ y
        if(is.null(newdata))
            stop("no 'newdata' specified, and object$y is NULL")
    }
    else if(is.character(newdata)) {
	if(!all(1 == nchar(newdata)))
	    stop("character argument must have *all* 1-character strings")
    } else if(!(is.factor(newdata) || is.numeric(newdata)))
        stop("'newdata' must be discrete t.s.: character, factor, or numeric (in 0:m1)")

    if(!is.factor(newdata)) # must make sure the integer conversion is ok
        newdata <- factor(newdata, levels = alphabet)
    ## newdata is now a factor
    n <- length(newdata)
    int.data <- as.integer(newdata) - 1L

    if(check.alphabet) {
	nABC <- levels(newdata) # possibly unsorted
	alpha.len <- length(nABC)
	if(alpha.len > object$alpha.len)
	    stop("alphabet of 'newdata' is larger than the vlmc fit 'object' one")
	is.smaller <- alpha.len < object$alpha.len
	if(is.smaller && !allow.subset)
	    stop("alphabet of 'newdata' is smaller than the vlmc fit 'object' one")
	if(any(nchar(nABC) > 1)) {
	    nABC <- abbreviate(nABC, 1)
	    if(any(nchar(nABC) > 1))
		nABC <-
		    if(alpha.len <= 10) as.character(0:(alpha.len - 1))
		    else letters[1:alpha.len]
	}
	Alpha <- paste(nABC, collapse = "")
        if(Alpha != if(is.smaller) substr(object$alpha, 1,alpha.len)
        else object$alpha)
	    warning(paste("alphabet of newdata '", Alpha,
			  "' differs from that of the vlmc fit '",
                          object$alpha,"'", sep=""))
    }
    m <- as.integer(object $ alpha.len)

    ##-- consider allowing MULTIPLE types simultaneously --
    ## "kind" coding for call to .C():
    kind <-
        as.integer(if(type == "ALL") 1 + 4 # Probs + ID => get everything
                   else switch(type,
                               probs = 1,
                               class =, response = 2,
                               id.node = 4,
                               depth = 8)
                   )
    res <- flags <- integer(n); res[] <- NA
    do.probs <- type %in% c("probs","ALL")
    probs <-
        if(do.probs)
            matrix(as.double(NA), n, m,
                   dimnames = list(as.character(newdata), alphabet))
###%-- TODO : first row of prob[], instead of NA, use marginals Pr[ Y[i] = k ]!!
        else double(0)

    ## This gives the prediction Probabilities / Class / Context.Nr / Depth
    r <- .C(predict_vlmc_p,
            vlmc.vec	= as.integer(ivlmc),
            size	= length(ivlmc),
            alpha.len	= m,
            Data	= int.data,
            data.len	= n,
            kind	= kind,

            ## Output (one of these, depending on 'kind'):
            res		= res,
            flags	= flags,
            probs	= probs,

            NAOK	= TRUE
            ## , DUP      = FALSE ## CRAN-forbidden now
            )[c("res", "probs", "flags")]
    if(type == "probs")
        r[["probs"]] ## was  structure(r[["probs"]], flags = r[["flags"]])
    else if(type == "ALL") {
        names(r)[1] <- "ID"
        structure(c(r,
                    list(ctxt = id2ctxt(r $ID, m = m, alpha = object$alpha),
                         fitted = alphabet[max.col(r $probs)],
                         alpha = object$alpha, alpha.len = m)),
                  class = "predict.vlmc")
    }
    else structure(if (type == "class")
                        factor(alphabet[1+r[["res"]]], levels=alphabet)
                   else r[["res"]])# , flags = r[["flags"]]
}

fitted.vlmc <- function(object, ...) predict(object, type = "class", ...)

print.predict.vlmc  <- function(x, digits = getOption("digits"), ...)
{
    if(!inherits(x,"predict.vlmc") ||
       is.null(x$probs) || is.null(x$ID) || is.null(x$ctxt))
        stop("not a valid 'predict.vlmc' object")
    Fprob <- function(x) as.character(fractions(x))
    colnames(x$probs) <- paste("Pr[X=",colnames(x$probs),"]")
    print(noquote(cbind(fit = x$fitted,
                        Fprob(x$probs),
                        id  = x $ ID,
                        flags= x $ flags,
                        ctxt = x $ ctxt)))
    invisible(x)
}


if(FALSE) {
    ## NOTE:  base R
    ##  family.R has  the following deviance for the binomial family :
    dev.resids <- function(y, mu, wt)
	2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) +
		  (1 - y) * log(ifelse(y == 1, 1, (1 - y)/(1 - mu))))
}

## OLD (and identical, but much slower):
deviance.vlmc <- function(object, ...)
{
    ## Purpose: Deviance =  - 2 * log likelihood(.)
    ## Author: Martin Maechler, Date: 26 Apr 2000, 15:30
    dn <- dimnames(pr <- predict(object))
    n <- nrow(pr)
    -2 * sum(log(pr[cbind(1:n, match(dn[[1]],dn[[2]]))[-1,]]))
}
## Fast
deviance.vlmc <- function(object, ...) -2 * entropy(object)

## Author: Martin Maechler, Date: 26 Apr 2000, 15:30

## Using "classwise" (a la multinom()) + the same "type"s as residuals.glm :
residuals.vlmc <-
    function(object, type = c("classwise", "deviance",
                     "pearson", "working", "response", "partial"),
             y = object$y, ...)
{
    type <- match.arg(type) # i.e. default is "classwise"
    if(!is.character(y))
        stop("y argument must be a proper discrete time series (as char.)")
    if(!is.numeric(m <- object$alpha.len))
        stop("no valid $alpha.len in first argument")

    switch(type,
           classwise =,
           deviance = {
               n <- nrow(pr <- predict(object, new = y))# n x m  matrix
               ij <- cbind(1:n, match(y, colnames(pr)))
           },
           pearson =,
           working =,
           response =
           n <- length(mu <- predict(object, new = y, type = "response")),
                                        # in {0:(m-1)}
           ## mu differs from glm(., binomial) for m = 2 {have no link!}

           partial = stop("partial residuals not yet available"))

    if(type == "classwise") {
        Y <- matrix(0, n, m)
        Y[ij] <- 1
        return(Y - pr)
    }
    else { ## type != "classwise"
        ##if(m > 2 && !quiet) ## dubious for m > 2!
        ##warning("vlmc( m > 2 ): mostly only 'classwise' residuals make sense")
        yi <- alpha2int(y, object$alpha)
        switch(type,
               deviance = {
                   ## The absolute deviance residuals for general |alphabet| :
                   dr <- sqrt(-2*log(pr[ij]))
                   mu <- yi[max.col(pr)]# << some randomness here (sign only)!
                   ifelse(yi > mu, dr, -dr)
               },
               pearson =,
               working =,
               response = yi - mu)
    }
}

## This is how  residuals.glm(.)  works  {for dev.resids(), see above }
if(FALSE){
    mu	<- (object$fitted.values)#.Alias
    wts <- (object$prior.weights)#.Alias
    switch(type,
	   deviance = if(object$df.res > 0) {
	       d.res <- sqrt(pmax((object$family$dev.resids)(y, mu, wts), 0))
	       ifelse(y > mu, d.res, -d.res)
	   } else rep(0, length(mu)),
	   pearson = object$residuals * sqrt(object$weights),
	   working = object$residuals,
	   response = y - mu,
	   partial = object$residuals + predict(object,type="terms")
	   )
}
