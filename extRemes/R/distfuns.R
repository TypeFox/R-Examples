revd <- function(n, loc=0, scale=1, shape=0, threshold=0, type=c("GEV","GP")) {
   type <- match.arg(type)
   type <- tolower(type)
   if(type == "gev") z <- rexp(n)
   else if(type=="gp") {
	z <- runif(n)
	loc <- threshold
   }
   else stop("revd: invalid type argument.")
   out <- numeric(n)+NA
   if(min(scale) < 0) stop("revd: scale parameter(s) must be positively valued.")
   if(length(loc)==1) loc <- rep(loc, n)
   if(length(scale)==1) sc <- rep(scale, n)
   else sc <- scale
   if(length(shape)==1) sh <- rep(shape, n)
   else sh <- shape
   if(length(loc) != n || length(sc) != n || length(sh) != n) stop("revd: parameters must have length equal to 1 or n.")
   id <- sh == 0
   if(any(id)) {
	if(type == "gev") out[id] <- loc[id] - sc[id]*log(z[id])
	else if(type=="gp") out[id] <- loc[id] + rexp(n, rate=1/sc[id])
   }
   if(any(!id)) out[!id] <- loc[!id] + sc[!id]*(z[!id]^(-sh[!id]) - 1)/sh[!id]
   return(out)
} # end of 'revd' function.

devd <- function(x, loc=0, scale=1, shape=0, threshold=0, log=FALSE, type=c("GEV","GP")) {
   type <- match.arg(type)
   type <- tolower(type)
   if(any(scale < 0)) stop("devd: invalid scale parameter(s).  Must be > 0.")
   n <- length(x)

   if(missing(loc)) loc <- 0
   else if(is.null(loc)) loc <- 0

   if(length(loc)==1) loc <- rep(loc,n)

   if(length(scale)==1) sc <- rep(scale,n)
   else sc <- scale
   if(length(shape)==1) sh <- rep(shape,n)
   else sh <- shape
   if(length(threshold)==1) th <- rep(threshold,n)
   else th <- threshold

   d <- (x - loc)/sc
   r <- log(1/sc)
   id <- sh == 0
   if(type=="gev") {
	if(any(id)) d[id] <- r[id] - d[id] - exp(-d[id])
	if(any(!id)) {
	   d[!id] <- 1 + sh[!id]*d[!id]
	   good <- (!id & (d > 0)) | is.na(d)
	   if(any(good)) {
		sc <- sc[good]
		sh <- sh[good]
		r <- r[good]
		z <- d[good]
		d[good] <- r - z^(-1/sh) - (1/sh + 1)*log(z)
	   }
	   if(any(!id & !good)) d[!id & !good] <- -Inf
	}
   } else if(type=="gp") {
	good <- (d > 0 & ((1 + sh*d) > 0)) | is.na(d)
	if(any(id) && any(good)) d[id][good[id]] <- r[id][good[id]] - d[id][good[id]]
	if(any(!id) && any(good)) d[!id][good[!id]] <- r[!id][good[!id]] - (1/sh[!id][good[!id]] + 1)*log(1 + sh[!id][good[!id]]*d[!id][good[!id]])
	if(any(!good)) d[!good] <- -Inf
   } else stop("devd: invalid type argument.")
   if(!log) d <- exp(d)
   return(d)
} # end of 'devd' function.

pevd <- function(q, loc=0, scale=1, shape=0, threshold=0, lambda=1, npy,
		    type=c("GEV", "GP", "PP", "Gumbel", "Frechet", "Weibull", "Exponential", "Beta", "Pareto"),
		    lower.tail=TRUE, log.p = FALSE) {

   type <- match.arg(type)
   type <- tolower(type)

   if(is.element(type, c("exponential","gumbel")) && shape != 0) stop("pevd: invalid type argument for choice of threshold.")
   if(scale <= 0) stop("pevd: invalid scale argument.  Must be > 0.")
   if(length(loc) > 1 || length(scale) > 1 || length(shape) > 1) stop("pevd: invalid parameter arguments.  Each must have length 1 only.")
   if(is.element(type, c("weibull","frechet"))) type <- "gev"
   else if(is.element(type, c("exponential", "beta","pareto"))) type <- "gp"

   if(type == "pp") {

	# scale <- scale - shape * (threshold - loc)
	type <- "gev"

   }

   if(is.element(type, c("gev","gumbel","frechet","weibull"))) {

	q <- (q - loc)/scale

	if(!log.p) {

	    if(shape==0 || type=="gumbel") p <- exp(-exp(-q))
	    else p <- exp(-pmax(1 + shape*q, 0)^(-1/shape))

	} else {

	    if(lower.tail) {

	        if(shape==0 || type=="gumbel") p <- -exp(-q)
                else p <- -pmax(1 + shape * q, 0)^(-1/shape)

	    } else {

		if(shape==0 || type=="gumbel") p <- log(1 - exp(-exp(-q)))
                else p <- log(1 - exp(-pmax(1 + shape*q, 0)^(-1/shape)))

	    }

	}

   } else if(is.element(type, c("gp","exponential","beta","pareto"))) {

	q <- pmax(q - threshold, 0)/scale

	if(!log.p) {

	    if(shape==0 || type=="exponential") p <- 1 - exp(-q)
	    else p <- 1 - pmax(1 + shape*q, 0)^(-1/shape)

	} else {

	    if(!lower.tail) {

	        if(shape==0 || type=="exponential") p <- log(1 - exp(-q))
                else p <- log(1 - pmax(1 + shape*q, 0)^(-1/shape))

	    } else {

		if(shape==0 || type=="exponential") p <- -q
                else p <- (-1 / shape) * log( pmax(1 + shape*q, 0))

	    }

	}

   } else stop("pevd: invalid type argument.")

#  else if(type=="pp") {
#
# 	lam <- 1 - exp(-(1 + shape*(threshold - loc)/scale)^(-1/shape)/npy)
#	scale <- scale + shape*(threshold - loc)
#	z <- pmax(1 + shape*(q - threshold)/scale, 0)
#	p <- 1 - (z^(-1/shape))
#
#   } 

   if(!lower.tail && !log.p) p <- 1 - p
   names(p) <- NULL
   return(p)

} # end of 'pevd' function. 
 
qevd <- function(p, loc=0, scale=1, shape=0, threshold=0,
	    type=c("GEV", "GP", "PP", "Gumbel", "Frechet", "Weibull", "Exponential", "Beta", "Pareto"),
	    lower.tail=TRUE) {

   type <- match.arg(type)
   type <- tolower(type)

   if(scale <= 0) stop("qevd: invalid scale argument.  Must be > 0.")

   if(min(p, na.rm=TRUE) <= 0 || max(p, na.rm=TRUE) >= 1) stop("qevd: invalid p argument.  Must have 0 < p < 1.")

   if(length(loc) > 1 || length(scale) > 1 || length(shape) > 1) stop("qevd: invalid parameter arguments.  Each must have length 1 only.")

   if(type=="pp") scale <- scale + shape * (threshold - loc)

   if(is.element(type, c("gev","gumbel","frechet","weibull"))) {

	if(!lower.tail) p <- 1 - p

        if(shape == 0 || type == "gumbel") q <- loc - scale*log(-log(p))
        else q <- loc + scale*((-log(p))^(-shape) - 1)/shape

   } else if(is.element(type, c("pp", "gp","beta","exponential","pareto"))) {

	if(lower.tail) p <- 1 - p

	if(shape == 0 || type == "exponential") q <- threshold - scale * log(p)
	else q <- threshold + scale * (p^(-shape) - 1)/shape

   } else stop("qevd: invalid type argument.")

   return(q)

} # end of 'qevd' function.

rlevd <- function(period, loc=0, scale=1, shape=0, threshold=0,
	    type=c("GEV", "GP", "PP", "Gumbel", "Frechet", "Weibull", "Exponential", "Beta", "Pareto"),
	    npy=365.25, rate=0.01) {

    if(any(period <= 1)) stop("rlevd: invalid period argument.  Must be greater than 1.")

    type <- match.arg(type)
    type <- tolower(type)

    if(missing(loc)) loc <- 0
    else if(is.null(loc)) loc <- 0

    if(is.element(type, c("gumbel","weibull","frechet"))) {
        if(type=="gumbel" && shape != 0) {
    	    warning("rlevd: shape is not zero, but type is Gumbel.  Re-setting shape parameter to zero.")
    	    shape <- 0
    	    type <- "gev"
        } else if(type=="gumbel") type <- "gev" 
        else if(type=="frechet" && shape <= 0) {
	    if(shape==0) {
	        warning("rlevd: shape is zero, but type is Frechet!  Re-setting type to Gumbel.")
	        shape <- 0
	    } else {
	        warning("rlevd: type is Frechet, but shape < 0.  Negating shape to force it to be Frechet.")
	        shape <- -shape
	    }
            type <- "gev"
        } else if(type=="frechet") type <- "gev" 
        else if(type=="weibull" && shape >= 0) {
	    if(shape==0) {
	        warning("rlevd: shape is zero, but type is Weibull!  Re-setting type to Gumbel.")
	        shape <- 0
	    } else {
                warning("rlevd: type is Weibull, but shape > 0.  Negating shape to force it to be Weibull.")
                shape <- -shape
            }
            type <- "gev"
        } else if(type=="weibull") type <- "gev"
    } # end of if type is from the GEV family but not GEV stmts.

    if(is.element(type, c("beta","pareto","exponential"))) {
        if(type=="exponential" && shape != 0) {
	    warning("rlevd: shape is not zero, but type is Exponential.  Re-setting shape parameter to zero.")
            shape <- 0
            type <- "gp"
        } else if(type=="exponential") type <- "gp" 
        else if(type=="beta" && shape >= 0) {
	    if(shape==0) {
	        warning("rlevd: shape is zero, but type is Beta!  Re-setting type to Exponential.")
	        shape <- 0
            } else {
                warning("rlevd: type is Beta, but shape > 0.  Negating shape to force it to be Beta.")
                shape <- -shape
            }
            type <- "gp"
        } else if(type=="beta") type <- "gp" 
        else if(type=="pareto" && shape <= 0) {
            if(shape==0) {
	        warning("rlevd: shape is zero, but type is Pareto!  Re-setting type to Exponential.")
	        shape <- 0
            } else {
                warning("rlevd: type is Pareto, but shape < 0.  Negating shape to force it to be Pareto.")
                shape <- -shape
            }
            type <- "gp"
        } else if(type=="pareto") type <- "gp"
    } # end of if type is from GP family, but not GP stmts.

    if(is.element(type, c("gev", "pp"))) {
	p <- 1 - 1/period
	res <- qevd(p=p, loc=loc, scale=scale, shape=shape, type="GEV")
    } else if(type=="gp") {
	m <- period * npy * rate
	if(shape==0) res <- threshold + scale * log(m)
	else res <- threshold + (scale/shape) * (m^shape - 1)
    }
    names(res) <- as.character(period)
    return(res)
} # end of 'rlevd' function.

