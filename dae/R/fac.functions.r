#"as.numfac" <- function(factor) {as.numeric(as.vector(factor))}

"as.numfac" <- function(factor)
{ if (!is.numeric(factor))
  { levs <- levels(factor)
    if (any(is.na(suppressWarnings(as.numeric(levs[!is.na(levs)])))))
      warning("Some levels do not have values that are numeric in nature")
    #see factor help
    factor <- as.numeric(levels(factor))[factor]
  }
  return(factor)
}

"mpone" <- function(factor) {2*as.numeric(factor)-3}

"fac.recode" <- function(factor, newlevels, ...)
#function to form a new factor by changing the levels of factor
{ lev.fac <- levels(factor)
  nlev <- length(lev.fac)
  if (nlev != length(newlevels))
  { stop("Must supply a new level for every level in the supplied factor")}
  new.fac <- as.character(factor)
  old.fac <- new.fac
  for (i in 1:nlev)
  { new.fac[old.fac == lev.fac[i]] <- newlevels[i]}
  new.fac <- factor(new.fac, ...)
  return(new.fac)
}

"fac.combine" <- function(factors, order="standard", combine.levels=FALSE, sep=",", ...)
{
#
# test arguments
#
	if (mode(factors) != "list" | length(factors) == 1) stop("Must supply a list")
	nfac <- length(factors)
	if (!all(sapply(factors, inherits, what="factor"))) 
		stop("All elements of list must be factors or ordereds")
	if (var(sapply(factors, length)) != 0) stop("All factors must be of the same length")
	which.ord <- pmatch(casefold(order), c("standard", "yates"), nomatch="")
	if (which.ord == "")	stop("order must be either standard or yates")
# standard order
	if (which.ord == "1") counter <- nfac:1
# Yates order
	else if (which.ord== "2") counter <- 1:nfac
#
# compute new factor
#		
	radix <- 1
	new.fac <- rep(1, length(factors[[1]])) 
	for (i in counter)
#reassign factor so unused levels removed
	{ f <- factor(factors[[i]])
    nlev <- length(levels(f))
    new.fac <- (as.numeric(f)-1) * radix + new.fac
    if (combine.levels)
      if (i == counter[1])
          radix.lev <- paste(levels(f))
      else
      { if (which.ord == 1)
          radix.lev <- paste(rep(levels(f), each=radix), rep(radix.lev, times=nlev), sep=sep)
        else
          radix.lev <- paste(rep(radix.lev, times=nlev), rep(levels(f), each=radix), sep=sep)
      }
		radix <- radix * nlev
	}
  if (combine.levels)
  { obslevs <- unique(new.fac)
    obslevs <- obslevs[order(obslevs)] 
  	new.fac <- factor(new.fac, labels=radix.lev[obslevs], ...)
  }	
 	else
  	new.fac <- factor(new.fac, ...)
	return(new.fac)
}

"fac.divide" <- function(combined.factor, factor.names, order="standard")
{
#
# test arguments
#
#	if (mode(factor.names) != "list" | length(factor.names) == 1) 
#    stop("Must supply a list of more than one component")
	if (mode(factor.names) != "list") 
    stop("Must supply a list of factor names")
	nfac <- length(factor.names)
	which.ord <- pmatch(casefold(order), c("standard", "yates"), nomatch="")
	if (which.ord == "")	stop("order must be either standard or yates")
# standard order
	if (which.ord == "1") counter <- nfac:1
# Yates order
	else if (which.ord== "2") counter <- 1:nfac
#
# compute new factor
#		
	kombined.factor <- as.numeric(combined.factor)
	for (i in counter)
	{ klev <- length(factor.names[[i]])
    if (is.numeric(factor.names[[i]]))
    { if (klev != 1)
      { lev <- factor.names[[i]]
        val <- (kombined.factor - 1) %% klev + 1
		   factor.names[[i]] <- factor(lev[val])
      }
      else
      {	klev <- factor.names[[i]]
		    factor.names[[i]] <- factor((kombined.factor - 1) %% klev + 1)
	    }
    }
    else
	  {	factor.names[[i]] <- factor((kombined.factor - 1) %% klev + 1, 
                                          labels = factor.names[[i]])
	  }
		kombined.factor <- (kombined.factor - 1) %/% klev + 1
	}
  new.factors <- data.frame(factor.names)
}

"fac.nested" <- function(nesting.fac, levels=NA, labels=NA, ...)
{
	n <- length(nesting.fac)
	reps <- table(nesting.fac)
	levs <- levels(nesting.fac)
	no.lev.between <- length(levs)
	no.lev.within <- max(table(nesting.fac))
	nested.fac <- c(rep(1, n))
  for(i in 1:no.lev.between)
   nested.fac[nesting.fac == levs[i]] <- 1:reps[i]
  	if (length(levels) == 1 && is.na(levels)) levels <- 1:no.lev.within
	if (length(labels) == 1 && is.na(labels)) labels <- as.character(levels)
	nested.fac <- factor(nested.fac, levels=levels, labels=labels, ...)
	return(nested.fac)
}