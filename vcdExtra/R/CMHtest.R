# Cochran-Mantel-Haenszel tests for ordinal factors in contingency tables

# The code below follows Stokes, Davis & Koch, (2000). 
#   "Categorical Data Analysis using the SAS System", 2nd Ed.,
#   pp 74--75, 92--101, 124--129.

# Ref: Landis, R. J., Heyman, E. R., and Koch, G. G. (1978), 
#		Average Partial Association in Three-way Contingency Tables: 
#		A Review and Discussion of Alternative Tests, 
#		International Statistical Review, 46, 237-254. 

# See: https://onlinecourses.science.psu.edu/stat504/book/export/html/90
# http://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_freq_a0000000648.htm

# DONE: this should be the main function, handling 2-way & higher-way tables
#  With strata, use apply() or recursion over strata
# DONE: With strata, calculate overall CMH tests controlling for strata
# FIXED: rmeans and cmeans tests were labeled incorrectly

CMHtest <- function(x, ...)
  UseMethod("CMHtest")

CMHtest.formula <-
function(formula, data = NULL, subset = NULL, na.action = NULL, ...)
{

  m <- match.call(expand.dots = FALSE)
  edata <- eval(m$data, parent.frame())

  fstr <- strsplit(paste(deparse(formula), collapse = ""), "~")
  vars <- strsplit(strsplit(gsub(" ", "", fstr[[1]][2]), "\\|")[[1]], "\\+")
  varnames <- vars[[1]]

  condnames <- if (length(vars) > 1) vars[[2]] else NULL

  dep <- gsub(" ", "", fstr[[1]][1])
  if (!dep %in% c("","Freq")) {
     if (all(varnames == ".")) {
       varnames <- if (is.data.frame(data))
         colnames(data)
       else
         names(dimnames(as.table(data)))
       varnames <- varnames[-which(varnames %in% dep)]
     }

    varnames <- c(varnames, dep)
  }

  if (inherits(edata, "ftable") || inherits(edata, "table") || length(dim(edata)) > 2) {
    condind <- NULL
    dat <- as.table(data)
    if(all(varnames != ".")) {
      ind <- match(varnames, names(dimnames(dat)))
      if (any(is.na(ind)))
        stop(paste("Can't find", paste(varnames[is.na(ind)], collapse=" / "), "in", deparse(substitute(data))))

      if (!is.null(condnames)) {
        condind <- match(condnames, names(dimnames(dat)))
        if (any(is.na(condind)))
          stop(paste("Can't find", paste(condnames[is.na(condind)], collapse=" / "), "in", deparse(substitute(data))))
        ind <- c(condind, ind)
      }
      dat <- margin.table(dat, ind)
    }
    CMHtest.default(dat, 
                   strata = if (is.null(condind)) NULL else match(condnames, names(dimnames(dat))), ...)
  } else {
      m <- m[c(1, match(c("formula", "data", "subset", "na.action"), names(m), 0))]
      m[[1]] <- as.name("xtabs")
      m$formula <-
          formula(paste(if("Freq" %in% colnames(data)) "Freq",
                        "~",
                        paste(c(varnames, condnames), collapse = "+")))
      tab <- eval(m, parent.frame())
      CMHtest.default(tab, ...)
  }
}

CMHtest.default <- function(x, strata = NULL, rscores=1:R, cscores=1:C, 
	types=c("cor", "rmeans", "cmeans", "general"),
	overall=FALSE, details=overall,  ...)
{

	snames <- function(x, strata) {
	  	sn <- dimnames(x)[strata]
	  	dn <- names(sn)
	  	apply(expand.grid(sn), 1, function(x) paste(dn, x, sep=":", collapse = "|"))
	}

  ## check dimensions
  L <- length(d <- dim(x))
  if(any(d < 2L)) stop("All table dimensions must be 2 or greater")
  if(L > 2L & is.null(strata)) strata <- 3L:L
  if(is.character(strata)) strata <- which(names(dimnames(x)) == strata)
  if(L - length(strata) != 2L) stop("All but 2 dimensions must be specified as strata.")

  ## rearrange table to put primary dimensions first
  x <- aperm(x, c(setdiff(1:L, strata), strata))

	d <- dim(x)
	R <- d[1]
	C <- d[2]

  # handle strata
  if (!is.null(strata)) {
    sn <- snames(x, strata)
    res <- c(apply(x, strata, CMHtest2, rscores=rscores, cscores=cscores, 
			types=types,details=details, ...))
    # DONE: fix names if there are 2+ strata
    names(res) <- sn
    for (i in seq_along(res)) res[[i]]$stratum <- sn[i]
    # DONE: Calculate generalized CMH, controlling for strata
		if (overall) {
			if (!details) warning("Overall CMH tests not calculated because details=FALSE")
			else {
				resall <- CMHtest3(res, types=types)
				res$ALL <- resall
				}
			}
		return(res)
	}
	else CMHtest2(x, rscores=rscores, cscores=cscores, 
		types=types,details=details, ...)	
}

# handle two-way case, for a given stratum
#  DONE:  now allow rscores/cscores == 'midrank' for midrank scores
#  DONE:  allow rscores/cscores=NULL for unordered factors, where ordinal
#     scores don't make sense
#  DONE: modified to return all A matrices as a list
#  DONE: cmh() moved outside

CMHtest2 <- function(x, stratum=NULL, rscores=1:R, cscores=1:C, 
	types=c("cor", "rmeans", "cmeans", "general"),
	details=FALSE, ...) {

	# left kronecker product
	lkronecker <- function(x, y, make.dimnames=TRUE, ...)
		kronecker(y, x, make.dimnames=make.dimnames, ...)
		
	# midrank scores (modified ridits) based on row/column totals
	midrank <- function (n) {
		cs <- cumsum(n)
		(2*cs - n +1) / (2*(cs[length(cs)]+1))
	}

    L <- length(d <- dim(x))
	R <- d[1]
	C <- d[2]
	
	if (is.character(rscores) && rscores=="midrank") rscores <- midrank(rowSums(x))
	if (is.character(cscores) && cscores=="midrank") cscores <- midrank(colSums(x))

	nt <- sum(x)
	pr <- rowSums(x) / nt
	pc <- colSums(x) / nt
	
	m <- as.vector(nt * outer(pr,pc))  # expected values under independence
	n <- as.vector(x)                  # cell frequencies
	
	V1 <- (diag(pr) - pr %*% t(pr))
	V2 <- (diag(pc) - pc %*% t(pc))
	V <- (nt^2/(nt-1)) * lkronecker(V1, V2, make.dimnames=TRUE)
	
	if (length(types)==1 && types=="ALL") types <- c("general", "rmeans", "cmeans", "cor" )
	types <- match.arg(types, several.ok=TRUE)
	# handle is.null(rscores) etc here
	if (is.null(rscores)) types <- setdiff(types, c("cmeans", "cor"))
	if (is.null(cscores)) types <- setdiff(types, c("rmeans", "cor"))

	table <- NULL
	Amats <- list()
	if("cor" %in% types) {
		A <- lkronecker( t(rscores), t(cscores) )
		df <- 1
		table <- rbind(table, cmh(n, m, A, V, df))
		Amats$cor <- A
		}
	if("rmeans" %in% types) {
		A <- lkronecker( cbind(diag(R-1), rep(0, R-1)), t(cscores))
		df <- R-1
		table <- rbind(table, cmh(n, m, A, V, df))
		Amats$rmeans <- A
		}
	if("cmeans" %in% types) {
		A <- lkronecker( t(rscores), cbind(diag(C-1), rep(0, C-1)))
		df <- C-1
		table <- rbind(table, cmh(n, m, A, V, df))
		Amats$cmeans <- A
		}
	if ("general" %in% types) {
		A <- lkronecker( cbind(diag(R-1), rep(0, R-1)), cbind(diag(C-1), rep(0, C-1)))
		df <- (R-1)*(C-1)
		table <- rbind(table, cmh(n, m, A, V, df))
		Amats$general <- A
		}

	colnames(table) <- c("Chisq", "Df", "Prob")
	rownames(table) <- types
	xnames <- names(dimnames(x))
	result <- list(table=table, names=xnames, rscores=rscores, cscores=cscores, stratum=stratum )
	if (details) result <- c(result, list(A=Amats, V=V, n=n, m=m))
	class(result) <- "CMHtest"
	result
}

# do overall test, from a computed CMHtest list
CMHtest3 <- function(object,
	types=c("cor", "rmeans", "cmeans", "general")) 
{
	nstrat <- length(object)   # number of strata

	# extract components, each a list of nstrat terms
	n.list <- lapply(object, function(s) s$n)
	m.list <- lapply(object, function(s) s$m)
	V.list <- lapply(object, function(s) s$V)
	A.list <- lapply(object, function(s) s$A)
	nt <- sapply(lapply(object, function(s) s$n), sum)
	Df <- object[[1]]$table[,"Df"]

	if (length(types)==1 && types=="ALL") types <- c("general", "rmeans", "cmeans", "cor" )
	types <- match.arg(types, several.ok=TRUE)

	table <- list()
	for (type in types) {
		AVA <- 0
		Anm <- 0
		for (k in 1:nstrat) {
			A <- A.list[[k]][[type]]
			V <- V.list[[k]]
			n <- n.list[[k]]
			m <- m.list[[k]]
			AVA <- AVA + A %*% V %*% t(A)
			Anm <- Anm + A %*% (n-m)
			}
		Q <- t(Anm) %*% solve(AVA) %*% Anm
		df <- Df[type]
		pvalue <- pchisq(Q, df, lower.tail=FALSE)
		table <- rbind(table, c(Q, df, pvalue))
	}
	rownames(table) <- types
	colnames(table) <- c("Chisq", "Df", "Prob")
	xnames <- object[[1]]$names
	result=list(table=table, names=xnames, stratum="ALL")
	class(result) <- "CMHtest"
	result
}



# basic CMH calculation
cmh <- function(n, m,A, V, df) {
	AVA <- A %*% V %*% t(A)
	Q <- t(n-m) %*% t(A) %*% solve(AVA) %*% A %*% (n-m)
	pvalue <- pchisq(Q, df, lower.tail=FALSE)
	c(Q, df, pvalue)
}

# DONE: incorporate stratum name in the heading
# TODO: handle the printing of pvalues better

print.CMHtest <- function(x, digits = max(getOption("digits") - 2, 3), ...) {
	heading <- "Cochran-Mantel-Haenszel Statistics"
	if (!is.null(x$names)) heading <- paste(heading, "for", paste(x$names, collapse=" by "))
	if (!is.null(x$stratum)) heading <- paste(heading, 
				ifelse(x$stratum=="ALL", "\n\tOverall tests, controlling for all strata", paste("\n\tin stratum", x$stratum)))
	# TODO: determine score types (integer, midrank) for heading
	
	df <- x$table
	types <- rownames(df)
	labels <- list(cor="Nonzero correlation", rmeans="Row mean scores differ",
			cmeans="Col mean scores differ", general="General association")
	labels <- unlist(labels[types])  # select the labels for the types
	df <- data.frame("AltHypothesis"=as.character(labels), df, stringsAsFactors=FALSE)
	cat(heading,"\n\n")
	print(df, digits=digits, ...)
	cat("\n")
	
	invisible(x)
}

