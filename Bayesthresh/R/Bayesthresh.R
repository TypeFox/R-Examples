##### Arquivo modificado em 20 de março de 2013.
#####

### Funções obtidas da library lme4 (Douglas Bates)
### para a inserção da fórmula

### Utilities for parsing the mixed model formula

findbars <- function(term)
	### Return the pairs of expressions that separated by vertical bars
{
	if (is.name(term) || !is.language(term)) return(NULL)
	if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
	if (!is.call(term)) stop("term must be of class call")
	if (term[[1]] == as.name('|')) return(term)
	if (length(term) == 2) return(findbars(term[[2]]))
	c(findbars(term[[2]]), findbars(term[[3]]))
}

nobars <- function(term)
	### Return the formula omitting the pairs of expressions that are
	### separated by vertical bars
{
	if (!('|' %in% all.names(term))) return(term)
	if (is.call(term) && term[[1]] == as.name('|')) return(NULL)
	if (length(term) == 2) {
		nb <- nobars(term[[2]])
		if (is.null(nb)) return(NULL)
		term[[2]] <- nb
		return(term)
	}
	nb2 <- nobars(term[[2]])
	nb3 <- nobars(term[[3]])
	if (is.null(nb2)) return(nb3)
	if (is.null(nb3)) return(nb2)
	term[[2]] <- nb2
	term[[3]] <- nb3
	term
}

subbars <- function(term)
	### Substitute the '+' function for the '|' function
{
	if (is.name(term) || !is.language(term)) return(term)
	if (length(term) == 2) {
		term[[2]] <- subbars(term[[2]])
		return(term)
	}
	stopifnot(length(term) >= 3)
	if (is.call(term) && term[[1]] == as.name('|'))
		term[[1]] <- as.name('+')
	for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
	term
}

subnms <- function(term, nlist)
	### Substitute any names from nlist in term with 1
{
	if (!is.language(term)) return(term)
	if (is.name(term)) {
		if (any(unlist(lapply(nlist, get("=="), term)))) return(1)
		return(term)
	}
	stopifnot(length(term) >= 2)
	for (j in 2:length(term)) term[[j]] <- subnms(term[[j]], nlist)
	term
}

slashTerms <- function(x)
	### Return the list of '/'-separated terms in an expression that
	### contains slashes
{
	if (!("/" %in% all.names(x))) return(x)
	if (x[[1]] != as.name("/"))
		stop("unparseable formula for grouping factor")
	list(slashTerms(x[[2]]), slashTerms(x[[3]]))
}

makeInteraction <- function(x)
	### from a list of length 2 return recursive interaction terms
{
	if (length(x) < 2) return(x)
	trm1 <- makeInteraction(x[[1]])
	trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
	list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
}

expandSlash <- function(bb)
	### expand any slashes in the grouping factors returned by findbars
{
	if (!is.list(bb)) return(expandSlash(list(bb)))
	## I really do mean lapply(unlist(... - unlist returns a
	## flattened list in this case
	unlist(lapply(bb, function(x) {
				  if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
					  return(lapply(unlist(makeInteraction(trms)),
									function(trm) substitute(foo|bar,
															 list(foo = x[[2]],
																  bar = trm))))
				  x
}))
}

### Utilities used in lmer, glmer and nlmer

createCm <- function(A, s)
	### Create the nonzero pattern for the sparse matrix Cm from A.
	### ncol(A) is s * ncol(Cm).  The s groups of ncol(Cm) consecutive
	### columns in A are overlaid to produce Cm.
{
	stopifnot(is(A, "dgCMatrix"))
	s <- as.integer(s)[1]
	if (s == 1L) return(A)
	if ((nc <- ncol(A)) %% s)
		stop(gettextf("ncol(A) = %d is not a multiple of s = %d",
					  nc, s))
	ncC <- as.integer(nc / s)
	TA <- as(A, "TsparseMatrix")
	as(new("dgTMatrix", Dim = c(nrow(A), ncC),
		   i = TA@i, j = as.integer(TA@j %% ncC), x = TA@x),
	   "CsparseMatrix")
}

### FIXME: somehow the environment of the mf formula does not have
### .globalEnv in its parent list.  example(Mmmec, package = "mlmRev")
### used to have a formula of ~ offset(log(expected)) + ... and the
### offset function was not found in eval(mf, parent.frame(2))
lmerFrames <- function(mc, formula, contrasts, vnms = character(0))
	### Create the model frame, X, Y, wts, offset and terms

	### mc - matched call of calling function
	### formula - two-sided formula
	### contrasts - contrasts argument
	### vnms - names of variables to be included in the model frame
{
	mf <- mc
	m <- match(c("data", "subset", "weights", "na.action", "offset"),
			   names(mf), 0)
	mf <- mf[c(1, m)]

	## The model formula for evaluation of the model frame.  It looks
	## like a linear model formula but includes any random effects
	## terms and any names of parameters used in a nonlinear mixed model.
	frame.form <- subbars(formula)      # substitute `+' for `|'
	if (length(vnms) > 0)               # add the variables names for nlmer
		frame.form[[3]] <-
			substitute(foo + bar,
					   list(foo = parse(text = paste(vnms, collapse = ' + '))[[1]],
							bar = frame.form[[3]]))

	## The model formula for the fixed-effects terms only.
	fixed.form <- nobars(formula)       # remove any terms with `|'
	if (inherits(fixed.form, "name"))   # RHS is empty - use `y ~ 1'
		fixed.form <- substitute(foo ~ 1, list(foo = fixed.form))

	## attach the correct environment
	environment(fixed.form) <- environment(frame.form) <- environment(formula)

	## evaluate a model frame
	mf$formula <- frame.form
	mf$drop.unused.levels <- TRUE
	mf[[1]] <- as.name("model.frame")
	fe <- mf                            # save a copy of the call
	mf <- eval(mf, parent.frame(2))

	## evaluate the terms for the fixed-effects only (used in anova)
	fe$formula <- fixed.form
	fe <- eval(fe, parent.frame(2)) # allow model.frame to update them

	## response vector
	Y <- model.response(mf, "any")
	## avoid problems with 1D arrays, but keep names
	if(length(dim(Y)) == 1) {
		nm <- rownames(Y)
		dim(Y) <- NULL
		if(!is.null(nm)) names(Y) <- nm
	}
	mt <- attr(fe, "terms")

	## Extract X checking for a null model. This check shouldn't be
	## needed because an empty formula is changed to ~ 1 but it can't hurt.
	X <- if (!is.empty.model(mt))
		model.matrix(mt, mf, contrasts) else matrix(,NROW(Y),0)
	storage.mode(X) <- "double"      # when ncol(X) == 0, X is logical
	fixef <- numeric(ncol(X))
	names(fixef) <- colnames(X)
	dimnames(X) <- NULL

	## Extract the weights and offset.  For S4 classes we want the
	## `not used' condition to be numeric(0) instead of NULL
	wts <- model.weights(mf); if (is.null(wts)) wts <- numeric(0)
	off <- model.offset(mf); if (is.null(off)) off <- numeric(0)

	## check weights and offset
	if (any(wts <= 0))
		stop(gettextf("negative weights or weights of zero are not allowed"))
	if(length(off) && length(off) != NROW(Y))
		stop(gettextf("number of offsets is %d should equal %d (number of observations)",
					  length(off), NROW(Y)))

	## remove the terms attribute from mf
	attr(mf, "terms") <- mt
	list(Y = Y, X = X, wts = as.double(wts), off = as.double(off), mf = mf, fixef = fixef)
}

##' Is f1 nested within f2?
##'
##' Does every level of f1 occur in conjunction with exactly one level
##' of f2? The function is based on converting a triplet sparse matrix
##' to a compressed column-oriented form in which the nesting can be
##' quickly evaluated.
##'
##' @param f1 factor 1
##' @param f2 factor 2

##' @return TRUE if factor 1 is nested within factor 2

isNested <- function(f1, f2)
{
	f1 <- as.factor(f1)
	f2 <- as.factor(f2)
	stopifnot(length(f1) == length(f2))
	sm <- as(new("ngTMatrix",
				 i = as.integer(f2) - 1L,
				 j = as.integer(f1) - 1L,
				 Dim = c(length(levels(f2)),
						 length(levels(f1)))),
			 "CsparseMatrix")
	all(diff(sm@p) < 2)
}

##' dimsNames and devNames are in the package's namespace rather than
##' in the function lmerFactorList because the function sparseRasch
##' needs to access them.

dimsNames <- c("nt", "n", "p", "q", "s", "np", "LMM", "REML",
			   "fTyp", "lTyp", "vTyp", "nest", "useSc", "nAGQ",
			   "verb", "mxit", "mxfn", "cvg")
dimsDefault <- list(s = 1L)             # identity mechanistic model
#                    mxit= 300L,         # maximum number of iterations
#                    mxfn= 900L, # maximum number of function evaluations
#                    verb= 0L,           # no verbose output
#                    np= 0L,             # number of parameters in ST
#                    LMM= 0L,            # not a linear mixed model
#                    REML= 0L,         # glmer and nlmer don't use REML
#                    fTyp= 2L,           # default family is "gaussian"
#                    lTyp= 5L,           # default link is "identity"
#                    vTyp= 1L, # default variance function is "constant"
#                    useSc= 1L, # default is to use the scale parameter
#                    nAGQ= 1L,                  # default is Laplace
#                    cvg = 0L)                  # no optimization yet attempted

devNames <- c("ML", "REML", "ldL2", "ldRX2", "sigmaML",
			  "sigmaREML", "pwrss", "disc", "usqr", "wrss",
			  "dev", "llik", "NULLdev")


##' Create model matrices from r.e. terms.
##'
##' Create the list of model matrices from the random-effects terms in
##' the formula and the model frame.
##'
##' @param formula model formula
##' @param fr: list with '$mf': model frame; '$X': .. matrix
##' @param rmInt logical scalar - should the `(Intercept)` column
##'        be removed before creating Zt
##' @param drop logical scalar indicating if elements with numeric
##'        value 0 should be dropped from the sparse model matrices
##'
##' @return a list with components named \code{"trms"}, \code{"fl"}
##'        and \code{"dims"}
lmerFactorList <- function(formula, fr, rmInt, drop)
{
	mf <- fr$mf
	## record dimensions and algorithm settings

	## create factor list for the random effects
	bars <- expandSlash(findbars(formula[[3]]))
	if (!length(bars)) stop("No random effects terms specified in formula")
	names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
	fl <- lapply(bars,
				 function(x)
				 {
					 ff <- eval(substitute(as.factor(fac)[,drop = TRUE],
										   list(fac = x[[3]])), mf)
					 im <- as(ff, "sparseMatrix") # transpose of indicators
					 ## Could well be that we should rather check earlier .. :
					 if(!isTRUE(validObject(im, test=TRUE)))
						 stop("invalid conditioning factor in random effect: ", format(x[[3]]))

					 mm <- model.matrix(eval(substitute(~ expr, # model matrix
														list(expr = x[[2]]))),
										mf)
					 if (rmInt) {
						 if (is.na(icol <- match("(Intercept)", colnames(mm)))) break
						 if (ncol(mm) < 2)
							 stop("lhs of a random-effects term cannot be an intercept only")
						 mm <- mm[ , -icol , drop = FALSE]
					 }
					 ans <- list(f = ff,
								 A = do.call(rBind,
											 lapply(seq_len(ncol(mm)), function(j) im)),
								 Zt = do.call(rBind,
											  lapply(seq_len(ncol(mm)),
													 function(j) {im@x <- mm[,j]; im})),
								 ST = matrix(0, ncol(mm), ncol(mm),
											 dimnames = list(colnames(mm), colnames(mm))))
					 if (drop) {
						 ## This is only used for nlmer models.
						 ## Need to do something more complicated for A
						 ## here.  Essentially you need to create a copy
						 ## of im for each column of mm, im@x <- mm[,j],
						 ## create the appropriate number of copies,
						 ## prepend matrices of zeros, then rBind and drop0.
						 ans$A@x <- rep(0, length(ans$A@x))
						 ans$Zt <- drop0(ans$Zt)
					 }
					 ans
				 })
	dd <-
		VecFromNames(dimsNames, "integer",
					 c(list(n = nrow(mf), p = ncol(fr$X), nt = length(fl),
							q = sum(sapply(fl, function(el) nrow(el$Zt)))),
					   dimsDefault))
	## order terms by decreasing number of levels in the factor but don't
	## change the order if this is already true
	nlev <- sapply(fl, function(el) length(levels(el$f)))
	## determine the number of random effects at this point
	if (any(diff(nlev)) > 0) fl <- fl[rev(order(nlev))]
	## separate the terms from the factor list
	trms <- lapply(fl, "[", -1)
	names(trms) <- NULL
	fl <- lapply(fl, "[[", "f")
	attr(fl, "assign") <- seq_along(fl)
	## check for repeated factors
	fnms <- names(fl)
	if (length(fnms) > length(ufn <- unique(fnms))) {
		## check that the lengths of the number of levels coincide
		fl <- fl[match(ufn, fnms)]
		attr(fl, "assign") <- match(fnms, ufn)
	}
	names(fl) <- ufn
	## check for nesting of factors
	dd["nest"] <- all(sapply(seq_along(fl)[-1],
							 function(i) isNested(fl[[i-1]], fl[[i]])))

	list(trms = trms, fl = fl, dims = dd)
}

checkSTform <- function(ST, STnew)
	### Check that the 'STnew' argument matches the form of ST.
{
	stopifnot(is.list(STnew), length(STnew) == length(ST),
			  all.equal(names(ST), names(STnew)))
	lapply(seq_along(STnew), function (i)
		   stopifnot(class(STnew[[i]]) == class(ST[[i]]),
					 all.equal(dim(STnew[[i]]), dim(ST[[i]]))))
	all(unlist(lapply(STnew, function(m) all(diag(m) > 0))))
}

#######################################################################
########################################################################
#######################################################################
#######################################################################
#######################################################################

VecFromNames <- function(nms, mode = "numeric", defaults = list())
	### Generate a named vector of the given mode
{
	ans <- vector(mode = mode, length = length(nms))
	names(ans) <- nms
	ans[] <- NA
	if ((nd <- length(defaults <- as.list(defaults))) > 0) {
		if (length(dnms <- names(defaults)) < nd)
			stop("defaults must be a named list")
		stopifnot(all(dnms %in% nms))
		ans[dnms] <- as(unlist(defaults), mode)
	}
	ans
}

mkZt <- function(FL, start, s = 1L)
	### Create the standard versions of flist, Zt, Gp, ST, A, Cm,
	### Cx, and L. Update dd.
{
	dd <- FL$dims
	fl <- FL$fl
	asgn <- attr(fl, "assign")
	trms <- FL$trms
	ST <- lapply(trms, `[[`, "ST")
	Ztl <- lapply(trms, `[[`, "Zt")
	Zt <- do.call(rBind, Ztl)
	Zt@Dimnames <- vector("list", 2)
	#    Gp <- unname(c(0L, cumsum(sapply(Ztl, nrow))))
	#    .Call(mer_ST_initialize, ST, Gp, Zt)
	#    A <- do.call(rBind, lapply(trms, `[[`, "A"))
	#    rm(Ztl, FL)                         # because they could be large
	nc <- sapply(ST, ncol)         # of columns in els of ST
	#    Cm <- createCm(A, s)
	#    L <- .Call(mer_create_L, Cm)
	#    if (s < 2) Cm <- new("dgCMatrix")
	#    if (!is.null(start) && checkSTform(ST, start)) ST <- start

	nvc <- sapply(nc, function (qi) (qi * (qi + 1))/2) # no. of var. comp.
	### FIXME: Check number of variance components versus number of
	### levels in the factor for each term. Warn or stop as appropriate

	dd["np"] <- as.integer(sum(nvc))    # number of parameters in optimization
	dev <- VecFromNames(devNames, "numeric")
	fl <- do.call(data.frame, c(fl, check.names = FALSE))
	attr(fl, "assign") <- asgn

	list(ST = ST, Zt = Zt, dd = dd, dev = dev, flist = fl)
}


# For one variance component (vu*A) 
#
#  A = Variance and covariance matrix
#
#  vu = Prior of the variance component
#

Avu1 <- function(Zz, Dd, A, Vu, FL)
{
	V <- rbind(Zz, cbind(Dd, Vu*A))
	return(V)
}



# For more two variance compoents (vu*A)
Avuall <- function(Zz, Dd, A, Vu, FL)
{
	m <- dim(A)[1]
	n <- dim(A)[2]
	tc <- length(Vu)
	ifc <- unlist(lapply(FL$fl, function(x) length(levels(x))))
	il <-  1
	ic <- ifc[1]
	for(i in 2:tc){
		ic[i] <- ifc[i]+ic[i-1]
		il[i] <- ic[i]-ifc[i]+1
	}
	storage.mode(A) <- "double"
	Aux <- .Fortran("vua", as.double(Vu), A=A, as.integer(ic), as.integer(il), as.integer(m), 
					as.integer(n), as.integer(tc), PACKAGE="Bayesthresh")$A
	V <- rbind(Zz, cbind(Dd,Aux))
	return(V)
}


# Obtained an estimate for conditional variance 
#


UniComp <- function(naleat, ru, u, su, A = NULL, Vu = NULL, FL = NULL)
{
	c1 <- (naleat+2*ru)/2
	c2 <- (crossprod(u)+2*su)/2
	Vu <- rgamma(1,c1,c2)
	Vu <- as.double(1/Vu)
	return(Vu)
}

VarComp <- function(naleat, ru, u, su, A, Vu, FL)
{
	m <- dim(A)[1]
	tc <- length(Vu)
	ifc <- unlist(lapply(FL$fl, function(x) length(levels(x))))
	il <-  1
	ic <- ifc[1]
	for(i in 2:tc){
		ic[i] <- ifc[i]+ic[i-1]
		il[i] <- ic[i]-ifc[i]+1
	}
	storage.mode(A) <- "double"
	Vu <- .Fortran("part", A, as.integer(m), as.double(u), Vu=as.double(Vu),as.double(ru), 
				   as.double(su), as.integer(tc), as.integer(il),as.integer(ic), PACKAGE="Bayesthresh")$Vu
	return(Vu)
}



# WriT <- function(saida = NULL){}

# wriT <- function(theta,tau,vu,verossim,ef.iter)
# {
# Lout <- length(saida)
# write(saida, file="output.txt", ncolumns=Lout, append=TRUE)
#     op.mcmc <- NULL
#     op.mcmc$theta <- matrix(0,nrow=ef.iter, ncol=length(theta))
#     op.mcmc$tau <- matrix(0,nrow=ef.iter, ncol=length(tau))
#     op.mcmc$vu <- matrix(0,nrow=ef.iter, ncol=length(vu))
#     return(op.mcmc)
# }
# 
priors.control <- function(ru = 3, su = 5, dre = 20, dse = 5, method = c("AC", "MC", "NC")) 
{
	method <- match.arg(method)
	if(method == "AC" | method == "MC")
	{
		if(ru <= 0 || su <= 0) stop("Priors more zero")
		else
			prior <- list(ru=ru, su=su)
	}
	#    if(method == "MC")
	#      {
	#        if(se <= 0 || ru <= 0 || su <= 0) stop("Priores devem ser maior que zero")
	#         prior <- list(se=se, ru=ru, su=su)
	#      }
	else
	{
		if(ru <= 0 || su <= 0){
			if(dre <= 0 || dse <= 0){
				stop("Priors more zero")
			}
		}
		prior <-list(ru=ru, su=su, dre = dre, dse = dse)
	}
	return(prior)
}



algorithm <- function(method = c("AC", "MC", "NC") , link = c("Gaussian","t"))
{
	method <- match.arg(method)
	link <- match.arg(link)
	if(method == c("AC")){
		if(link == c("Gaussian"))list("ACGaussian", "AC")
		else
			list("ACt", "AC")
	}
	else if(method == c("MC")){
		if(link == c("Gaussian"))list("MCGaussian","MC")
		else
			list("MCt", "MC")
	}
	else{
		if(link == c("Gaussian"))list("NCGaussian", "NC")
		else
			list("NCt", "NC")
	}
}

#############################################################################
###
### Albert and Chib algorithm with Gaussian distribution for latent variable
###
############################################################################

ACGaussian <- function(Y,X,Z2,W,A,escala, ru, su, ntrat, nfixo, naleat,
					   burn, jump, ef.inter, Write = FALSE, FL,nomes,nomes2)
{

	### Initial random value

	N <- length(Y)
	L <- c(1:N)
	L <- sapply(L,function(x)qnorm(pnorm(0)+pnorm(0)*runif(1)))
	theta <- ginv(crossprod(W))%*%t(W)%*%L
	ct <- length(theta)
	e <- L - W%*%theta
	ve <- 1
	vu <- unlist(lapply(FL$fl, function(x) length(levels(x))))
	taunovo <- c(0,sort(runif((length(escala))-2,0,3)))
	tau <- taunovo
	delta2  <- 1
	vero <- rep(1,N)
	verossim <- 0  
	Ww <- crossprod(W)
	tW <- t(W)
	iA <- solve(A)
	rowmcmc <- 1
	output <- rep(0,length(c(theta,tau,vu,verossim)))

	Lout  <- length(output)
	sumsquare <- output

	if(length(vu) == 1){
		aVu <- match.fun(Avu1)
		Wpart <- match.fun(UniComp)
	}
	else {
		aVu <- match.fun(Avuall)
		Wpart <- match.fun(VarComp)
	}
	outp.mcmc <- NULL
	if(Write == TRUE){
		#         WFun <- match.fun(WriT)
		#     }
		#     else {
		wriT <- function(theta,tau,vu,nrow=ef.inter)
		{
			#     Lout <- length(saida)
			#     write(saida, file="output.txt", ncolumns=Lout, append=TRUE)
			op.mcmc <- NULL
			op.mcmc$theta <- matrix(0,nrow=nrow, ncol=length(theta))
			op.mcmc$tau <- matrix(0,nrow=nrow, ncol=length(tau))
			op.mcmc$vu <- matrix(0,nrow=nrow, ncol=length(vu))
			return(op.mcmc)
		}
		outp.mcmc <- wriT(theta,tau,vu,nrow=ef.inter)
		#         write(nomes2, file="output.txt", ncolumns=length(nomes2))
	}


	### Auxiliar step for construction W matrix

	Zz  <- cbind(1000*diag(nfixo), matrix(0,nfixo,naleat))
	Dd  <- matrix(0,naleat,nfixo)

	verossim <- 0

	### Initialize sample posterior
	iter <- burn + jump*ef.inter

	for(cont in 1:iter){

		### Conditional for theta
		V <- aVu(Zz, Dd, A, vu, FL)
		S <- solve(V)
		Iw <- solve(Ww + S)
		theta <- Iw%*%tW%*%L
		var <- ve*Iw
		thetaest <- rmvnorm(1,theta,var, method="chol")
		theta <- thetaest
		beta <- theta[1:nfixo]
		u <- theta[(nfixo+1):ct]

		### Conditional for variance of random effects

		vu <- Wpart(naleat, ru, u, su, A, vu, FL)

		### Conditional for L

		m <- tcrossprod(W,theta)

		### Cummulative density for latent variable
		probacum  <- .C("acumN", as.integer(Y), L=as.double(L), as.integer(escala), 
						as.integer(max(escala)), as.integer(min(escala)), as.integer(length(escala)), 
						as.integer(N), as.double(vero), as.double(m), as.double(tau), as.double(ve), 
						verossim=as.double(verossim), PACKAGE="Bayesthresh")[c("L","verossim")]

		L    <- probacum$L
		verossim <- probacum$verossim
		#      if(verossim == 0) verossim <- 1e-320

		### Conditional for thresholds parameters
		min <- as.vector(tapply(L,Y,min))
		max <- as.vector(tapply(L,Y,max))


		for(j in 2:(length(escala)-1)) {
			tau[j] <- runif(1,max[j],min[j+1])
		}


		### Print output

		if(cont > burn && (cont-burn)%%jump == 0) {
			saida <- c(theta,tau,vu,verossim)
			output <- output + saida
			sumsquare <- sumsquare + (saida^2)
			if(Write==TRUE){
				outp.mcmc[[1]][rowmcmc,] <- theta
				outp.mcmc[[2]][rowmcmc,] <- tau
				outp.mcmc[[3]][rowmcmc,] <- vu
				rowmcmc <- rowmcmc + 1
			}
		}
	}
	list(Sum=output, SumSquare = sumsquare,outp.mcmc=outp.mcmc)
}


###############################################################################
###
###
### Albert and Chib algorithm with t-Student distribution for latent variable
###
###
##############################################################################


ACt <- function(Y,X,Z2,W,A,escala, ru, su, ntrat, nfixo, naleat, burn, jump, 
				ef.inter, Write = FALSE, FL,nomes,nomes2)
{

	### Arbitrary initial values
	N <- length(Y)
	L <- qnorm(pnorm(0)+pnorm(0)*runif(N))
	theta <- ginv(crossprod(W))%*%t(W)%*%L
	ct <- length(theta)
	e <- L - W%*%theta
	ve <- 1
	vu <- unlist(lapply(FL$fl, function(x) length(levels(x))))
	taunovo <- c(0,sort(runif((length(escala))-2,0,3)))
	tau <- taunovo
	pnuc <- 0.5
	pnu <- 0.5
	nu <- 10
	lambda <- rep(1,N)
	R <- diag(N)
	dw1 <- dim(W)[1]
	dw2 <- dim(W)[2]

	### Object necessary for estimate of likelihood
	vero <- rep(1,N)
	verossim <- 0  

	output <- rep(0,length(c(theta,tau,vu,verossim)))
	Lout  <- length(output)
	sumsquare <- output

	if(length(vu) == 1){
		aVu <- match.fun(Avu1)
		Wpart <- match.fun(UniComp)
	}
	else {
		aVu <- match.fun(Avuall)
		Wpart <- match.fun(VarComp)
	}
	outp.mcmc <- NULL
	if(Write == TRUE){
		#         WFun <- match.fun(WriT)
		#     }
		#     else {
		wriT <- function(theta,tau,vu,nrow=ef.inter)
		{
			#     Lout <- length(saida)
			#     write(saida, file="output.txt", ncolumns=Lout, append=TRUE)
			op.mcmc <- NULL
			op.mcmc$theta <- matrix(0,nrow=nrow, ncol=length(theta))
			op.mcmc$tau <- matrix(0,nrow=nrow, ncol=length(tau))
			op.mcmc$vu <- matrix(0,nrow=nrow, ncol=length(vu))
			return(op.mcmc)
		}
		outp.mcmc <- wriT(theta,tau,vu,nrow=ef.inter)
		#         write(nomes2, file="output.txt", ncolumns=length(nomes2))
	}
	### Auxiliary step of W matrix
	Zz  <- cbind(1000*diag(nfixo), matrix(0,nfixo,naleat))
	Dd  <- matrix(0,naleat,nfixo)
	verossim <- 0
	tW <- t(W)
	iA <- solve(A)

	### Sampling of posteriori
	iter <- burn + jump*ef.inter
	rowmcmc <- 1

	for(cont in 1:iter){

		### Conditional for theta
		V <- aVu(Zz, Dd, A, vu, FL)
		S <- ve*(solve(V))
		Iw <- solve(crossprod(W,R)%*%W + S)
		theta <- Iw%*%tW%*%R%*%L
		var <- ve*Iw
		thetaest <- rmvnorm(1,theta,var, method="chol")
		theta <- thetaest
		beta <- theta[1:nfixo]
		u <- theta[(nfixo+1):ct]

		### Conditional for variance of random effects

		vu <- Wpart(naleat, ru, u, su, A, vu, FL)

		### Sampling lambda (elementwise)
		storage.mode(W) <- "double"
		lambda <- .Fortran("gslambda", lambda=as.double(lambda), as.double(nu), as.double(L), 
						   W, as.double(theta), as.integer(dw1), as.integer(dw2), PACKAGE="Bayesthresh")$lambda
		R <- diag(lambda)

		### Passo Metropolis-Hastings para amostrar nu
		nuc <- 2
		while(nuc < 3) {
			nuc <- rpois(1,nu)
		}

		#      pnu  <- N*((nu/2)*log(nu/2)-lgamma(nu/2))+sum((nu/2-1)*log(lambda)+(-(nu/2)*lambda))-2*log(1+nu)
		#      pnuc <- N*((nuc/2)*log(nuc/2)-lgamma(nuc/2))+sum((nuc/2-1)*log(lambda)+(-(nuc/2)*lambda))-2*log(1+nuc)
		pnu  <- N*((nu/2)*(nu/2)-lgamma(nu/2))+sum((nu/2-1)*(lambda)+(-(nu/2)*lambda))-2*(1+nu)
		pnuc <- N*((nuc/2)*(nuc/2)-lgamma(nuc/2))+sum((nuc/2-1)*(lambda)+(-(nuc/2)*lambda))-2*(1+nuc)


		#      pnu  <- exp(pnu)
		#      pnuc <- exp(pnuc)
		dnu  <- dpois(nu,nuc)
		dnuc <- dpois(nuc,nu)

		num   <- pnuc*dnu
		denom <- pnu*dnuc


		ifelse(num > denom, alfa <- 1, alfa <- num/denom)
		if(runif(1) < alfa) nu <- nuc

		### Condicional para L
		mw <- tcrossprod(W,theta)
		dp <- sqrt(ve/lambda)

		probacum  <- .C("acumt", as.integer(Y), L=as.double(L), as.integer(escala), as.integer(max(escala)),
						as.integer(min(escala)), as.integer(length(escala)), as.integer(N), as.double(vero),
						as.double(mw), as.double(tau), as.double(dp), verossim=as.double(verossim),
						PACKAGE="Bayesthresh")[c("L","verossim")]

		L <- probacum$L
		verossim <- probacum$verossim
		#      if(verossim == 0) verossim <- 1e-320

		### Amostragem dos thresholds
		mini <-  as.vector(tapply(L,Y,min))
		maxi <-  as.vector(tapply(L,Y,max))
		for(j in 2:(length(escala)-1)) {
			tau[j] <- runif(1, maxi[j], mini[j+1])
		}


		### Gera a saida no arquivo cadeia.txt
		if(cont > burn && (cont-burn)%%jump == 0) {
			saida <- c(theta,tau,vu,verossim)
			output <- output + saida
			sumsquare <- sumsquare + (saida^2)
			if(Write==TRUE){
				outp.mcmc[[1]][rowmcmc,] <- theta
				outp.mcmc[[2]][rowmcmc,] <- tau
				outp.mcmc[[3]][rowmcmc,] <- vu
				rowmcmc <- rowmcmc + 1
			}
		}
	}
	list(Sum=output, SumSquare = sumsquare,outp.mcmc=outp.mcmc)
}


#############################################################
###
###
### Algoritimo MC e variavel latente com distribuicao normal
###
###
############################################################ 


MCGaussian <- function(Y,X,Z2,W,A,escala, ru, su, ntrat, nfixo, naleat, burn, jump, 
					   ef.inter, Write = FALSE, FL,nomes,nomes2)
{

	### Valores arbitrarios iniciais
	N <- length(Y)
	L <- qnorm(pnorm(0)+pnorm(0)*runif(N))
	theta <- ginv(crossprod(W))%*%t(W)%*%L
	ct <- length(theta)
	e <- L - W%*%theta
	ve <- 1
	vu <- unlist(lapply(FL$fl, function(x) length(levels(x))))
	taunovo <- c(0,sort(runif((length(escala))-2,0.5,3.5)),1000)
	tau <- taunovo
	ftau <- runif(N,0.5,1.5)
	ftaunovo <- ftau
	dw1 <- dim(W)[1]
	dw2 <- dim(W)[2]

	### Iteracoes totais

	iter <- burn + jump*ef.inter
	rowmcmc <- 1

	### Objetos necessarios para o calculo da verossimilhanca
	vero <- rep(1,N)
	lvero <- log(vero)
	verossim <- 1


	### Matrix com os objetos p1, p2, ..., pn.
	ce <- length(escala)
	mp  <- matrix(0,ce-2,ce-1)
	nl  <- nrow(mp)
	nc  <- ncol(mp)

	### Desvio-padrao para a candidata
	dp <- 0.15
	d  <- rep(1,iter) # Vetor para atualizacao do dp

	output <- rep(0,length(c(theta,tau[-(length(tau))],vu,verossim)))
	Lout  <- length(output)
	sumsquare <- output

	if(length(vu) == 1){
		aVu <- match.fun(Avu1)
		Wpart <- match.fun(UniComp)
	}
	else {
		aVu <- match.fun(Avuall)
		Wpart <- match.fun(VarComp)
	}
	outp.mcmc <- NULL
	if(Write == TRUE){
		#         WFun <- match.fun(WriT)
		#     }
		#     else {
		wriT <- function(theta,tau,vu,nrow=ef.inter)
		{
			#     Lout <- length(saida)
			#     write(saida, file="output.txt", ncolumns=Lout, append=TRUE)
			op.mcmc <- NULL
			op.mcmc$theta <- matrix(0,nrow=nrow, ncol=length(theta))
			op.mcmc$tau <- matrix(0,nrow=nrow, ncol=(length(tau)-1))
			op.mcmc$vu <- matrix(0,nrow=nrow, ncol=length(vu))
			return(op.mcmc)
		}
		outp.mcmc <- wriT(theta,tau,vu,nrow=ef.inter)
		#         write(nomes2, file="output.txt", ncolumns=length(nomes2))
	}

	### Passo auxiliar na construcao da matriz W
	Zz  <- cbind(1000*diag(nfixo), matrix(0,nfixo,naleat))
	Dd  <- matrix(0,naleat,nfixo)
	verossim <- 0
	tW <- t(W)
	Ww <- crossprod(W)
	iA <- solve(A)


	### Laco de Amostragem da Posteriori


	for(cont in 1:iter){

		### Condicional para theta
		V <- aVu(Zz, Dd, A, vu, FL)
		S <- ve*(solve(V))
		Iw <- solve(Ww + S)
		theta <- Iw%*%tW%*%L
		var <- ve*Iw
		thetaest <- rmvnorm(1,theta,var, method="chol")
		theta <- thetaest
		beta <- theta[1:nfixo]
		u <- theta[(nfixo+1):ct]
		m <- tcrossprod(W,theta)

		### Condicional para variancia dos efeitos aleatorios

		vu <- Wpart(naleat, ru, u, su, A, vu, FL)

		############# AMOSTRAGEM DOS THRESHOLDS ###########


		taunovo <- .C("taunN", taunovo=as.double(taunovo), as.double(tau), as.double(dp), as.integer(ce),
					  PACKAGE="Bayesthresh")$taunovo


		### Matrix com as probabilidades de aceitacao do tau/taunovo 

		mp <- matrix(.C("calcpN", as.double(tau), as.double(taunovo), as.double(dp),
						p=as.double(vec(mp)), as.integer(nl), PACKAGE="Bayesthresh")$p, nrow = ce-2)


		### Funcao para estimar a probabilidade de tau/taunovo

		probtau <- .C("probtauN", as.integer(Y), ftau=as.double(ftau), ftaunovo=as.double(ftaunovo),
					  as.integer(escala), as.integer(max(escala)), as.integer(min(escala)), as.integer(ce),
					  as.integer(N), as.double(m), as.double(tau), as.double(taunovo),
					  PACKAGE="Bayesthresh")[c("ftau","ftaunovo")]

		ftau     <- probtau$ftau
		ftaunovo <- probtau$ftaunovo

		mpvet <- mp[,3]-mp[,4]
		mpvet[mpvet <= 0.00] <- 1E-45

		alfa  <- exp(sum(log(mp[,1]-mp[,2]) - log(mpvet)) + sum(log(ftaunovo) - log(ftau)))
		R     <- min(1,alfa)

		if(runif(1) < R) {
			tau <- taunovo

			### Funcao para a densidade acumulada da variavel Latente

			probacum  <- .C("acumN", as.integer(Y), L=as.double(L), as.integer(escala), as.integer(max(escala)),
							as.integer(min(escala)), as.integer(ce), as.integer(N), vero=as.double(vero),
							as.double(m), as.double(tau), as.double(ve), verossim=as.double(verossim),
							PACKAGE="Bayesthresh")[c("L","verossim")]

			L    <- probacum$L
			verossim <- probacum$verossim
			#   if(verossim == 0) verossim <- 1e-320
		}


		### Atualizacao do desvio-padrao da candidata

		d[cont] <- c(tau[2])
		ifelse(cont > 20, dp <- sd(d[(cont-20):cont]), dp <- 0.25)


		### Gera a saida no arquivo cadeia.txt
		if(cont > burn && (cont-burn)%%jump == 0) {
			saida <- c(theta,tau[-(length(tau))],vu,verossim)
			output <- output + saida
			sumsquare <- sumsquare + (saida^2)
			if(Write==TRUE){
				outp.mcmc[[1]][rowmcmc,] <- theta
				outp.mcmc[[2]][rowmcmc,] <- tau[-(length(tau))]
				outp.mcmc[[3]][rowmcmc,] <- vu
				rowmcmc <- rowmcmc + 1
			}
		}
	}
	list(Sum=output, SumSquare = sumsquare,outp.mcmc=outp.mcmc)
}

########################################################
###
###
### Algoritimo MC e variavel latente com distribuicao t
###
###
########################################################


MCt <- function(Y,X,Z2,W,A,escala, ru, su, ntrat, nfixo, naleat, burn, jump,
				ef.inter, Write = FALSE, FL, nomes, nomes2)
{

	### Valores arbitrarios iniciais
	ce <- length(escala)
	N <- length(Y)
	L <- qnorm(pnorm(0)+pnorm(0)*runif(N))
	theta <- ginv(crossprod(W))%*%t(W)%*%L
	ct <- length(theta)
	e <- L - W%*%theta
	ve <- 1
	vu <- unlist(lapply(FL$fl, function(x) length(levels(x))))
	taunovo <- c(0,sort(runif((length(escala))-2,0.5,3.5)),1000)
	tau <- taunovo
	ftau <- rep(1,N)
	ftaunovo <- ftau
	lambda <- ftau
	nu <- 10
	nuc <- ce
	pnu <- 0.5
	pnuc <- pnu
	R <- diag(N)
	d <- NULL
	dw1 <- dim(W)[1]
	dw2 <- dim(W)[2]

	### Iteracoes totais

	iter <- burn + jump*ef.inter
	rowmcmc <- 1

	### Objetos necessarios para o calculo da verossimilhanca
	vero <- rep(1,N)
	verossim <- 1

	mp  <- matrix(0,ce-2,ce-1)
	nl  <- nrow(mp)
	nc  <- ncol(mp)

	### Desvio-padrao para a candidata
	dpc <- 0.15
	dc  <- rep(1,iter) # Vetor para atualizacao do dp

	output <- rep(0,length(c(theta,tau[-(length(tau))],vu,verossim)))
	Lout  <- length(output)
	sumsquare <- output


	### Para vu*A
	if(length(vu) == 1){
		aVu <- match.fun(Avu1)
		Wpart <- match.fun(UniComp)
	}
	else {
		aVu <- match.fun(Avuall)
		Wpart <- match.fun(VarComp)
	}
	outp.mcmc <- NULL
	if(Write == TRUE){
		#         WFun <- match.fun(WriT)
		#     }
		#     else {
		wriT <- function(theta,tau,vu,nrow=ef.inter)
		{
			#     Lout <- length(saida)
			#     write(saida, file="output.txt", ncolumns=Lout, append=TRUE)
			op.mcmc <- NULL
			op.mcmc$theta <- matrix(0,nrow=nrow, ncol=length(theta))
			op.mcmc$tau <- matrix(0,nrow=nrow, ncol=(length(tau)-1))
			op.mcmc$vu <- matrix(0,nrow=nrow, ncol=length(vu))
			return(op.mcmc)
		}
		outp.mcmc <- wriT(theta,tau,vu,nrow=ef.inter)
		#         write(nomes2, file="output.txt", ncolumns=length(nomes2))
	}

	### Passo auxiliar na construcao da matriz W
	Zz  <- cbind(1000*diag(nfixo), matrix(0,nfixo,naleat))
	Dd  <- matrix(0,naleat,nfixo)
	verossim <- 0
	tW <- t(W)
	Ww <- crossprod(W)
	iA <- solve(A)


	### Laco de Amostragem da Posteriori


	for(cont in 1:iter){

		### Condicional para theta
		V <- aVu(Zz, Dd, A, vu, FL)
		S        <- ve*solve(V)
		Iw       <- solve(tW%*%R%*%W + S)
		theta    <- Iw%*%tW%*%R%*%L
		var      <- ve*Iw
		thetaest <- rmvnorm(1,theta, var, method="chol")
		theta    <- thetaest
		beta     <- theta[1:nfixo]
		u        <- theta[(nfixo+1):length(theta)]
		mw       <- tcrossprod(W,theta)

		### Condicional para variancia de blocos

		vu <- Wpart(naleat, ru, u, su, A, vu, FL)

		############## AMOSTRAGEM DOS THRESHOLDS ###############


		taunovo <- .C("taunt", taunovo=as.double(taunovo), as.double(tau), as.double(dpc),
					  as.integer(ce-1), PACKAGE="Bayesthresh")$taunovo

		### Matrix com as probabilidades de aceitacao do tau/taunovo 

		mp <- matrix(.C("calcpt", as.double(tau), as.double(taunovo), as.double(dpc),
						p=as.double(vec(mp)), as.integer(nl), PACKAGE="Bayesthresh")$p, nrow = ce-2)

		### Amostra lambda (GS elemento a elemento)

		storage.mode(W) <- "double"
		lambda <- .Fortran("gslambda", lambda=as.double(lambda), as.double(nu), as.double(L), 
						   W, as.double(theta), as.integer(dw1), as.integer(dw2), PACKAGE="Bayesthresh")$lambda
		R      <- diag(lambda)
		dp     <- sqrt(ve/lambda)

		### Funcao para estimar a probabilidade de tau/taunovo - TESTADA E APROVADA!!!

		probtau <- .C("probtaut", as.integer(Y), ftau=as.double(ftau), ftaunovo=as.double(ftaunovo),
					  as.integer(escala), as.integer(max(escala)), as.integer(min(escala)), as.integer(ce), 
					  as.integer(N), as.double(mw), as.double(tau), as.double(taunovo), 
					  as.double(dp), PACKAGE="Bayesthresh")[c("ftau","ftaunovo")]

		ftau     <- probtau$ftau
		ftaunovo <- probtau$ftaunovo

		mpvet <- mp[,3]-mp[,4]
		mpvet[mpvet <= 0.00] <- 1E-45

		alfa  <- exp(sum(log((mp[,1]-mp[,2])) - log(mpvet)) + sum(log(ftaunovo) - log(ftau)))

		D <- min(1,alfa)

		if(runif(1) < D) {
			tau <- taunovo

			### Amostra nu (Metropolis-Hastings)
			nuc <- 2
			while(nuc < 3){
				nuc <- rpois(1,nu)
			}
			# 
			# pnu   <- exp(N*((nu/2)*log(nu/2)-lgamma(nu/2))+sum(((nu/2)-1)*log(lambda)+(-(nu/2)*lambda)) - 2*log(1+nu))
			# pnuc  <- exp(N*((nuc/2)*log(nuc/2)-lgamma(nuc/2))+sum(((nuc/2)-1)*log(lambda)+(-(nuc/2)*lambda)) - 2*log(1+nuc))

			pnu   <- N*((nu/2)*(nu/2)-lgamma(nu/2))+sum(((nu/2)-1)*(lambda)+(-(nu/2)*lambda)) - 2*(1+nu)
			pnuc  <- N*((nuc/2)*(nuc/2)-lgamma(nuc/2))+sum(((nuc/2)-1)*(lambda)+(-(nuc/2)*lambda)) - 2*(1+nuc)


			dnu   <- dpois(nu,nu)
			dnuc  <- dpois(nuc,nuc)

			num   <- pnuc*dnu
			denom <- pnu*dnuc

			ifelse((num > denom), alfa <- 1, (alfa <- as.double(num/denom)))
			if(runif(1) < alfa) { nu <- nuc }


			### Funcao para a densidade acumulada da variavel Latente
			probacum  <- .C("acumt", as.integer(Y), L=as.double(L), as.integer(escala), as.integer(max(escala)), 
							as.integer(min(escala)), as.integer(length(escala)), as.integer(N), as.double(vero),
							as.double(mw), as.double(tau), as.double(dp), verossim=as.double(verossim),
							PACKAGE="Bayesthresh")[c("L","verossim")]

			L        <- probacum$L
			verossim <- probacum$verossim
			#    if(verossim == 0) verossim <- 1e-320
		}


		### Atualizacao do desvio-padrao da candidata

		d[cont] <- c(tau[3])
		ifelse(cont > 20, dpc <- sd(d[(cont-20):cont]), dpc <- 0.25)


		### Gera a saida no arquivo cadeia.txt
		if(cont > burn && (cont-burn)%%jump == 0) {
			saida <- c(theta,tau[-(length(tau))],vu,verossim)
			output <- output + saida
			sumsquare <- sumsquare + (saida^2)
			if(Write==TRUE){
				outp.mcmc[[1]][rowmcmc,] <- theta
				outp.mcmc[[2]][rowmcmc,] <- tau[-(length(tau))]
				outp.mcmc[[3]][rowmcmc,] <- vu
				rowmcmc <- rowmcmc + 1
			}
		}
	}
	list(Sum=output, SumSquare = sumsquare,outp.mcmc=outp.mcmc)
}

#############################################################
###
###
### Algoritimo NC e variavel latente com distribuicao normal
###
###
############################################################# 

NCGaussian <- function(Y,X,Z2,W,A,escala, ru, su, dre, dse, ntrat, nfixo, naleat, burn,
					   jump, ef.inter, Write = FALSE, FL,nomes,nomes2)
{
	### Valores arbitrarios iniciais
	ce <- length(escala)
	N <- length(Y)
	L <- qnorm(pnorm(0)+pnorm(0)*runif(N))
	theta <- ginv(crossprod(W))%*%t(W)%*%L
	ct <- length(theta)
	e <- L - W%*%theta
	ve <- 1
	vu <- unlist(lapply(FL$fl, function(x) length(levels(x))))
	taunovo <- c(0,sort(runif((length(escala))-2,0.5,3.5)),1000)
	nk <- tapply(L,Y,length)
	tau <- taunovo
	ftau <- rep(1,N)
	ftaunovo <- ftau
	delta2 <- 0.5
	delta <- sqrt(delta2)
	pnovo    <- sort(nk[c(-1,-2)]/10)
	p        <- pnovo
	nomes[length(nomes)+1] <- c("delta")

	### Iteracoes totais

	iter <- burn + jump*ef.inter
	rowmcmc <- 1


	### Objetos necessarios para o calculo da verossimilhanca
	vero <- rep(1,N)
	verossim <- 1

	output <- rep(0,length(c(theta,tau[-(length(tau))],vu,verossim,delta2)))
	Lout  <- length(output)
	sumsquare <- output

	if(length(vu) == 1){
		aVu <- match.fun(Avu1)
		Wpart <- match.fun(UniComp)
	}
	else {
		aVu <- match.fun(Avuall)
		Wpart <- match.fun(VarComp)
	}
	outp.mcmc <- NULL
	if(Write == TRUE){
		#         WFun <- match.fun(WriT)
		#     }
		#     else {
		wriT <- function(theta,tau,vu,nrow=ef.inter)
		{
			#     Lout <- length(saida)
			#     write(saida, file="output.txt", ncolumns=Lout, append=TRUE)
			op.mcmc <- NULL
			op.mcmc$theta <- matrix(0,nrow=nrow, ncol=length(theta))
			op.mcmc$tau <- matrix(0,nrow=nrow, ncol=(length(tau)-1))
			op.mcmc$vu <- matrix(0,nrow=nrow, ncol=length(vu))
			return(op.mcmc)
		}
		outp.mcmc <- wriT(theta,tau,vu,nrow=ef.inter)
		#         write(nomes2, file="output.txt", ncolumns=length(nomes2))
	}
	### Passo auxiliar na construcao da matriz W
	Zz  <- cbind(1000*diag(nfixo), matrix(0,nfixo,naleat))
	Dd  <- matrix(0,naleat,nfixo)
	tW <- t(W)
	Ww <- crossprod(W)
	iA <- solve(A)


	### Laco de Amostragem da Posteriori


	for(cont in 1:iter){

		### Condicional para theta
		V <- aVu(Zz, Dd, A, vu, FL)
		Iv       <- solve(V)
		S        <- delta2*Iv
		Iw       <- solve(Ww + S)
		theta    <- Iw%*%tW%*%L
		vari     <- delta2*Iw
		theta    <- rmvnorm(1,theta, vari)
		beta     <- theta[1:nfixo]
		u        <- theta[(nfixo+1):(length(theta))]

		### Condicional para variancia de blocos

		vu <- Wpart(naleat, ru, u, su, A, vu, FL)

		### Parametros para acelerar a convergencia (reparametrizacao)
		npWT <- L-tcrossprod(W,theta)
		np1 <- crossprod(npWT)
		Ivtheta <- theta%*%Iv
		np2 <- tcrossprod(Ivtheta,theta)

		delta2 <- 1/rgamma(1,(N+naleat+ce+2*dre)/2,(np1+2*dse+np2)/2)

		delta  <- sqrt(delta2)

		### Amostragem dos parametros de limiar

		nkp   <- 0.8*nk[2:(ce-1)]
		pnovo <- c(rdiric(1,shape=nkp))

		for(j in 2:(ce-1)) taunovo[j] <- sum(pnovo[1:(j-1)])


		lfp     <- nkp*log(p)
		lfpnovo <- nkp*log(pnovo)

		mw  <- tcrossprod(W,theta)


		probtau <- .C("probtauN", as.integer(Y), ftau=as.double(ftau), ftaunovo=as.double(ftaunovo), 
					  as.integer(escala), as.integer(max(escala)), as.integer(min(escala)), as.integer(ce),
					  as.integer(N), as.double(mw), as.double(tau), as.double(taunovo), as.double(delta), 
					  PACKAGE="Bayesthresh")[c("ftau","ftaunovo")]

		ftau     <- probtau$ftau
		ftaunovo <- probtau$ftaunovo

		lftau     <- log(ftau)
		lftaunovo <- log(ftaunovo)

		lftau[is.infinite(lftau)] <- log(1e-15)
		lftaunovo[is.infinite(lftaunovo)] <- log(1e-20)


		alfa      <- exp(sum(lftaunovo-lftau)-sum(lfpnovo-lfp))

		if(runif(1) < (min(1,alfa))){
			p   <- pnovo
			tau <- taunovo

			### Condicional para L

			### Funcao para a densidade acumulada da variavel Latente

			probacum  <- .C("acumN", as.integer(Y), L=as.double(L), as.integer(escala), as.integer(max(escala)),
							as.integer(min(escala)), as.integer(ce), as.integer(N), vero=as.double(vero),
							as.double(mw),	as.double(tau), as.double(delta), verossim=as.double(verossim),
							PACKAGE="Bayesthresh")[c("L","verossim")]


			L    <- probacum$L
			verossim <- probacum$verossim
			#     if(verossim == 0) verossim <- 1e-320
		}

		### Gera a saida no arquivo output.txt
		if(cont > burn && (cont-burn)%%jump == 0) {
			saida <- c(theta,tau[-(length(tau))],vu,verossim,delta2)
			output <- output + saida
			sumsquare <- sumsquare + (saida^2)
			if(Write==TRUE){
				outp.mcmc[[1]][rowmcmc,] <- theta
				outp.mcmc[[2]][rowmcmc,] <- tau[-(length(tau))]
				outp.mcmc[[3]][rowmcmc,] <- vu
				rowmcmc <- rowmcmc + 1
			}
		}
	}
	list(Sum=output, SumSquare = sumsquare,outp.mcmc=outp.mcmc)
}

#######################################################
###
###
### Algoritimo NC e variavel latente com distribuicao t
###
###
####################################################### 

NCt <- function(Y,X,Z2,W,A,escala, ru, su, dre, dse, ntrat, nfixo, naleat, burn, 
				jump, ef.inter, Write = FALSE,FL,nomes,nomes2)
{
	### Valores arbitrarios iniciais
	ce <- length(escala)
	N <- length(Y)
	L <- qnorm(pnorm(0)+pnorm(0)*runif(N))
	theta <- ginv(crossprod(W))%*%t(W)%*%L
	ct <- length(theta)
	e <- L - W%*%theta
	ve <- 1
	vu <- unlist(lapply(FL$fl, function(x) length(levels(x))))
	taunovo <- c(0,sort(runif((length(escala))-2,0.5,3.5)),1000)
	nk <- tapply(L,Y,length)
	tau <- taunovo
	ftau <- rep(1,N)
	ftaunovo <- ftau
	delta2 <- 0.5
	delta <- sqrt(delta2)
	pnovo    <- sort(nk[c(-1,-2)]/10)
	p        <- pnovo
	lambda <- ftau
	nu <- 10
	nuc <- ce
	pnu <- 0.5
	pnuc <- pnu
	R <- diag(lambda)
	nomes[length(nomes)+1] <- c("delta")

	### Iteracoes totais

	iter <- burn + jump*ef.inter
	rowmcmc <- 1


	### Objetos necessarios para o calculo da verossimilhanca
	vero <- rep(1,N)
	verossim <- 1

	output <- rep(0,length(c(theta,tau[-(length(tau))],vu,verossim,delta2)))
	Lout  <- length(output)
	sumsquare <- output

	if(length(vu) == 1){
		aVu <- match.fun(Avu1)
		Wpart <- match.fun(UniComp)
	}
	else {
		aVu <- match.fun(Avuall)
		Wpart <- match.fun(VarComp)
	}
	outp.mcmc <- NULL
	if(Write == TRUE){
		wriT <- function(theta,tau,vu,nrow=ef.inter)
		{
			op.mcmc <- NULL
			op.mcmc$theta <- matrix(0,nrow=nrow, ncol=length(theta))
			op.mcmc$tau <- matrix(0,nrow=nrow, ncol=(length(tau)-1))
			op.mcmc$vu <- matrix(0,nrow=nrow, ncol=length(vu))
			return(op.mcmc)
		}
		outp.mcmc <- wriT(theta,tau,vu,nrow=ef.inter)
	}

	### Passo auxiliar


	Wm  <- cbind(1000*diag(nfixo), matrix(0,nfixo,naleat))
	Dd  <- matrix(0,naleat,nfixo)
	Wt  <- t(W)
	Ag  <- solve(A)
	verossim <- 0

	### Laco de Amostragem da Posteriori (Gibbs Sampling)


	for(cont in 1:iter){

		### Condicional para theta
		V <- aVu(Wm, Dd, A, vu, FL)
		Iv       <- ginv(V)
		S        <- delta2*Iv
		Wrw      <- solve(crossprod(W,R)%*%W + S)
		theta    <- Wrw%*%Wt%*%R%*%L
		vari     <- delta2*Wrw
		theta    <- rmvnorm(1,theta, vari, method="chol")
		beta     <- theta[1:nfixo]
		u        <- theta[(nfixo+1):(length(theta))]

		### Condicional para variancia de blocos

		vu <- Wpart(naleat, ru, u, su, A, vu, FL)

		### Parametros para acelerar a convergencia (reparametrizacao)

		npWT <- L-tcrossprod(W,theta)
		np1 <- crossprod(npWT)
		Ivtheta <- theta%*%Iv
		np2 <- tcrossprod(Ivtheta,theta)

		delta2 <- 1/rgamma(1,(N+naleat+ce+2*dre)/2,(np1+2*dse+np2)/2)

		delta  <- sqrt(delta2)

		### Amostra lambda (GS elemento a elemento)

		lambda <- rgamma(N, (nu+1)/2, (nu+(npWT)^2)/2)
		R  <- diag(lambda)

		### Amostra nu (Metropolis-Hastings)
		nuc <- 2
		while(nuc < 3){
			nuc <- rpois(1,nu)
		}

		# pnu   <- exp(N*((nu/2)*log(nu/2)-lgamma(nu/2))+sum(((nu/2)-1)*log(lambda)+(-(nu/2)*lambda)) - 2*log(1+nu))
		# pnuc  <- exp(N*((nuc/2)*log(nuc/2)-lgamma(nuc/2))+sum(((nuc/2)-1)*log(lambda)+(-(nuc/2)*lambda)) - 2*log(1+nuc))
		pnu   <- N*((nu/2)*(nu/2)-lgamma(nu/2))+sum(((nu/2)-1)*(lambda)+(-(nu/2)*lambda)) - 2*(1+nu)
		pnuc  <- N*((nuc/2)*(nuc/2)-lgamma(nuc/2))+sum(((nuc/2)-1)*(lambda)+(-(nuc/2)*lambda)) - 2*(1+nuc)

		dnu   <- dpois(nu,nu)
		dnuc  <- dpois(nuc,nuc)

		num   <- pnuc*dnu
		denom <- pnu*dnuc

		ifelse((num > denom), alfa <- 1, (alfa <- as.double(num/denom)))
		if(runif(1) < alfa)  nu <- nuc 

		### Amostragem dos parametros de limiar

		nkp   <- 0.8*nk[2:(ce-1)]
		pnovo <- c(rdiric(1,shape=nkp))

		for(j in 2:(ce-1)) taunovo[j] <- sum(pnovo[1:(j-1)])

		lfp     <- nkp*log(p)
		lfpnovo <- nkp*log(pnovo)

		mw  <- tcrossprod(W,theta)

		dp  <- sqrt(delta/lambda)

		probtau <- .C("probtaut", as.integer(Y), ftau=as.double(ftau), ftaunovo=as.double(ftaunovo), 
					  as.integer(escala), as.integer(max(escala)), as.integer(min(escala)), as.integer(ce),
					  as.integer(N), as.double(mw), as.double(tau), as.double(taunovo), as.double(dp),
					  PACKAGE="Bayesthresh")[c("ftau","ftaunovo")]

		ftau     <- probtau$ftau
		ftaunovo <- probtau$ftaunovo

		lftau     <- log(ftau)
		lftaunovo <- log(ftaunovo)

		lftau[is.infinite(lftau)] <- log(1e-15)
		lftaunovo[is.infinite(lftaunovo)] <- log(1e-20)


		alfa      <- exp(sum(lftaunovo-lftau)-sum(lfpnovo-lfp))

		mini <- min(1,alfa)

		if(runif(1) < mini){
			p   <- pnovo
			tau <- taunovo

			### Condicional para L
			### Funcao para a densidade acumulada da variavel Latente

			probacum  <- .C("acumt", as.integer(Y), L=as.double(L), as.integer(escala), as.integer(max(escala)),
							as.integer(min(escala)), as.integer(length(escala)), as.integer(N), as.double(vero), 
							as.double(mw), as.double(tau), as.double(dp), verossim=as.double(verossim),
							PACKAGE="Bayesthresh")[c("L","verossim")]

			L        <- probacum$L
			verossim <- probacum$verossim
			#    if(verossim == 0) verossim <- 1e-320

		}

		### Gera a saida no arquivo cadeia.txt
		if(cont > burn && (cont-burn)%%jump == 0) {
			saida <- c(theta,tau[-(length(tau))],vu,verossim,delta2)
			output <- output + saida
			sumsquare <- sumsquare + (saida^2)
			if(Write==TRUE){
				outp.mcmc[[1]][rowmcmc,] <- theta
				outp.mcmc[[2]][rowmcmc,] <- tau[-(length(tau))]
				outp.mcmc[[3]][rowmcmc,] <- vu
				rowmcmc <- rowmcmc + 1
			}
		}
	}
	list(Sum=output, SumSquare = sumsquare,outp.mcmc=outp.mcmc)
}
### Funcao que retorna a matriz de delineamento dos efeitos aleatorios

MatrixW <- function(FL)
{
	dimen <- NULL
	nL <- FL$dims[1]
	W <- t(as.matrix(FL$trms[[1]]$Zt))
	if(nL > 1){
		for(i in 2:FL$dims[1]){
			Ww <- t(as.matrix(FL$trms[[i]]$Zt))
			W <- cbind(W,Ww)
		}
	}
	return(W)
}


### Promove a saida


outp <- function(formula,res,method,nomesX,nomesZ,ntau,NomesVu, inter)
{
	resi <- NULL
	ans <- NULL
	ifelse(method$algorithm != "NC", lik <- res[nrow(res),], lik <- res[(nrow(res))-1,])
	colnames(lik) <- c("Post. mean","Post.std.dev")
	rownames(lik) <- c("")
	dim1 <- (length(nomesX)+length(nomesZ)+length(ntau)+1)
	dimCompVar <- c(dim1:(dim1+length(NomesVu)))
	compVar <- res[dimCompVar,]
	#     compVar[,2] <- sqrt(compVar[,1] )
	colnames(compVar) <- c("Post.variance", "Post.std.dev")
	rownames(compVar) <- c(NomesVu,"Residuals")
	ifelse(method$algorithm != "NC", resi <- c(1,1), resi <- res[nrow(res),])
	compVar[nrow(compVar),] <- resi
	EfFixef <- res[1:(length(nomesX)),]
	colnames(EfFixef) <- c("Estimate", "Std. Dev")
	rownames(EfFixef) <- nomesX
	ans$compVar <- compVar
	ans$NomesZ <- nomesZ
	ans$algorithm <- method$algorithm
	EfRandom <- res[((length(nomesX))+1):((length(nomesX))+length(nomesZ)),]
	ans$EfRandom <- data.frame(EfRandom)
	Ntau <- res[(dim1-(length(ntau))):(dim1-1),1][-1]
	names(Ntau) <- ntau[-1]
	ans$method <- method$method
	ans$link <- method$link
	ans$Ntau <- Ntau
	ans$EfFixef <- EfFixef
	ans$lik <- lik
	ans$inter <- inter
	ans$formula <- formula
	return(ans)
}


summary.Bayesthresh <- function(object, ...)
{
	if(!inherits(object, "Bayesthresh"))
		stop("Use an object of class Bayesthresh")
	else{
		veros <- as.numeric(object$lik)
		Postmean <- veros[1]
		Poststd <- veros[2]
		vero_dat <- data.frame(Postmean, Poststd)
		colnames(vero_dat) <- c("Post. mean", "Post.std.dev")
		rownames(vero_dat) <- c(" ")
		cat("Threshold model with algorithm",object$algorithm, "and link", object$link, "\n")
		cat("Formula:", deparse(object$formula), "\n")

		cat("\nDeviance:", -2*veros[1], "\n")

		cat("\nMarginal Log-likelihood:\n")
		print(vero_dat)
		#    cat("\nThresholds:\n")
		#    print(object$Ntau)
		cat("\nRandom effects:\n")
		print(object$compVar)
		cat("\nFixed effects:\n")
		print(object$EfFixef)
		cat("\nIteraction Control:\n")
		cat(paste("Burn =",object$inter$burn),",", paste("Jump =",object$inter$jump),",", paste("Iteraction =",object$inter$ef.iter),fill=TRUE)
		cat("Time elapsed", object$final.time, "seconds",  "\n")
	}
}



print.Bayesthresh <- function(x,...)
{
	summary.Bayesthresh(x,...)
}



### The main function
Bayesthresh <-function(formula, data, subset, na.action, A=NULL, algor = list(algorithm = "NC", link = "Gaussian"),
					   Write = FALSE, priors = list(ru=10, su=2, dre=20, dse=5),
					   burn = 50, jump = 2, ef.iter = 4000, model=TRUE)
	### Linear Mixed-Effects in R with threshold models
{
	mc <- match.call()
	stopifnot(length(formula <- as.formula(formula)) == 3)
	fr <- lmerFrames(mc, formula, contrasts) # model frame, X, etc.
	FL <- lmerFactorList(formula, fr, 0L, 0L)
	Zl <- MatrixW(FL)
	X <- fr$X
	Y <- fr$Y
	nomesX <- names(fr$fixef)
	nomesZ <- colnames(Zl)
	escala <- as.numeric(levels(factor(Y)))
	NomesVu <- names(FL$fl)
	ntau <- paste("tau", 1:(length(escala)-1), sep=".")
	nomes <- c(nomesX,nomesZ,ntau,NomesVu, c("verossim"))
	nomes2 <- c(nomesX,nomesZ,NomesVu)
	W <- cbind(as.matrix(X),as.matrix(Zl))
	ntrat <- ncol(as.matrix(X))
	nfixo <- ntrat
	naleat <- ncol(as.matrix(Zl))
	if(!is.null(A)){
		if(dim(A)[1] != dim(A)[2]) stop("non-square matrix A")
	}
	else
		A <- diag(naleat)


		alg <- algorithm(algor$algorithm, algor$link)
		alg <- as.character(call(alg[[1]]))
		FUN <- match.fun(alg)
		prior <- priors.control(priors$ru, priors$su, priors$dre, priors$dse, algor$algorithm)

		initial.time <- proc.time()

		if(algor$algorithm == "NC"){
			result <- FUN(Y,X,Zl,W,A,escala, prior$ru, prior$su,
						  prior$dre, prior$dse, ntrat, nfixo,
						  naleat, burn, jump, ef.iter, Write,FL,nomes,nomes2)
		}
		else
			result <- FUN(Y,X,Zl,W,A,escala, prior$ru, prior$su, ntrat, nfixo,
						  naleat, burn, jump, ef.iter, Write, FL,nomes,nomes2)

		Mean <- result$Sum/ef.iter
		stadev <- sqrt((result$SumSquare - ((result$Sum)^2)/ef.iter)/ef.iter)
		res <- data.frame(Mean,stadev)
		colnames(res) <- c("Post.mean", "Post.std.dev")
		inter <- list(burn=burn, jump=jump, ef.iter=ef.iter)
		saida <- outp(formula,res,algor,nomesX,nomesZ,ntau,NomesVu, inter)
		saida$fl <- unlist(lapply(FL$fl, function(x) length(levels(x))))
		saida$X <- X
		saida$Zl <- Zl
		saida$Y <- Y
		saida$NamesTheta <- c(nomesX,nomesZ)
		saida$NomesVu <- NomesVu
		saida$outp.mcmc <- result$outp.mcmc
		saida$Write <- Write
		#  cat("Threshold model with algorithm", algor$algorithm, "and link", algor$link, "\n")
		#  cat("Formula:", deparse(formula), "\n")
		#  cat(paste("Burn =",inter$burn),",", paste("Jump =",inter$jump),",", paste("Iteraction =",inter$ef.iter),fill=TRUE)
		final.time <- proc.time()-initial.time
		final.time <- final.time[3]
		#  cat("Time elapsed", final.time, "seconds",  "\n")
		saida$final.time <- final.time
		class(saida) <- c("Bayesthresh")
		invisible(saida)
}


