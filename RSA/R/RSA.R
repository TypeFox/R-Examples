#' @title Performs several RSA model tests on a data set with two predictors
#'
#' @description
#' Performs several RSA model tests on a data set with two predictors
#'
#' @details
#' Even if the main variables of the model are normally distirbuted, their squared terms and interaction terms are necessarily non-normal. By default, the RSA function uses a scaled test statistic (\code{test="Satorra-Bentler"}) and robust standard errors (\code{se="robust"}), which are robust against violations of the normality assumption. 
#'
#' \emph{Why does my standard polynomial regression give different p-values and SEs than the RSA package? Shouldn't they be the same?} This is due to the robust standard errors employed in the RSA package. If you set \code{estimator="ML"} and \code{se="standard"}, you get p-values that are very close to the standard approach. (They might still not be identical because the standard regression approach usually uses an OLS estimator and RSA uses an ML estimator).
#'
#' You can also fit \strong{binary outcome variables} with a probit link function. For that purpose, the response variable has to be defined as "ordered", and the \code{lavaan} estimator changed to "WLSMV": \code{r1 <- RSA(Z.binary ~ X*Y, dat, ordered="Z.binary", estimator="WLSMV")} (for more details see the help file of the \code{\link{sem}} function in the \code{lavaan} package.). The results can also be plotted with probabilities on the z axis using the probit link function: \code{plot(r1, link="probit", zlim=c(0, 1), zlab="Probability")}. \code{lavaan} at the moment only supports a probit link function for binary outcomes, not a logit link.
#'
#' @export
#' @import lavaan
#' @import ggplot2
#' @import lattice
#' @import RColorBrewer
#' @param formula A formula in the form \code{z ~ x*y}, specifying the variable names used from the data frame, where z is the name of the response variable, and x and y are the names of the predictor variables.
#' @param data A data frame with the variables
#' @param center Should predictor variables be centered on \emph{each variable's} sample mean before analyses? You should think carefully about this option, as different centering of the predictor variables can affect the commensurability of the predictor scales.
#' @param scale Should predictor variables be scales on the SD of \emph{each variable} before analyses? You should think carefully about this option, as different scaling of the predictor variables can affect the commensurability of the predictor scales.
#' @param na.rm Remove missings before proceeding?
#' @param add Additional syntax that is added to the lavaan model. Can contain, for example, additional constraints, like "p01 == 0; p11 == 0"
#' @param out.rm Should outliers according to Bollen & Jackman (1980) criteria be excluded from the analyses? In large data sets this analysis is the speed bottleneck. If you are sure that no outliers exist, set this option to FALSE for speed improvements.
#' @param breakline Should the breakline in the unconstrained absolute difference model be allowed (the breakline is possible from the model formulation, but empirically rather unrealistic ...). Defaults to \code{FALSE}
#' @param verbose Should additional information during the computation process be printed?
#' @param models A vector with names of all models that should be computed. Should be any from \code{c("absdiff", "absunc", "diff", "mean", "additive", "IA", "SQD", "RR", "SRR", "SRRR", "SSQD", "SRSQD", "full", "null", "onlyx", "onlyy", "onlyx2", "onlyy2")}. For \code{models="all"}, all models are computed, for \code{models="default"} all models besides absolute difference models are computed.
#' @param cubic Should a cubic model with the additional terms Y^3, XY^2, YX^2, and X^3 be included? WARNING: This is experimental, and not all functions will treat the cubic extension properly yet.
#' @param control.variables A string vector with variable names from \code{data}. These variables are added as linear predictors to the model (in order "to control for them"). No interactions with the other variables are modeled. WARNING: This feature is not implemented yet!
#' @param estimator Type of estimator that should be used by lavaan. Defaults to "MLR", which provides robust standard errors, a robust scaled test statistic, and can handle missing values. If you want to reproduce standard OLS estimates, use \code{estimator="ML"} and \code{se="standard"}
#' @param se Type of standard errors. This parameter gets passed through to the \code{sem} function of the \code{lavaan} package. See options there. By default, robust SEs are computed. If you use \code{se="boot"}, \code{lavaan} provides CIs and p-values based on the bootstrapped standard error. If you use \code{confint(..., method="boot")}, in contrast, you get CIs and p-values based on percentile bootstrap (see also \code{\link{confint.RSA}}).
#' @param missing Handling of missing values. By default (\code{NA}), Full Information Maximum Likelihood (FIML) is employed in case of missing values. If cases with missing values should be excluded, use \code{missing = "listwise"}.
#' @param ... Additional parameters passed to the \code{lavaan} \code{\link{sem}} function.
#'
#'
#' @seealso \code{\link{demoRSA}}, \code{\link{plotRSA}}, \code{\link{RSA.ST}}, \code{\link{confint.RSA}}, \code{\link{compare}}
#'
#' @examples
#' # Compute response surface from a fake data set
#' set.seed(0xBEEF)
#' n <- 300
#' err <- 15
#' x <- rnorm(n, 0, 5)
#' y <- rnorm(n, 0, 5)
#' df <- data.frame(x, y)
#' df <- within(df, {
#' 	diff <- x-y
#' 	absdiff <- abs(x-y)
#' 	SD <- (x-y)^2
#' 	z.diff <- diff + rnorm(n, 0, err)
#' 	z.abs <- absdiff + rnorm(n, 0, err)
#' 	z.sq <- SD + rnorm(n, 0, err)
#' 	z.add <- diff + 0.4*x + rnorm(n, 0, err)
#' 	z.complex <- 0.4*x + - 0.2*x*y + + 0.1*x^2 - 0.03*y^2 + rnorm(n, 0, err)
#' })
#' \dontrun{
#' r1 <- RSA(z.sq~x*y, df)
#' summary(r1)
#' compare(r1)
#' plot(r1)
#' plot(r1, model="SRSQD")
#' plot(r1, model="full", type="c")
#' getPar(r1, "coef")	# print model parameters including SE and CI
#' RSA.ST(r1)	# get surface parameters
#'
#' # Motive congruency example
#' data(motcon)
#' r.m <- RSA(postVA~ePow*iPow, motcon)
#'
#' # Get boostrapped CIs with 10 bootstrap samples (usually this should be set to 5000 or higher),
#' # only from the SSQD model
#' c1 <- confint(r.m, model="SSQD", method="boot", R=10)
#' 
#' # Plot the final model
#' plot(r.m, model="RR", xlab="Explicit power motive", 
#' 		ylab="Implicit power motive", zlab="Affective valence")
#' }

RSA <- function(formula, data=NULL, center=FALSE, scale=FALSE, na.rm=FALSE, 
	out.rm=TRUE, breakline=FALSE, models="default", cubic=FALSE, 
	verbose=TRUE, add = "", estimator="MLR",
	se = "robust", missing=NA, ..., control.variables=c()) {


	if (length(control.variables) > 0) stop("Control.variables feature not implemented yet!")

	validmodels <- c("absdiff", "absunc", "diff", "mean", "additive", "IA", "SQD", "SRRR", "SRR", "RR", "SSQD", "SRSQD", "full", "null", "onlyx", "onlyy", "onlyx2", "onlyy2", "weak", "strong")
	
	if (length(models)==1 & models[1]=="all") {models <- validmodels}
	#if (length(models)==1 & models[1]=="default") {models <- c("diff", "mean", "additive", "IA", "SQD", "SRRR", "SRR", "RR", "SSQD", "SRSQD", "full", "null", "onlyx2", "onlyy2", "onlyx", "onlyy", "weak", "strong")}
	if (length(models)==1 & models[1]=="default") {models <- c("additive", "IA", "SQD", "SRRR", "SRR", "RR", "SSQD", "SRSQD", "full", "null", "onlyx2", "onlyy2", "onlyx", "onlyy")}
	if (any(!models %in% validmodels))
		stop("Unknown model name provided in parameter 'models'.")
	

	
	# set all result objects to NULL as default
	s.NULL <- s.full <- s.IA <- s.diff <- s.mean <- s.absdiff <- s.additive <- s.SQD <- s.SSQD <- s.SRSQD <- s.absunc <- s.cubic <- s.RR <- s.SRR <- s.SRRR <- s.onlyx <- s.onlyy <- s.onlyx2 <- s.onlyy2 <- s.weak <- s.strong <-  NULL
	SRSQD.rot <- ""
	SRRR.rot <- ""
	
	add <- paste0("\n# User defined syntax:\n", add)
	
	DV <- all.vars(formula)[1]
	IV1 <- all.vars(formula)[2]
	IV2 <- all.vars(formula)[3]

	## Step 0a: Standardize values (if requested) and calculate higher order terms
	df <- data[, c(DV, IV1, IV2, control.variables)]	# reduce data frame to actually used variables
	df[, IV1] <- scale(df[, IV1], center=center, scale=scale)
	df[, IV2] <- scale(df[, IV2], center=center, scale=scale)
		
	df <- add.variables(formula, data.frame(data.matrix(df)))
	
	# give warnings if the zero point is outside of data range
	if (0 < min(df[, IV1], na.rm=TRUE) | 0 > max(df[, IV1], na.rm=TRUE)) 
		warning(paste("The numerical zero point is outside of the range of variable", IV1, ". Please consider re-centering the variable."))
	if (0 < min(df[, IV2], na.rm=TRUE) | 0 > max(df[, IV2], na.rm=TRUE)) 
		warning(paste("The numerical zero point is outside of the range of variable", IV2, ". Please consider re-centering the variable."))
		
	# give warning if one variable has a much higher range than the other variable
	if ((max(df[, IV1], na.rm=TRUE) - min(df[, IV1], na.rm=TRUE)) / (max(df[, IV2], na.rm=TRUE) - min(df[, IV2], na.rm=TRUE)) > 2)
		warning("Predictor variables have a very different range (by factor 2 or larger)- please check scaling of variables.")
	
	
	
	# set defaults for missing option
	if (is.na(missing)) {
		if (any(is.na(df))) {
			missing <- "fiml"
			warning("There are missing values in your data set. Model is computed with option `missing = 'fiml'`. This is only valid if the data are missing completely at random (MCAR) or missing at random (MAR)! If you want to exclude NA, use `missing = 'listwise'`", call.=FALSE)
		} else {
			missing <- "listwise"
		}
	}
	
	
	
	IV12 <- paste0(IV1, "2")
	IV22 <- paste0(IV2, "2")
	IV13 <- paste0(IV1, "3")
	IV23 <- paste0(IV2, "3")
	IV_IA <- paste0(IV1, "_", IV2)
	IV_IA2 <- paste0(IV1, "_", IV2, "2")
	IV_IA3 <- paste0(IV1, "2", "_", IV2)
	W_IV1 <- paste0("W_", IV1)
	W_IV2 <- paste0("W_", IV2)

	# define control variable
	CV <- ifelse(length(control.variables > 0), paste0(" + ", paste(control.variables, collapse=" + ")), "")

	## Run polynomial regression as a OLS linear model
	addcubic <- ""
	if (cubic==TRUE) addcubic <- paste0(" + ", paste(IV13, IV23, IV_IA2, IV_IA3, sep=" + "))
	f <- paste0(paste0(DV, " ~ ", paste(IV1, IV2, IV12, IV_IA, IV22, sep=" + ")), addcubic, CV)
	lm.full <- lm(f, df, na.action=na.exclude)
	
	# ---------------------------------------------------------------------
	#  Mark outliers and influential cases 
	
	# define the defaults
	if (is.null(out.rm) || (typeof(out.rm) == "logical" && out.rm == TRUE)) {
		out.rm <- "bj1980"
	}
	if ((typeof(out.rm) == "logical" && out.rm == FALSE)) {
		out.rm <- "none"
	}
	out.rm <- match.arg(out.rm, c("bj1980", "robust", "none"))
	df$out <- FALSE

	#...according to Bollen & Jackman, 1980
	if (out.rm == "bj1980") {
		inf <- influence.measures(lm.full)
		df$out <- apply(inf$is.inf[, c("dffit", "cook.d", "hat")], 1, sum) == 3
		n.out <- sum(na.omit(df$out) == TRUE)
		if (verbose==TRUE & n.out>0) {
			warning(paste("Removed", n.out, "multivariate outlier(s) according to Bollen & Jackman (1980) criteria. Outliers are in row(s):", paste(which(df$out == TRUE) , collapse=", ")))
		}
	}
	#...outpro function from WRS
	if (out.rm == "robust") {
		stop("Robust outlier detection not implemented yet.")
		# out.pro <- outpro(df[, c(DV, IV1, IV2)])
# 		n.out <- out.pro$n.out
# 		df$out[out.pro$out.id] <- TRUE
# 		if (verbose==TRUE & n.out>0) {
# 			warning(paste("Removed", n.out, "multivariate outlier(s) using a robust procedure (Wilcox, 2012). Outliers are in row(s):", paste(which(df$out == TRUE) , collapse=", ")))
# 		}
	}
	
	# Rows with outliers have a NA in $out. Should they be kept or removed? I will keep them, otherwise it would not make sense to have fiml estimation, as all NAs are automatically processed as outliers.
	# FIXME: is this a good default choice?
	df$out[is.na(df$out)] <- FALSE


	## Test all models
	
# suppress some types of lavaan warning, which cannot be ruled out analytically ...
withCallingHandlers({	
	
	# Standard full polynomial of second degree
	poly <- paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + b3*", IV12, " + b4*", IV_IA, " + b5*", IV22, CV)
	
	if ("null" %in% models) {
		s.NULL <- sem(paste0(DV, "~ 1 + 0*", IV1, " + 0*", IV2, " + 0*", IV12, " + 0*", IV_IA, " + 0*", IV22, CV), data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	if ("additive" %in% models) {
		if (verbose==TRUE) print("Computing additive model (additive) ...")
		m.additive <-  paste(poly,
			"b3==0",
			"b4==0",
			"b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
		add, sep="\n")
		s.additive <- sem(m.additive, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	if ("onlyx2" %in% models) {
		if (verbose==TRUE) print("Computing x + x^2 model (onlyx2) ...")
		m.onlyx2 <-  paste(poly,
			"b2==0",
			"b4==0",
			"b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
		add, sep="\n")
		s.onlyx2 <- sem(m.onlyx2, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	if ("onlyy2" %in% models) {
		if (verbose==TRUE) print("Computing y + y^2 model (onlyy2) ...")
		m.onlyy2 <-  paste(poly,
			"b1==0",
			"b3==0",
			"b4==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
		add, sep="\n")
		s.onlyy2 <- sem(m.onlyy2, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}

	if ("onlyx" %in% models) {
		if (verbose==TRUE) print("Computing x model (onlyx) ...")
		m.onlyx <-  paste(poly,
			"b2==0",
			"b3==0",
			"b4==0",
			"b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
		add, sep="\n")
		s.onlyx <- sem(m.onlyx, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	if ("onlyy" %in% models) {
		if (verbose==TRUE) print("Computing y model (onlyy) ...")
		m.onlyy <-  paste(poly,
			"b1==0",
			"b3==0",
			"b4==0",
			"b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
		add, sep="\n")
		s.onlyy <- sem(m.onlyy, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	

	if ("diff" %in% models) {
		if (verbose==TRUE) print("Computing difference model (diff) ...")
		m.diff <- paste(poly,
			"b3==0",
			"b4==0",
			"b5==0",
			"b1 == -b2",
			"a1 := b1+b2",
			"a2 := 0",
			"a3 := b1-b2",
			"a4 := 0",
			add, sep="\n")
			s.diff <- sem(m.diff, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
		#summary(s.diff, fit.measures=TRUE)
	}
	
	if ("mean" %in% models) {
		if (verbose==TRUE) print("Computing mean model (mean) ...")
		m.mean <- paste(poly,
			"b3==0",
			"b4==0",
			"b5==0",
			"b1 == b2",
			"a1 := b1+b2",
			"a2 := 0",
			"a3 := b1-b2",
			"a4 := 0",
			add, sep="\n")
			s.mean <- sem(m.mean, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
		#summary(s.mean, fit.measures=TRUE)
	}

	if ("IA" %in% models) {
		if (verbose==TRUE) print("Computing interaction model (IA)...")
		m.IA <- paste(poly,
			"b3==0",
			"b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)",
			"Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"p10 := Y0 - p11*X0",			
			"p20 := Y0 - p21*X0",
			"PA1.curv := b3 - b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 - b4*p21 + b5*(p21^2)",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
		add, sep="\n")
		
			s.IA <- sem(m.IA, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	if ("SQD" %in% models) {
		if (verbose==TRUE) print("Computing squared difference model (SQD) ...")
		m.SQD <- paste(poly,
			"b1==0",
			"b2==0",
			"b3==b5",
			"b3+b4+b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"PA1.curv := b3 - b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 - b4*p21 + b5*(p21^2)",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			add, sep="\n")			
			s.SQD <- sem(m.SQD, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	if ("SSQD" %in% models) {
		if (verbose==TRUE) print("Computing shifted squared difference model (SSQD) ...")
		m.SSQD <- paste(poly,
			"b1==-b2",
			"b3==b5",
			"b3+b4+b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"PA1.curv := b3 - b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 - b4*p21 + b5*(p21^2)",
			"C := b1 / (2*b3)",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			add, sep="\n")			
			s.SSQD <- sem(m.SSQD, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	if (any(models %in% c("RR"))) {
		if (verbose==TRUE) print("Computing rising ridge model (RR) ...")
		m.RR <- paste(poly,
			"b1==b2",
			"b3==b5",
			"b3+b4+b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"meaneffect := b1+b2",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"PA1.curv := b3 - b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 - b4*p21 + b5*(p21^2)",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			
			add, sep="\n")
			s.RR <- sem(m.RR, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	if (any(models %in% c("SRR"))) {
		if (verbose==TRUE) print("Computing shifted rising ridge model (SRR) ...")
		m.SRR <- paste(poly,
			"b3==b5",
			"b3+b4+b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"PA1.curv := b3 - b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 - b4*p21 + b5*(p21^2)",
			"meaneffect := a1",
			"C := (b1-b2) / (4*b3)",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			
			add, sep="\n")
			s.SRR <- sem(m.SRR, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	
	if (any(models %in% c("SRRR"))) {
			if (verbose==TRUE) print("Computing rotated and shifted rising ridge model (SRRR), up ...")
			m.SRRR.up <- paste(paste(poly, " + start(0.01)*", IV12, " + start(0.01)*", IV22),
				"b3 > 0.000001",
				"b5 > 0.000001",
				"b4^2 == 4*b3*b5",
				"a1 := b1+b2",
				"a2 := b3+b4+b5",
				"a3 := b1-b2",
				"a4 := b3-b4+b5",
				"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
				"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
				"PA1.curv := b3 - b4*p11 + b5*(p11^2)",
				"PA2.curv := b3 - b4*p21 + b5*(p21^2)",
				"meaneffect := (b2*b4 - 2*b1*b5) / b4",
				"C := (-2*b1*b5 - b2*b4) / (4*b4*b5)",
				"S := (-b4) / (2*b5)",
				"a4.rescaled := b3/S^2 - b4/S + b5",
				# eigenvalues
				"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
				"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
		
				add, sep="\n")
				s.SRRR.up <- sem(m.SRRR.up, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)	
				
			if (verbose==TRUE) print("Computing rotated and shifted rising ridge model (SRRR), down ...")
			m.SRRR.down <- paste(paste(poly, " + start(-0.01)*", IV12, " + start(-0.01)*", IV22),
			#m.SRRR <- paste(paste(poly, " + start(-0.001)*", IV22),
				"b3 < -0.000001",
				"b5 < -0.000001",
				"b4^2 == 4*b3*b5",
				"a1 := b1+b2",
				"a2 := b3+b4+b5",
				"a3 := b1-b2",
				"a4 := b3-b4+b5",
				"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
				"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
				"PA1.curv := b3 - b4*p11 + b5*(p11^2)",
				"PA2.curv := b3 - b4*p21 + b5*(p21^2)",
				"meaneffect := (b2*b4 - 2*b1*b5) / b4",
				"C := (-2*b1*b5 - b2*b4) / (4*b4*b5)",
				"S := (-b4) / (2*b5)",
				"a4.rescaled := b3/S^2 - b4/S + b5",
				# eigenvalues
				"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
				"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
	
				add, sep="\n")
				s.SRRR.down <- sem(m.SRRR.down, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)	
				
			if (inspect(s.SRRR.up, "converged") == FALSE & inspect(s.SRRR.down, "converged") == TRUE) {
				SRRR.rot <- "down"
			} else 
			if (inspect(s.SRRR.up, "converged") == TRUE & inspect(s.SRRR.down, "converged") == FALSE) {
				SRRR.rot <- "up"
			} else 
			if (inspect(s.SRRR.up, "converged") == TRUE & inspect(s.SRRR.down, "converged") == TRUE) {
				SRRR.rot <- ifelse(fitMeasures(s.SRRR.up, "chisq") > fitMeasures(s.SRRR.down, "chisq"), "down", "up")
			} else {
				if (verbose==TRUE) print("Warning: SRRR model has not converged (neither up nor down curvature)")
			}
			if (SRRR.rot == "up") {
				s.SRRR <- s.SRRR.up
			} else 
			if (SRRR.rot == "down") {
				s.SRRR <- s.SRRR.down
			}
			if (verbose == TRUE) print(paste0("Direction of SRRR curvature: ", SRRR.rot))
			
	}
	
	
	if (any(models %in% c("SRSQD"))) {
		if (verbose==TRUE) print("Computing rotated squared difference model (SRSQD), up ...")
		m.SRSQD.up <- paste(paste(poly, " + start(0.001)*", IV22),
			"b1 == (b2*b4)/(2*b5)",
			"b3 > 0.000001",
			"b5 > 0.000001",
			"b4^2 == 4*b3*b5",	# this is a different (but algebraically equivalent) formulation of the constraints
			
			"C := -.5*(b2/b5)",
			"S := (-b4) / (2*b5)",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"a4.rescaled := b3/S^2 - b4/S + b5",
			"X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)",
			"Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"p10 := Y0 - p11*X0",			
			"p20 := Y0 - p21*X0",
			"PA1.curv := b3 - b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 - b4*p21 + b5*(p21^2)",
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			add, sep="\n")
			s.SRSQD.up <- sem(m.SRSQD.up, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
			
			
		if (verbose==TRUE) print("Computing rotated squared difference model (SRSQD), down ...")
		m.SRSQD.down <- paste(paste(poly, " + start(-0.001)*", IV22),
			"b1 == (b2*b4)/(2*b5)",
			"b3 < -0.000001",
			"b5 < -0.000001",
			"b4^2 == 4*b3*b5",
		
			"C := -.5*(b2/b5)",
			"S := (-b4) / (2*b5)",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"a4.rescaled := b3/S^2 - b4/S + b5",
			"X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)",
			"Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p10 := Y0 - p11*X0",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"p20 := Y0 - p21*X0",
			"PA1.curv := b3 - b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 - b4*p21 + b5*(p21^2)",
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			add, sep="\n")
			s.SRSQD.down <- sem(m.SRSQD.down, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	
			if (inspect(s.SRSQD.up, "converged") == FALSE & inspect(s.SRSQD.down, "converged") == TRUE) {
				SRSQD.rot <- "down"
			} else 
			if (inspect(s.SRSQD.up, "converged") == TRUE & inspect(s.SRSQD.down, "converged") == FALSE) {
				SRSQD.rot <- "up"
			} else 
			if (inspect(s.SRSQD.up, "converged") == TRUE & inspect(s.SRSQD.down, "converged") == TRUE) {
				SRSQD.rot <- ifelse(fitMeasures(s.SRSQD.up, "chisq") > fitMeasures(s.SRSQD.down, "chisq"), "down", "up")
			} else {
				if (verbose==TRUE) warning("Warning: SRSQD model has not converged (neither up nor down curvature)")
			}
			if (SRSQD.rot == "up") {
				s.SRSQD <- s.SRSQD.up
			} else 
			if (SRSQD.rot == "down") {
				s.SRSQD <- s.SRSQD.down
			}
			if (verbose == TRUE) print(paste0("Direction of SRSQD curvature: ", SRSQD.rot))

		
	}
	
	
	if ("full" %in% models) {
		if (verbose==TRUE) print("Computing polynomial model (full) ...")
		m.full <-  paste(poly,
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)",
			"Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p10 := Y0 - p11*X0",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"p20 := Y0 - p21*X0",
			"PA1.curv := b3 - b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 - b4*p21 + b5*(p21^2)",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			# specific tests for fit pattern
			"weakcondition    := b3*b5",				# must be > 0
			"strongcondition1 := (b2*b4)/(2*b5) - b1",	# must not be different from 0
			"strongcondition2 := 2*sqrt(b3*b5)  - b4",	# must not be different from 0
			add,
			sep="\n"
		)
		s.full <- sem(m.full, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}	
	
	if ("weak" %in% models) {
		if (verbose==TRUE) print("Computing weak fit pattern ...")
		m.weak <-  paste(poly,
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)",
			"Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p10 := Y0 - p11*X0",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"p20 := Y0 - p21*X0",
			"PA1.curv := b3 - b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 - b4*p21 + b5*(p21^2)",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			# constraints for weak pattern
			"b3*b5 > 0",	# must be > 0
			add,
			sep="\n"
		)
		s.weak <- sem(m.weak, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	
	if ("strong" %in% models) {
		if (verbose==TRUE) print("Computing strong fit pattern ...")
		m.strong <-  paste(poly,
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"PA1.curv := b3 - b4*p11 + b5*(p11^2)",
			"PA2.curv := b3 - b4*p21 + b5*(p21^2)",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			# constraints for strong pattern
			"b3*b5 > 0.000001",			# must be > 0
			"(b2*b4) == 2*b1*b5",	# must be 0
			"4*b3*b5  == b4^2",	# must be 0
			add,
			sep="\n"
		)
		s.strong <- sem(m.strong, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	if (cubic==TRUE) {
		if (verbose==TRUE) print("Computing full cubic model (cubic) ...")
		m.cubic <-  paste(paste0(poly, " + b9*", IV13, " + b10*", IV_IA2, " + b11*", IV_IA3, " + b12*", IV23),
			"u1 := b1 + b2",				# linear part of LOC
			"u2 := b3 + b4 + b5",			# quadratic part of LOC
			"u3 := b9 + b10 + b11 + b12",	# cubic part of LOC
			"v1 := b1 - b2",				# linear part of LOIC
			"v2 := b3 - b4 + b5",			# quadratic part of LOIC
			"v3 := b9 + b10 - b11 - b12",	# cubic part of LOIC: If v3 != 0, then there is an enhancement effect (i.e., the slope is different on both sides of the optimum)
			add,
			sep="\n"
		)
		#print(m.cubic)
		s.cubic <- sem(m.cubic, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)		
	}
	
	#m.absdiff.JRE <-  paste(
	#	paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + 0*", IV12, " + 0*", IV_IA, " + 0*", IV22, " + 0*W.JRE + b7*W.JRE_", IV1, " + b8*W.JRE_", IV2),
	#	"b1 == -b2",
	#	"b7 == -b8",
	#	"b7 == -2*b1",
	#	add, sep="\n")
	#s.absdiff.JRE <-  sem(m.absdiff.JRE, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	#summary(s.absdiff.JRE, fit.measures=TRUE)
	
	# the unconstrained absolute difference model - Edwards (2002) formula
	#m.absunc.JRE <-  paste(
	#	paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + 0*", IV12, " + 0*", IV_IA, " + 0*", IV22, " + b6*W.JRE + b7*W.JRE_", IV1, " + b8*W.JRE_", IV2),
	#	add, sep="\n")
	#s.absunc.JRE <-  sem(m.absunc.JRE, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	#summary(s.absunc.JRE, fit.measures=TRUE)
	
	
	if ("absdiff" %in% models) {
		if (verbose==TRUE) print("Computing constrained absolute difference model (absdiff) ...")
		m.absdiff <-  paste(
			paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + b6*W + b7*W_", IV1, " + b8*W_", IV2),
			"b1 == 0",
			"b2 == 0",
			"b6 == 0",
			"b7 == -b8",
			add, sep="\n")
			
			s.absdiff <- sem(m.absdiff, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
	if ("absunc" %in% models) {
		# the unconstrained absolute difference model - new formula
		if (verbose==TRUE) print("Computing unconstrained absolute difference model (absunc) ...")
		m.absunc <-  paste(
			paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + b6*W + b7*W_", IV1, " + b8*W_", IV2),
			ifelse(breakline==FALSE, "b6==0", ""),
			add, sep="\n")
			
			s.absunc <- sem(m.absunc, data=df[df$out==FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se=se, estimator=estimator, missing=missing, ...)
	}
	
},	  # end of "withCallingHandlers"

# suppress several types of warning
  warning=function(w) {
	   W <- as.character(w$call)
	   if (
		   (W[1] == "sqrt" & W[2] == "diag(def.cov)" & grepl("NaNs", w$message)) |
		   (W[1] == "sqrt" & W[2] == "b3 * b5") |
		   (W[1] == "nlminb" & W[2] == "x.par") |
		   (W[2] %in% c("m.SRRR.up", "m.SRRR.down", "m.SRSQD.up", "m.SRSQD.down") & grepl("model has NOT converged", w$message))
		  ) {invokeRestart("muffleWarning")}
} )


	# ---------------------------------------------------------------------
	# Sanity check: Check results for convergence problems
	# Sometimes the SRRR or the SRSQD model find a bad solution and have a higher chi2 than their nested models (which is theoretically not possible)
	chisq1 <- plyr::ldply(list(full=s.full, SRRR=s.SRRR, SRR=s.SRR, RR=s.RR, SQD=s.SQD), function(x) {
		chi <- -1
		if (!is.null(x)) {
			if (inspect(x, "converged")==TRUE) chi <-  fitMeasures(x, "chisq")
		}
		return(chi)
	})
	chisq1 <- chisq1[chisq1[, 2]>=0, ]
	if (nrow(chisq1)>1) {
		chisq1$lag <- c(diff(chisq1[, 2], lag=1), NA)
		if (any(chisq1$lag < 0, na.rm=TRUE)) {
			warning(paste0("There are convergence problems with model ", chisq1[which(chisq1$lag < 0), ".id"], ". Its chi-square value is higher than that of a nested model, which is theoretically not possible. Please inspect the results with care, using the compare()-function"))
		}
	}
	
	chisq2 <- plyr::ldply(list(full=s.full, SRRR=s.SRRR, SRSQD=s.SRSQD, SSQD=s.SSQD, SQD=s.SQD), function(x) {
		chi <- -1
		if (!is.null(x)) {
			if (inspect(x, "converged")==TRUE) chi <-  fitMeasures(x, "chisq")
		}
		return(chi)
		
	})
	chisq2 <- chisq2[chisq2[, 2]>=0, ]
	
	if (nrow(chisq1)>1) {
		chisq2$lag <- c(diff(chisq2[, 2], lag=1), NA)
		if (any(chisq2$lag < 0, na.rm=TRUE)) {
			warning(paste0("There are convergence problems with model ", chisq2[which(chisq2$lag < 0), ".id"], ". Its chi-square value is higher than that of a nested model, which is theoretically not possible. Please inspect the results with care, using the compare()-function"))
		}
	}

	
	# ---------------------------------------------------------------------
	# Build results object
	modellist <- list(null=s.NULL, full=s.full, IA=s.IA, diff=s.diff, mean=s.mean, absdiff=s.absdiff, additive=s.additive, SQD=s.SQD, SRRR=s.SRRR, SRR=s.SRR, RR=s.RR, SSQD=s.SSQD, SRSQD=s.SRSQD, absunc=s.absunc, cubic=s.cubic, onlyx=s.onlyx, onlyy=s.onlyy, onlyx2=s.onlyx2, onlyy2=s.onlyy2, weak=s.weak, strong=s.strong)
	
	res <- list(
		models = modellist, 
		SRSQD.rot = SRSQD.rot, SRRR.rot = SRRR.rot, LM=summary(lm.full), formula=formula, 
		data=df, out.rm = out.rm, outliers = which(df$out == TRUE), DV=DV, IV1=IV1, IV2=IV2, IV12=IV12, IV22=IV22, 
		IV_IA=IV_IA, W_IV1=W_IV1, W_IV2=W_IV2, IV13=IV13, IV23=IV23, IV_IA2=IV_IA2, IV_IA3=IV_IA3, 
		r.squared = summary(lm.full)$r.squared)
	
	attr(res, "class") <- "RSA"
	return(res)
}
