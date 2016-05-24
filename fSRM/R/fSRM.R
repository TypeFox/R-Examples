#' @title Run a Social Relations Model with roles ("Family SRM")
#' @aliases fSRM
#' @aliases print.fSRM
#'
#' @description
#' Run a Social Relations Model with roles ("Family SRM")
#'
#' @details
#' A model can be rerun with additional syntax using the \code{add} function:
#' \code{s1 <- fSRM(dep1/dep2 ~ actor*partner | fam, dat2)}
#' \code{s2 <- add(s1, "Ac ~~ Pm")}
#' A model can be run with new parameters using the update function:
#' \code{s1 <- fSRM(dep1/dep2 ~ actor*partner | fam, dat2)}
#' \code{s2 <- update(s1, diff=TRUE)}
#'
#' @export
#' @import lavaan
#' @import plyr
#' @import reshape2
#' @param formula A formula that defines the variable names. Should be in one of following formats: (1) Single manifest dependent variable: DV ~ actor.id * partner.id | group.id, (2) Multiple indicators for dependent variable: DV1/DV2/DV3 ~ actor.id * parter.id | group.id.
#' @param data A data frame with all variables defined by \code{formula}. Must be in long format where each row describes one directed dyadic relationship.
#' @param drop In three-member families at least one component has to be dropped. \code{drop} defines which one: "none": drop nothing; "family" - drop family effect; "GR" - drop generalized reciprocities; "actor" - drop actor factors and actor-partner covariances; "partner" - drop partner effects and actor-partner covariances; "default": drop nothing in >= 4 members and drop family effect with 3 members. Although usually not necessary, the drop parameter can also be applied to >= 4 member families.
#' @param add Additional lavaan syntax pasted at the end of the generated model. Can contain, for example, user specified error correlations.
#' @param IGSIM Define intragenerational similarity correlations. Must be a list where the levels of actor.id and partner.id are combined, e.g.: \code{IGSIM=list(c("m", "f"), c("c", "y"))}. Here "m"other and "f"ather are defined as one generation, and "y"ounger and "o"lder as the other generation.
#' @param syntax In that variable the user can directly provide a lavaan model syntax. Then no automatical model syntax is generated.
#' @param add.variable Not yet fully implemented: Add external variables to the model syntax.
#' @param ... Additional arguments passed to the \code{sem} function of \code{lavaan}
#' @param means Should the structured means of the SRM factors be calculated?
#' @param pairwise Compute pairwise comparison of actor and partner means between all roles? Only works when \code{means} is also set to TRUE.
#' @param group Variable name indicating group membership. This can be used two compare two classes of families (e.g., clinical families vs. control families). If this variable is provided, it must contain exactly two levels.
#' @param diff Compare groups with the delta method? You need to specify a group identifier in parameter \code{group}. If \code{diff = TRUE} and \code{means = FALSE}, only variances are compared between groups. If \code{diff = TRUE} and \code{means = TRUE}, variances and means are compared between groups.
#' @param noNegVar Should variance estimates be constrained to be positive?
#' @param missing Handling of missing values. By default (\code{NA}), Full Information Maximum Likelihood (FIML) is employed in case of missing values. If families with missing values should be excluded, use \code{missing = "listwise"}
#' @param rolesEqual Maximal constraints: Do roles matter at all? If this parameter is set to TRUE, it is a model with no mean difference, the actor variances equal, partner variances equal, relationship variances equal, and the respective reciprocities equal (Thanks to a suggestion of David Kenny). Model comparisons via \code{anova} can show whether roles matter at all.

## OLD PARAMETERS, NOT CURRENTLY USED
# @param err Defines the type of correlations between error terms. err = 1: Correlate same items BETWEEN ALL RATERS (e.g., Dyadic Data Analysis, Kenny, Kashy, & Cook, 2000); err = 2: Correlate same items WITHIN RATERS (e.g., Branje et al., 2003, Eichelsheim)
# @param self Should self-ratings be included in the analysis (if present in the data set)?
# @param selfmode Defines the style how the selfratings are combined with the latent actor and partner effects. If \code{selfmode="cor"} they are correlated (as in REFERENCE), if \code{selfmode="kq"} the k and q paths are calculated (see Kenny & West, 2010)

#' @details
#' The \code{fSRM} function relies on the \code{lavaan} package for computation: A syntax for the SRM with roles is generated and then passed to the \code{lavaan} function. Hence, many options of the \code{lavaan} function can be used out-of-the-box (additional parameters are passed to the \code{lavaan} function through the \code{...} operator). For example, one can deal with missing values. The default behavior is to exclude families with missing values (listwise deletion). Set \code{fSRM(..., missing="fiml")} for ML / FIML estimation. Or, you can request bootstrapped standard errors with \code{fSRM(..., se="boot")}.

#' You can test for a very restricted model by constraining the roles to be equal ("Do roles matter at all?"). Therefore, compare a model with free roles (\code{m1 <- fSRM(..., means=TRUE, rolesEqual = FALSE)}) with a model with equal roles (\code{m2 <- fSRM(..., means=TRUE, rolesEqual=TRUE)}) using \code{anova(m1$fit, m2$fit)} (Thanks to David Kenny for the suggestion).

#' For plotting relative variances and mean structure, see \code{\link{plot.fSRM}}.

#' @references
#' Kenny, D. A., & West, T. V. (2010). Similarity and Agreement in Self-and Other Perception: A Meta-Analysis. Personality and Social Psychology Review, 14(2), 196-213. doi:10.1177/1088868309353414

#' @seealso \code{\link{plot.fSRM}}, \code{\link{equalMeans}}, \code{\link{mod}}

#' @examples
#' # Example from Dyadic Data Analysis
#' data(two.indicators)

#' # 4 persons, 1 indicator
#' f4.1 <- fSRM(dep1 ~ actor.id*partner.id | family.id, two.indicators)
#' f4.1		# by default, one-sided p-values and CIs are printed for variances
#' print(f4.1, var.onesided=FALSE)	# Show two-sided p-values and CIs for variances
#' plot(f4.1)	# plot relative variances
#' plot(f4.1, bw=TRUE)
#' 
#' \dontrun{
#' # 4 persons, 2 indicators
#' f4.2 <- fSRM(dep1/dep2 ~ actor.id*partner.id | family.id, two.indicators)
#' f4.2
#' plot(f4.2, bw=TRUE)
#' plot(f4.2, bw=TRUE, onlyStable=TRUE)
#' 
#' # 4 persons, 1 indicator, mean structure
#' f4.1.m <- fSRM(dep1 ~ actor.id*partner.id | family.id, two.indicators, means=TRUE)
#' f4.1.m
#' plot(f4.1.m, means=TRUE)		# plot mean structure
#' 
#' # 4 persons, 2 indicators, mean structure
#' f4.2.m <- fSRM(dep1/dep2 ~ actor.id*partner.id | family.id, two.indicators, means=TRUE)
#' f4.2.m
#' 
#' # ... add intragenerational similarity (now results are identical to Cook, 2000)
#' f4.ig <- fSRM(dep1/dep2 ~ actor.id*partner.id | family.id, two.indicators, 
#'               IGSIM=list(c("m", "f"), c("c", "y")))
#' f4.ig
#' 
#' ## ======================================================================
#' ## Wald-test for equality of means
#' ## ======================================================================
#' 
#' f4.1.m <- fSRM(dep1 ~ actor.id*partner.id | family.id, two.indicators, means=TRUE)
#' f4.1.m
#' equalMeans(f4.1.m)
#' 
#' # construct a 3-person data set
#' two.indicators3 <- 
#' two.indicators[two.indicators$actor.id != "y" & two.indicators$partner.id != "y", ]
#' f3.2.m <- fSRM(dep1/dep2 ~ actor.id*partner.id | family.id, two.indicators3, means=TRUE)
#' f3.2.m
#' equalMeans(f3.2.m)
#' 
#' 
#' 
#' # ---------------------------------------------------------------------
#' # Run analyses based on Eichelsheim, V. I., Buist, K. L., Dekovic, M., 
#'   Cook, W. L., Manders, W., Branje, S. J. T., Frijns, T., et al. (2011). 
#'   Negativity in problematic and nonproblematic families: 
#'   A multigroup social relations model analysis with structured means. 
#'   Journal of Family Psychology, 25, 152-156. doi:10.1037/a0022450
#' 
#' # The data set is a simulated data set which has 
#' # comparable properties as the original data set
#' 
#' data(two.groups)
#' str(two.groups)
#' 
#' E1 <- fSRM(neg ~ actor.id*partner.id | family.id, data=two.groups)
#' E1
#' 
#' # make group comparison: 
#' # group = 1: non-problematic families, group = 2: problematic families
#' # The data set must contain exactly two groups, otherwise an error is printed
#' E2 <- fSRM(neg ~ actor.id*partner.id | family.id, data=two.groups, group="group")
#' E2
#' 
#' # Compare means and differences between groups. Beware: This model takes *really* long ...
#' E3 <- fSRM(neg ~ actor.id*partner.id | family.id, 
#'            data=two.groups, group="group", means=TRUE, diff=TRUE)
#' E3
#' 
#' }

fSRM <-
function(formula=NULL, data, drop="default", add="", means=FALSE, pairwise=FALSE, diff=FALSE, IGSIM=list(), add.variable=c(), syntax="", group=NULL, noNegVar=FALSE, rolesEqual=FALSE, missing=NA, ...) {
	
	dots <- list(...)
	
	if (rolesEqual==TRUE & means==FALSE)
		stop("For the parameter `rolesEqual == TRUE` you also have to set `means=TRUE`.")
	
	# TODO: Re-introduce self-ratings? Preliminarily, fix it to FALSE
	self <- FALSE
	
	# save the function call for use in refitting
	call <- match.call(expand.dots = TRUE)
	
	# Parse the formula
	lhs <- strsplit(gsub(" ","",as.character(formula)[2], fixed=TRUE), "/", fixed=TRUE)[[1]]
	rhs <- strsplit(gsub(" ","",as.character(formula)[3], fixed=TRUE),"\\*|\\|", perl=TRUE)[[1]]
	
	# Retrieve the variable names from the formula
	var.id <- lhs
	actor.id <- rhs[1]
	partner.id <- rhs[2]
	group.id <- NULL
	if (length(rhs)>=3) {
		group.id <- rhs[3]
	} else {
		stop("For SRM with roles a group id has to be defined in the formula (after the | operator).")
	}
	
	# Restructure data format from long to wide
	
	fam0 <- list()
	for (v in c(var.id, add.variable)) {
		fam0[[v]] <- dcast(data[, c(var.id, actor.id, partner.id, group.id, add.variable)], formula(paste(group.id, "~", actor.id, "+", partner.id)), value.var=v)
		colnames(fam0[[v]])[-1] <- paste(colnames(fam0[[v]])[-1], v, sep="_")
	}
	
	fam <- merge.rec(fam0, by=group.id)
	
	# add group variable (for group comparison)
	if (!is.null(group)) {
		# check whether exactly two groups are present in the data set
		if (length(unique(data[, group])) != 2) stop("For group comparisons, the group variable must contain exactly two groups!")
			
		g2 <- ddply(data, group.id, function(x) x[1, group])
		colnames(g2) <- c(group.id, group)
		fam <- merge(fam, g2, by=group.id)
	}
	
	# remove all-NA columns
	NAcol <- which(apply(fam, 2, function(x) sum(is.na(x))) == nrow(fam))
	if(length(NAcol)>0)	{fam <- fam[, -NAcol]}
	#fam <- na.omit(fam)
	included <- fam[, group.id]
	
	# set defaults for missing option
	if (is.na(missing)) {
		if (any(is.na(fam))) {
			missing <- "fiml"
			warning("There are missing values in your data set. Model is computed with option `missing = 'fiml'`. This is only valid if the data are missing completely at random (MCAR) or missing at random (MAR)! If you want to exclude families with missing values, use `missing = 'listwise'`", call.=FALSE)
		} else {
			missing <- "listwise"
		}
	}
	
	roles <- sort(unique(data[, actor.id]))
	
	# define defaults for drop
	drop <- match.arg(drop, c("nothing", "family", "GR", "actor", "partner", "default"))
	if (drop == "default" & length(roles) == 3 & syntax=="") {
		message("Three-member families: Dropping family variance per default.")
		drop <- "family"
	}
	if (drop == "default" & length(roles) > 3 & syntax=="") {drop <- "nothing"}
	
	# Do some sanity checks
	if (length(roles) == 3 & drop == "nothing" & means == FALSE) {warning('Data set with 3-member-groups detected - model is not identified. Maybe you should remove the family effect (drop = "family") or one of the reciprocities?', call.=FALSE)}
	if (!identical(sort(unique(data[, actor.id])), sort(unique(data[, partner.id])))) {
		warning("Actor.id and Partner.id have different factor levels; results might be wrong!", call.=FALSE)
	}
	if (diff==TRUE & is.null(group)) stop("For comparing groups with the delta method you have to provide a `group` variable.")
	
	
	if (!is.null(group)) {
		groupnames <- as.character(unique(fam$group))
	} else {
		groupnames <- NULL
	}
	
	# if no syntax is directly provided:
	if (syntax == "") {
		syntax0 <- buildSRMSyntax(roles, var.id, drop=drop, err="default", IGSIM=IGSIM, means=means, pairwise=pairwise, diff=diff, groupnames=groupnames, self=self, add.variable=add.variable, noNegVar=noNegVar, rolesEqual=rolesEqual)
	
		syntax <- paste(syntax0, add, sep="\n")
	} else {
		print("Model syntax is directly specified; skipping buildfSRMSyntax")
	}
	
	#cat(syntax)
	
	# suppress some types of lavaan warning
	withCallingHandlers({	
		m <- lavaan(
				model		= syntax, 
				data		= fam,
				std.ov		= FALSE,
				orthogonal	= TRUE,
				fixed.x 	= FALSE,
				int.ov.free	= TRUE,
				int.lv.free = FALSE,
				auto.fix.first = FALSE,
				auto.fix.single = TRUE,
				auto.var 	= TRUE,
				auto.cov.lv.x = TRUE,
				auto.cov.y 	= TRUE, 
				group		= group, 
				missing		= missing, ...)
		},	  # end of "withCallingHandlers"

		# suppress two types of warning
		  warning=function(w) {
		   W <- as.character(w)
		   if (
			   grepl("some estimated variances are negative", w$message) |
			   grepl("covariance matrix of latent variables is not positive definite", w$message)
 			  ) {invokeRestart("muffleWarning")}
	})
		
	# does not work for boostrapped SEs? But we don't need it anyway ...
	# suppressWarnings(
	# 	SS <- standardizedSolution(m, type="std.all")
	# )
	
	res <- list(
		fit		= m,
		#SS		= SS,
		syntax	= syntax,
		roles	= roles,
		actor.id 	= actor.id,
		partner.id 	= partner.id,
		group.id 	= group.id,
		var.id	= var.id,
		drop	= drop,
		means	= means,
		pairwise = pairwise,
		diff	= diff,
		group	= group,
		groupnames = groupnames,
		latent	= ifelse(length(var.id) == 1, FALSE, TRUE),
		IGSIM	= IGSIM,
		self	= self,
		call	= call,
		data	= fam)
	
	attr(res, "class") <- "fSRM"

	return(res)
}
