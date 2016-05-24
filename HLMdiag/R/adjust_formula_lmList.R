# Adjusting formula for lmList
#
# When dealing with nested data, sometimes a factor will only take on one level
# within a group. This causes an error, since a factor must have 2 or more levels
# to be included as a predictor within a formula for a linear model. When a
# factor takes on only one level, the effect of that factor can then be included
# in the intercept of the model. This function takes a formula for use with \code{lmList}
# and will make the necessary adjustments to avoid errors while fitting the separate LS models. 
#
# @param formula a linear formula that is used by \code{lmList()} e.g. y ~ x1 + x2 + ... + xn | g
# @param data the model frame from the model fit by \code{lmer()}
# @return a list containing the adjusted formulas
# @author Adam Loy \email{loyad01@@gmail.com}
# @keywords models regression
problem_factor_groups <- function(formula, data){
	form <- formula(formula)
	model_frame <- data
	
	g <- deparse(form[[3]][[3]])
	ngroups <- length(unique(data[,g]))
	y <- deparse(form[[2]])
	xs <- setdiff(names(model_frame), c(deparse(form[[2]]), g))
	xs_frame <- model_frame[, xs]
	
	are_factors <- rapply(xs_frame, class)
	are_factors <- names(are_factors)[which(are_factors == "factor")]
	
	id_vars <- c(g, are_factors)
	molten_model_frame <- melt(model_frame, id = id_vars)
	
	one_level_factor <- matrix(NA, nrow = ngroups, ncol = 1 + length(are_factors))
	colnames(one_level_factor) <- c(g, are_factors)
	
	for(i in 1:length(are_factors)){
		cast_model_frame <- dcast(molten_model_frame, paste(g, "~", are_factors[i]), length)
		if(i == 1) one_level_factor[,1] <- cast_model_frame[,1]
		one_level_factor[,(i + 1)] <- apply(cast_model_frame[,-1], 1, function(x){
			empty_factors <- sum(x == 0)
			empty_factors > 0
		})
	}
	
	return(one_level_factor)
}


subbars <- function(term)
#@ Substitute the '+' function for the '|' function (from lmer)
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


adjust_formula_lmList <- function(formula, data){
	lm_form <- formula(formula)
	lm_form[[3]] <- lm_form[[3]][[2]]
	
	mf <- match.call()
	m <- match(c("data", "subset", "weights", "na.action", "offset"), names(mf), 0)
	mf <- mf[c(1, m)]
	mf$formula <- subbars(formula)
	mf$drop.unused.levels <- TRUE
	mf[[1]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	
	problem_mat <- data.frame(problem_factor_groups(formula, mf))
	
	problem_list <- split(problem_mat, problem_mat[,1])

	formula_list <- lapply(problem_list, function(x){
	#formula_list <- vector("list", length(problem_list))
	#for(i in 1:length(problem_list)){
	#	x <- problem_list[[i]]
		changes <- x[,-1]
		
		if(dim(x)[2] == 2){
			if(changes == 1){changes <- colnames(x)[2]
				changes <- paste(changes, ":.", sep = "")
			}
			else changes <- NULL
		} 
		
		else{
			nchanges <- sum(changes)
			
			if(nchanges == 0) changes <- NULL
			
			else{
				if(nchanges == 1){
					changes <- colnames(changes)[which(changes == 1)]
					changes <- paste(changes, ":.", sep = "")
				}
				
				else{
				#if(nchanges > 1){
					changes <- colnames(changes)[which(changes == 1)]
					changes <- paste(changes, ":.", sep = "")
					changes <- paste(changes, collapse = "-")
				}
			}
		}
		
		if(is.null(changes) != TRUE) {change_form <- update(lm_form, paste(". ~ .", changes, sep = "-"))}
		else{ change_form <- lm_form }
		#formula_list[[i]] <- change_form
	})
	
	return(formula_list)
}


pooledSD <- function(x, ...)
{
    stopifnot(is(x, "adjust_lmList"))
    sumsqr <- apply(sapply(x,
                           function(el) {
                               if (is.null(el)) {
                                   c(0,0)
                               } else {
                                   res <- resid(el)
                                   c(sum(res^2), length(res) - length(coef(el)))
                               }
                           }), 1, sum)
    if (sumsqr[2] == 0) {
        stop("No degrees of freedom for estimating std. dev.")
    }
    val <- sqrt(sumsqr[1]/sumsqr[2])
    attr(val, "df") <- sumsqr[2]
    val
}

#' @export
adjust_lmList <- function(object, data, pool){
  UseMethod("adjust_lmList", object)
}


#'Fitting Common Models via \code{lm}
#'
#'Separate linear models are fit via \code{lm} similar to \code{lmList},
#'however, \code{adjust_lmList} can handle models where a factor takes only one
#'level within a group. In this case, the \code{formula} is updated eliminating
#'the offending factors from the formula for that group as the effect is
#'absorbed into the intercept.
#'
#'
#'@aliases adjust_lmList adjust_lmList,formula,data.frame-method
#'@keywords models regression
#'@param formula a linear formula such as that used by \code{lmList}, e.g.
#'\code{y ~ x1 + ... + xn | g}, where \code{g} is a grouping factor.
#'@param data a data frame containing the variables in the model.
#'@param pool a logical value that indicates whether the pooled standard
#'deviation/error should be used.
#'@seealso \code{\link[lme4]{lmList}, \link[stats]{lm}}
#'@references Douglas Bates, Martin Maechler and Ben Bolker (2012). lme4:
#'Linear mixed-effects models using S4 classes. R package version 0.999999-0.
#'@examples
#'
#' data(Exam, package = 'mlmRev')
#' sepLM <- adjust_lmList(normexam ~ standLRT + sex + schgend | school, data = Exam)
#' confint(sepLM)
#'
adjust_lmList.formula <- function(object, data, pool){
# 	options(show.error.messages = FALSE)
# 	lmList_result <- try(lmList(formula = formula, data = data), silent = TRUE)
  lmList_result <- lme4::lmList(formula = object, data = data)
#   options(show.error.messages = TRUE)
	orig_names <- names(lmList_result)
	check_results <- unlist(lapply(lmList_result, is.null))
	form <- formula(object)
	g <- deparse(form[[3]][[3]]) 
	ngroups <- length(unique(data[,g]))
	
	if(sum(check_results) != 0){ #return(lmList_result)
	
	#else{
		new_formulas <- adjust_formula_lmList(object, data)
		problem_cases <- as.numeric(which(check_results == TRUE))
		#print(problem_cases)
		split_data <- split(data, data[,g])
		
		for(i in problem_cases){
			lmList_result[[i]] <- lm(formula = new_formulas[[i]], data = split_data[[i]])
		}
		
	}
	class(lmList_result) <- "adjust_lmList"
	if(missing(pool)) pool <- TRUE
	#lmList_result <- new("adjust_lmList", lmList_result, call = match.call(), pool = pool)
	attr(lmList_result, "names") <- orig_names
  attr(lmList_result,"call") <- match.call()
	lmList_result
}


# setClass("adjust_lmList", representation(call = "call", pool = "logical"), contains = "list")
# 
# setClass("adjust_lmList.confint", contains = "array")


# #' @export
# setMethod("adjust_lmList", signature(formula = "formula", data = "data.frame"),
# 	function(formula, data, pool){
# 	  #   options(show.error.messages = FALSE)
# 	  # 	lmList_result <- try(lmList(formula = formula, data = data), silent = TRUE)
# 	  lmList_result <- lmList(formula = formula, data = data)
# 	  #   options(show.error.messages = TRUE)
# 	  orig_names <- names(lmList_result)
# 		check_results <- unlist(lapply(lmList_result, is.null))
# 		form <- formula(formula)
# 		g <- deparse(form[[3]][[3]]) 
# 		ngroups <- length(unique(data[,g]))
# 	
# 		if(sum(check_results) != 0){ #return(lmList_result)
# 	
# 		#else{
# 			new_formulas <- adjust_formula_lmList(formula, data)
# 			problem_cases <- as.numeric(which(check_results == TRUE))
# 			#print(problem_cases)
# 			split_data <- split(data, data[,g])
# 			
# 			for(i in problem_cases){
# 				lmList_result[[i]] <- lm(formula = new_formulas[[i]], data = split_data[[i]])
# 			}
# 		
# 		}
# 		#class(lmList_result) <- "adjust_lmList"
# 		if(missing(pool)) pool <- TRUE
# 		lmList_result <- new("adjust_lmList", lmList_result, call = match.call(), pool = pool)
# 		attr(lmList_result, "names") <- orig_names
# 		lmList_result
# 	})

#' @export
#' @method coef adjust_lmList
#' @S3method coef adjust_lmList
# setMethod("coef", signature(object = "adjust_lmList"),
coef.adjust_lmList <- 	function(object, ...){
			coefs <- lapply(object, coef)
			non.null <- !unlist(lapply(coefs, is.null))
			if(sum(non.null) > 0){
				coefs <- lapply(coefs, function(x) as.data.frame(t(x)))
				coefs <- do.call('rbind.fill', coefs)
				rownames(coefs) <- names(object)
				#dimnames(coefs) <- list(names(object), names(coefs))
			#coefs <- as.data.frame(coefs)
			#effectNames <- names(coefs)
			}
			coefs
	}
# )


#' @export
#' @method print adjust_lmList
#' @S3method print adjust_lmList
# setMethod("show", signature(object = "adjust_lmList"), 
print.adjust_lmList <- 	function(x, ...){
		cat("Call:", deparse(attr(x, "call")), "\n")
		cat("Coefficients:\n")
		invisible(print(coef(x)))
	}
# )


#' @export
#' @method confint adjust_lmList
#' @S3method confint adjust_lmList
# setMethod("confint", signature(object = "adjust_lmList"),
confint.adjust_lmList <- function(object, parm, level = 0.95, pool = NULL, ...){
		#if(length(object < 1))
		#	return(new("adjust_lmList.confint", array(numeric(0), c(0,0,0))))
		coefs <- coef(object)
		val <- array(NA, c(dim(coefs)[2], 2, length(object)) )
		dimnames(val) <- list(names(coefs), 
			c(paste((1-level)/2 * 100, "%"), paste((1 - (1-level)/2) * 100, "%")), 
			names(object))
		#val <- aperm(val, c(3, 2, 1))l
		if(is.null(pool)) pool <- object@pool
		if(length(pool) > 0 && pool[1]){
			sd <- pooledSD(object)
			a <- (1 - level)/2
			fac <- sd * qt(c(a, 1 - a), attr(sd, "df"))
			parm <- names(coefs)
			for (i in seq_along(object)){
				cis <- coef(object[[i]]) + 
					sqrt(diag(summary(object[[i]], corr = FALSE)$cov.unscaled)) %o% fac
				for(j in names(coef(object[[i]]))){
					val[j, , i] <- cis[j,]
				}
			}
		}
		else{
			cis <- lapply(object, confint)
			for (i in seq_along(object)) {
				for(j in names(coef(object[[i]]))){
					val[j, , i] <- cis[[i]][j,]
				}
			}
		}
		class(val) <- "adjust_lmList.confint"
    val
	} #, valueClass = "adjust_lmList.confint"
# )


#' @export
#' @method plot adjust_lmList.confint
#' @S3method plot adjust_lmList.confint
# setMethod("plot", signature(x = "adjust_lmList.confint"),
plot.adjust_lmList.confint <-	function(x, y, ...){
# 		stopifnot(require("ggplot2"))
    group <- intervals <- NULL # Make codetools happy
		cis <- as(x, "array")
		df <- adply(cis, c(1, 2, 3), identity)
		colnames(df) <- c("what", "end", "group", "intervals")
		p <- ggplot(df, aes(x = group, y = intervals))
		p + geom_line() +
			geom_errorbar(aes(ymin = intervals, ymax = intervals)) + 
			facet_wrap( ~ what, nrow = 1, scales = "free")
# 			coord_flip() + ylab(NULL)
	}
# )



# setMethod("formula", signature(x = "adjust_lmList"),
#           function(x, ...) x@call[["formula"]])

# To ensure compatibility across the lme4 development versions
# modelFormula <- function(form)
# {
#   if (class(form) != "formula" || length(form) != 3)
#     stop("formula must be a two-sided formula object")
#   rhs <- form[[3]]
#   if (class(rhs) != "call" || rhs[[1]] != as.symbol('|'))
#     stop("rhs of formula must be a conditioning expression")
#   form[[3]] <- rhs[[2]]
#   list(model = form, groups = rhs[[3]])
# }

# To ensure compatibility across the lme4 development versions
# setMethod("lmList", signature(formula = "formula", data = "data.frame"),
#           function(formula, data, family, subset, weights,
#                    na.action, offset, pool, ...)

# lmList <- function (formula, data, family, subset, weights, na.action, offset, pool, ...) {
#             mCall <- mf <- match.call()           
#             m <- match(c("family", "data", "subset", "weights",
#                          "na.action", "offset"), names(mf), 0)
#             mf <- mf[c(1, m)]
#             ## substitute `+' for `|' in the formula
#             mf$formula <- subbars(formula) 
#             mf$x <- mf$model <- mf$y <- mf$family <- NULL
#             mf$drop.unused.levels <- TRUE
#             mf[[1]] <- as.name("model.frame")
#             frm <- eval(mf, parent.frame())
#             mform <- modelFormula(formula)
#             if (missing(family)) {
#               val <- lapply(split(frm, eval(mform$groups, frm)),
#                             function(dat, formula)
#                             {
#                               ans <- try({
#                                 data <- as.data.frame(dat)
#                                 lm(formula, data)
#                               })
#                               if (inherits(ans, "try-error"))
#                                 NULL
#                               else ans
#                             }, formula = mform$model)
#             } else {
#               val <- lapply(split(frm, eval(mform$groups, frm)),
#                             function(dat, formula, family)
#                             {
#                               ans <- try({
#                                 data <- as.data.frame(dat)
#                                 glm(formula, family, data)
#                               })
#                               if (inherits(ans, "try-error"))
#                                 NULL
#                               else ans
#                             }, formula = mform$model, family = family)
#             }
#             if (missing(pool)) pool <- TRUE
#             new("lmList", val, call = mCall, pool = pool)
#           }
