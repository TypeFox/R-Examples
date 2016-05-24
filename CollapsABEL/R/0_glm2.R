#' GLM with arbitrary column names
#' 
#' Substitute column names that are unsuitable for formulas
#' and substitute back when returning results.
#' 
#' @param dat data.frame. Souce data to build GLM upon.
#' @param y character. Column name of dependent variable.
#' @param xs character. Column names of independent variable.
#' @param ... passed to glm.
#' @return data.frame of coefficients.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
glm2 = function(dat, y, xs, ...) {
	name_orig = colnames(dat)
	name_legit = paste("x", seq_along(dat), sep = "")
	names_dat = data.frame(name_orig, name_legit, stringsAsFactors = FALSE)
	dat = setNames(dat, name_legit)
	y1 = changeByMap(y, names_dat)
	xs1 = changeByMap(xs, names_dat)
	tryCatch({
				glm_model = glm(
						as.formula(
								sprintf("%s~%s", y1, strConcat(xs1, "+"))),
						data = dat, 
						...
				)
				coefs = summary(glm_model)$coefficients
				rownames(coefs) = changeByMap(rownames(coefs), 
						names_dat, reverse = TRUE)
				idx = which(! xs %in% rownames(coefs))
				if(length(idx) == 0) {
					return(coefs)
				}
				row_names_to_add = xs[idx]
				new_rownames = c(rownames(coefs), row_names_to_add) 
				coefs = rbind(coefs, 
						matrix(rep(NA, 4 * length(row_names_to_add)), nrow = 4))
				rownames(coefs) = new_rownames
				coefs
			}, error = function(e) {
				n_row = 1 + length(xs)
				n_col = 4
				res = matrix(NA, n_row, n_col)
				colnames(res) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
				rownames(res) = c("(Intercept)", xs)
				res
			})
}

#' Perform glm iteratively through a number of independent variables with fixed dependent variables and covariates.
#' 
#' @param dat data.frame
#' @param y character. Name of dependent variable columns.
#' @param xs character. Names of independent variable columns. 
#' @param covars character. Names of covriate columns.
#' @param ... passed to glm.
#' @return matrix of coefficients
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
glmIter = function(dat, y, xs = NULL, covars = character(), ...) {
	cn = colnames(dat)
	stopifnot(y %in% cn)
	if(length(covars) != 0) {
		stopifnot(all(covars %in% cn))
	}
	if(!is.null(xs)) {
		stopifnot(all(xs %in% cn))
	} else {
		xs = cn[!(cn %in% y) & !(cn %in% covars)]
	}
	t(sapply(xs, function(x) {
				glm2(dat, y, c(covars, x), ...)[x, ]
			}))
}

#' Transform a vector by a mapping
#' 
#' The mapping is represented by a data.frame:
#' 1st column is the domain, 2st column is the range.
#' 
#'
#' @param old_vector vector of any type.
#' @param mapping_dat data.frame, first column must be the same type as the \code{old_vector}
#' @param reverse logical. Reverse domain and range if set to TRUE
#' @examples 
#' \dontrun{
#' names_dat = data.frame(c("a", "b", "c"), c("d", "e", "f"), stringsAsFactors=FALSE)
#' changeByMap(c("a", "a", "b"), names_dat) == c("d", "d", "e")
#' x = changeByMap(c(NA, "a", "b"), names_dat)
#' is.na(x[1])
#' }
#' @return The new vector (mapped from the old one).
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
changeByMap = function(old_vector, mapping_dat, reverse = FALSE) {
	if(reverse) {
		mapping_dat = mapping_dat[, c(2, 1)]
	}
	sapply(old_vector, function(elem) {
				if(is.na(elem)) {
					return(NA)
				}
				if(! elem %in% mapping_dat[, 1]) {
					elem
				} else {
					mapping_dat[mapping_dat[, 1] == elem, 2]
				}
			})
}

