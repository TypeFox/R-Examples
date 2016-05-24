#     This file is part of FAiR, a program to conduct Factor Analysis in R
#     Copyright 2008 Benjamin King Goodrich
#
#     FAiR is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Affero General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     FAiR is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Affero General Public License for more details.
#
#     You should have received a copy of the GNU Affero General Public License
#     along with FAiR.  If not, see <http://www.gnu.org/licenses/>.

## This file contains GUI-related functions

## NOTE: This file is meant to be read with 90 columns and 8 space tabs

## get an answer from the user; everything below goes through this function
FAiR_get_answer <-
function(text, yesno = FALSE, radio_items = NULL, select = 1,
	seq_args = NULL, check_args = NULL, message = NULL) {

	options(guiToolkit = "RGtk2") # must use RGtk2 for GUI
	the.env <- new.env()
	group <- ggroup(horizontal = FALSE)
	glabel(text, container = group)
	if(!is.null(radio_items) | yesno) {
		if(yesno) {
			gseparator(horizontal = FALSE, container = group)
			gobject <- gradio(c("Yes", "No"), selected = select,
					horizontal = TRUE, container = group)
		}
		else gobject <- gradio(radio_items, selected = select, container = group)
	}
	else if(!is.null(seq_args)) {
		gobject <- gslider(from = seq_args$from, to = seq_args$to, 
			by = seq_args$by, value = seq_args$value, 
			horizontal = seq_args$horizontal, container = group)
	}
	else if(!is.null(check_args)) {
		gseparator(horizontal = FALSE, container = group)
		gobject <- gcheckboxgroup(items = check_args$items, 
					checked = check_args$checked,
					container = group)
	}
	if(!is.null(message)) {
		gseparator(horizontal = TRUE, container = group)
		glabel(message, container = group)
	}
	flush.console()
	cat("Note: The GUI may hide behind your other windows at any time.\n")
	flush.console()
	gbasicdialog(title = if(!is.null(check_args))
			"Please select (multiple selections possible)" else
			"Please select one",
			widget = group, #parent = group,
			handler = function(h,...) {
					assign("out", value = svalue(gobject),
						envir = h$action)
			}, action = the.env, toolkit = guiToolkit("RGtk2"))
	if(exists("out", envir = the.env)) {
		out <- get("out", envir = the.env)
		return(out)
	}
	else stop("Cancelled", call. = FALSE)
}

## ask user to edit stuff
FAiR_edit <-
function(instructions, mat) {
	options(guiToolkit = "RGtk2") # must use RGtk2 for GUI
	the.env <- new.env()
	group <- ggroup(horizontal = FALSE)
	glabel(instructions, container = group, anchor = c(-1,-1))
	flush.console()
	cat("One moment please\n")
	flush.console()
	obj <- gdf(as.data.frame(mat), container = group, expand = TRUE)
	gbasicdialog(title = "Editing needed",
			widget = group, #parent = group,
			handler = function(h,...) {
				assign("out", value = obj[,1:ncol(obj),drop = FALSE],
					envir = h$action)
			}, action = the.env, toolkit = guiToolkit("RGtk2"))
	if(exists("out", envir = the.env)) {
		out <- get("out", envir = the.env)
		if(any(sapply(out, is.character))) {
			for(i in 1:ncol(out)) out[,i] <- as.numeric(out[,i])
		}
		return(as.matrix(out))
	}
	else stop("Cancelled", call. = FALSE)
}


## ask user which model to estimate
FAiR_get_model <-
function() {
	items  <- c("Semi-Exploratory Factor Analysis", "Exporatory Factor Analysis",
			"Confirmatory Factor Analysis")
	answer <- FAiR_get_answer(text = "Which model would you like to estimate?",
				radio_items = items)
	if(answer == items[1]) answer <- "SEFA"
	else if(answer == items[2]) answer <- "EFA"
	else answer <- "CFA"
	return(answer)
}

## ask user which algorithm to use to estimate EFA model
FAiR_get_algorithm <- 
function() {
	text  <- "What algorithm would you like to use?"
	items <- c("call factanal() [fast]",paste("use the factanal() objective function",
		   "with RGENOUD [slower but more reliable]"), paste("CFA with zeros in",
		   "upper triangle of coefficient matrix [slower and less reliable]"))
	answer <- FAiR_get_answer(text, radio_items = items)
	return(which(items == answer))
}

## ask user a yes or no question
FAiR_yesno <- 
function(question, selected = 1, message = NULL) {
	answer <- FAiR_get_answer(text = question, yesno = TRUE, 
				select = selected, message = message)
	return(answer == "Yes")
}

## ask user to supply a number
FAiR_get_number <- 
function(text, from, to, by, value, message = NULL, horizontal = TRUE) {
	return(FAiR_get_answer(text, seq_args = list(from = from, to = to, by = by,
							value = value, message = message, 
							horizontal = horizontal)))
}

## establish arguments for the mapping rules defined in FAiR/R/Classes.R mapping_rule()
FAiR_get_mapping_rule_args <- 
function(x, level, choice_mr, choice_zeros){
	factors <- ncol(x)
	text <- paste("Are you sure you would like to change the mapping rule at level",
			level, "?")
	items <- c("No, I just want plain zeros", 
		"At least one zero per factor on a high communality variable", 
		if(factors > 2) "Encourage cohyperplanarity (revisionist Yates 1987)",
# 		if(factors > 2) "Encourage cohyperplanarity via MCD estimator",
		paste(factors, "outcomes of complexity", factors - 1, 
			"for each factor (weaker Thurstone)"), if(factors > 3)
		paste(0.5 * factors *  (factors - 1), "outcomes each of complexity",
					factors - 2),
		"Unit complexity basis (revisionist Butler 1969)", 
		"Set the maximum complexity for each outcome")
	names(items) <- c("none", "communality", if(factors > 2) "quasi_Yates", 
# 			if(factors > 2) "mve", 
			"weak_Thurstone", if(factors > 3) "viral",
			"Butler", "row_complexity")
	if(choice_mr) answer <- FAiR_get_answer(text = text, radio_items = items)
	else          answer <- items[1]
	arg_list <- formals(mapping_rule)
  mark <- names(arg_list) %in% names(which(items == answer))
	if(any(mark)) arg_list[[mark]] <- TRUE
	if(names(which(items == answer)) == "row_complexity") {
		text <- paste("What complexity would you like to use for *most*",
				"outcomes at level", level, "?")
		row_complexity <- FAiR_get_number(text, from = 1, by = 1, to = factors,
						  value = factors - 1, message = message)
		row_complexity <- matrix(as.integer(row_complexity), nrow = nrow(x))
		rownames(row_complexity) <- rownames(x)
		colnames(row_complexity) <- "complexity"
		text <- paste("Please edit the following matrix to specify the\n",
				"complexity of each outcome at level", level)
		row_complexity <- FAiR_edit(text, row_complexity)
		invalid <- any(row_complexity > factors) | any(row_complexity < 1)
		while(invalid) {
			gmessage(paste("Each complexity must be between 1 and",
					factors, "inclusive. Try again."))
			row_complexity <- FAiR_edit(text, row_complexity)
			invalid <- any(row_complexity > factors) | any(row_complexity < 1)
		}
		uniques <- unique(c(row_complexity)) 
		if(length(uniques) == 1) arg_list$row_complexity <- row_complexity[1]
		else arg_list$row_complexity <- c(row_complexity)
	}
	else if(names(which(items == answer)) == "mve") { # disabled
		text <- paste("Note this mapping rule is painfully slow and essentially",
				"requires that all outcomes are scored so as to be",
				"positively correlated with the primary factors", 
				"Click 'Cancel' at the next opportunity if this",
				"bothers you.")
		gmessage(text)
	}

	text <- paste(" Are you sure you would like to change the number\n",
			"of zero coefficients per factor at level", level, "from",
			factors, "to\n something else?" )
	message <- paste(" If so, after pressing OK, please specify the required\n",
			"number of exact zeros for each factor, *inclusive* of\n",
			"any zeros previously specified in specific cells.")
	if(choice_zeros) zeros <- FAiR_yesno(text, message = message)
	else             zeros <- FALSE

	if(zeros) {
		zeros <- rep(NA_integer_, factors)
		for(i in 1:length(zeros)) {
			text <- paste("Required number of zeros for factor", i,
					"at level", level)
			zeros[i] <- FAiR_get_number(text=text, from = factors - 1, by = 1,
						to = nrow(x) - 3, value = factors)
		}
	}
	else zeros <- rep(factors, factors)
	arg_list$zeros <- zeros
	return(as.list(arg_list))
}

## establish bounds on the off-diagonals of a correlation matrix
FAiR_bounds_cormat <- 
function(x, level, choice) {
	factors <- ncol(x)
	if(factors == 1) stop("cannot set bounds on a 1x1 correlation matrix")

	Domains <- array(cbind( array(-1 + .Machine$double.eps, c(factors, factors)),
				array( 1 - .Machine$double.eps, c(factors, factors))),
			dim = c(factors, factors, 2))
	mark <- upper.tri(Domains[,,1])
	for(i in 1:2) {
		diag(Domains[,,i]) <- 1
		Domains[,,i][mark] <- Inf
	}
	mark <- lower.tri(Domains[,,1])
	if(choice) {
		mat  <- cbind(lower = Domains[,,1][mark], upper = Domains[,,2][mark])
		rn <- as.character(NULL)
		letter <- if(level == 1) "F" else "G"
		for(i in 2:(factors)) for(j in 1:(i - 1)) {
			rn <- c(rn, paste(letter, i, " <-> ", letter, j, sep = ""))
		}
		rownames(mat) <- rn
		text <- paste("Please edit the bounds on the correlations between\n",
				"primary factors at level", level, ". Press OK when",
				"finished")
		invalid <- TRUE
		while(invalid) {
			mat <- FAiR_edit(text, mat)
			if(all(mat[,1] >= -1) & all(mat[,2] <= 1) &
			   all(mat[,1] <= mat[,2]) ) break
			message <- paste("All bounds on correlations must be on the",
					"[-1.0, 1.0] interval and consistent. Try again.")
			gmessage(message)
		}
		Domains[,,1][mark] <- mat[,1]
		Domains[,,2][mark] <- mat[,2]
	}
	return(Domains)
}

## establish bounds on a primary pattern coefficient matrix
FAiR_bounds_coef <- 
function(x, level, choice) {
	thing <- if(level == 1) 1 else if(ncol(x) == 1) 2 else 3
	items <- c(paste("Change lower bound (currently ", 
			switch(thing, -1.5, -1.0, -1.5), ")", sep = ""), 
			paste("Change upper bound (currently ",
			switch(thing,  1.5,  1.0,  1.5), ")", sep = "") )
	names(items) <- c("lower", "upper")
	text <- paste("Are you sure would you like to place non-trivial bounds\non",
			"some of the (free) coefficients at level", level, "?")
	check_args <- list(items = items, checked = c(TRUE, FALSE))
	if(choice) bounds <- FAiR_get_answer(text, check_args = check_args)
	else       bounds <- as.character(NULL)

	if(length(bounds) && "lower" %in% names(bounds)) {
		text <- paste("What number would you like to use as a lower bound",
				"for *most* of the coefficients at level", level, "?")
		num <- FAiR_get_number(text, from = switch(thing, -3, -1.0, -3), to = 0, 
					by = .1, value = switch(thing, -1.5, -1.0, -1.5))
		lowers <- matrix(num, nrow = nrow(x), ncol = ncol(x))
		rownames(lowers) <- rownames(x)
		colnames(lowers) <- colnames(x)
		lowers[!is.na(x)] <- x[!is.na(x)]
		invalid <- TRUE
		text <- paste("Please edit this matrix to specify lower bounds",
				"on each coefficient at level", level, "\n",
				"Recall that coefficients are parameterized for",
				"*standardized* variables.\nPress OK when finished.")
		while(invalid) {
			lowers <- FAiR_edit(text, lowers)
# 			lowers <- edit(lowers, edit.row.names = TRUE)
			if( (ncol(x) == 1) && any(lowers < -1) ) {
				gmessage("All bounds must be greater than -1.0")
			}
			else break
		}
	}
	else lowers <- array(switch(thing, -1.5, -1.0, -1.5), dim(x))

	if(length(bounds) && "upper" %in% names(bounds)) {
		text <- paste("What number would you like to use as an upper bound",
				"for *most* of the coefficients at level", level, "?")

		num <- FAiR_get_number(text, from = 0, to = switch(thing, 3, 1.0, 3), 
					by = .1, value = switch(thing, 1.5, 1.0, 1.5))
		text <- paste("Please edit this matrix to specify upper bounds",
				"on each coefficient at level", level, "\n",
				"Recall that coefficients are parameterized for",
				"*standardized* variables.\nPress OK when finished.")
		uppers <- matrix(num, nrow = nrow(x), ncol = ncol(x))
		uppers[!is.na(x)] <- x[!is.na(x)]
		rownames(uppers) <- rownames(x)
		colnames(uppers) <- colnames(x)
		invalid <- TRUE
		while(invalid) {
			uppers <- FAiR_edit(text, uppers)
# 			uppers <- edit(uppers, edit.row.names = TRUE)
			if( (ncol(x) == 1) && any(lowers > 1) ) {
				gmessage("All bounds must be less than 1.0")
			}
			else if(any(lowers > uppers)) {
				gmessage(paste("All upper bounds must be greater than or",
						"equal to the lower bounds"))
			}
			else break
		}
	}
	else uppers <- array(switch(thing, 1.5, 1.0, 1.5), dim(x))

	Domains <- array(cbind(lowers, uppers), dim = c(dim(x), 2))
# 	rn <- as.character(NULL)
# 	for(i in 1:ncol(x)) rn <- c(rn, paste(rownames(x), "_", i, sep = ""))
# 	rownames(Domains) <- row_names
# 	Domains <- Domains[coef@free,]
	return(Domains)
}

## fix specific cells of a primary pattern coefficient matrix to specific numbers
FAiR_peg_coefficients <- 
function(coefs, level) {
	text <- paste("Please edit this *standardized* pattern matrix",
			"at level", level, "\nto designate particular coefficients",
			"as fixed parameters.\nUse finite numbers to fix a",
			"coefficient to that number;\nuse Inf if the coefficient",
			"will be a function of other coefficients.\nLeave",
			"any *unrestricted* coefficients as NA.\nPress OK",
			"when finished")
	rows <- nrow(coefs)
	cols <- ncol(coefs)
	coefs <- as.data.frame(coefs)
	for(i in 1:cols) coefs[,i] <- as.character(coefs[,i])
	coefs <- FAiR_edit(text, coefs)
	invalid <- cols == 1
	while(invalid) {
		if(all(coefs >= -1.0, na.rm = TRUE) & 
		   all(coefs <=  1.0, na.rm = TRUE)) break
		gmessage("All coefficients must be on the [-1.0,1.0] interval, try again")
		coefs <- FAiR_edit(text, coefs)
	}
	return(coefs[1:rows,1:cols,drop = FALSE])
}

## impose order restrictions
FAiR_order_FC <- 
function(orders, level, MARGIN) {
	rows <- nrow(orders)
	cols <- ncol(orders)
	the_max <- if(MARGIN == 1) cols else rows
	text <- paste(" This matrix is consistent with ALL configurations of",
			"factor\ncontributions at level", level, ". Edit it",
			"to restrict the rank orderings\nof factor contributions within",
			if(MARGIN == 1) "*rows*" else "*columns*", "using integers in",
			"[1,", the_max, "].\n\n",
			"Larger factor contributions correspond to smaller ranks;\n",
			"i.e., the biggest factor contribution should be ranked 1.\n",
			"Ranks need not be unique; i.e. if either of two factor\n",
			"contributions could be largest, rank them both as 2.\n",
			"When finished, press OK.")
	invalid <- TRUE
	if(level == 2) {
		rownames(orders) <- paste("F", 1:nrow(orders), sep = "")
		colnames(orders) <- paste("G", 1:ncol(orders), sep = "")
	}
	else    colnames(orders) <- paste("F", 1:ncol(orders), sep = "")

	while(invalid) {
		orders <- FAiR_edit(text, orders)
		if(all(c(orders) %in% c(1:the_max))) break
		message <- paste("all entries must be integers between 1 and", the_max)
		gmessage(message)
	}
	return(orders[1:rows,1:cols,drop = FALSE])
}

## designate best variables (1st-order factors) at level 1 (2)
FAiR_indicators_FC <-
function(indicators, level) {
	rows <- nrow(indicators)
	cols <- ncol(indicators)
	text <- paste("Change the appropriate cell of the following matrix to TRUE to",
			"indicate\nwhich", if(level == 1) "manifest variable" else
			"factor at level 1", "is the best indicator of a factor at level",
			level, ".\nIf you do not wish to take such a position for some",
			"factor,\nleave all cells in that column as FALSE. When",
			"finished, press OK.")
	invalid <- TRUE
	while(invalid) {
		indicators <- FAiR_edit(text, indicators)
		if(all(colSums(indicators) <= 1)) break
		message <- paste("there should be at most one TRUE per column, the rest",
					"must be FALSE; try again.")
		gmessage(message)
	}
	out <- rep(NA_integer_, ncol(indicators))
	for(i in 1:ncol(indicators)) {
		if(any(indicators[,i])) out[i] <- which(indicators[,i])
		else out[i] <- NA_integer_
	}
	return(out)
}

## set up blocking matrix
FAiR_blockanator <-
function(fixed, level) {
	rows <- nrow(fixed)
	cols <- ncol(fixed)
	mat  <- matrix(FALSE, nrow = nrow(fixed), ncol = ncol(fixed))
	mat[!is.na(fixed) & as.logical(fixed)] <- TRUE
	dimnames(mat) <- dimnames(fixed)
	text  <- paste("Indicate which coefficients at level", level, "must *not* be",
			"zero by changing\nthe appropriate cells to TRUE.",
			"Press OK when finished")
	mat <- FAiR_edit(text, mat)
	invalid <- !is.logical(mat)
	while(invalid) {
		if(is.logical(mat)) break
		gmessage("all cells must be TRUE or FALSE")
		mat <- !(!mat)
		dimnames(mat) <- dimnames(fixed)
		mat <- FAiR_edit(text, mat)
	}
	return(mat[1:rows,1:cols,drop = FALSE])
}

## establish equality constraints among cells in a primary pattern coefficient matrix
FAiR_equality_restrictions <- 
function(fixed, level, auto = FALSE) {
	out <- list()
	if(auto) return(out)
	count <- 1
	while(TRUE) {
	text <- paste(" Please change cells that must be equal to each other at level",
			level, "to TRUE and then press OK.\n",
			"You will have multiple opportunities to designate sets of",
			"equal coefficients.", "This is set", count, ".\n",
			"When finished, leave all cells as FALSE and press OK")
		mat <- matrix(FALSE, nrow = nrow(fixed), ncol = ncol(fixed),
				dimnames = dimnames(fixed))
		mat <- FAiR_edit(text, mat)
		trues <- sum(mat)
		if(trues == 1) {
			message <- paste("The number of cells that are TRUE must be",
				"greater than one (to indicate equality restrictions)",
				"or zero (to indicate you are finished imposing equality",
				"restrictions). Try again")
			gmessage(message)
			next
		}
		else if(trues == 0) break
		out[[count]] <- new("equality_restriction", free = which(mat)[1],
					fixed = which(mat)[-1], dims = dim(fixed),
					rownames = rownames(fixed), 
					level = as.integer(level))
		count <- count + 1
		cat("Thank you sir, may I have another set of equality restrictions?\n")
	}
	return(out)
}

FAiR_inequalities <-
function(level_1, level_2, factors, SEFA) {
	## Level 1
	constraints <- FAiR_constraints_1st(onecol = factors[1] == 1)
	items <- sapply(constraints, FUN = function(x) x[[2]])
	names(items) <- names(constraints)
	if(!SEFA) constraints$block_1st <- NULL
	checked <- rep(FALSE, length(items))
	text <- paste( "Please indicate which, if any, restrictions you would like",
			"to impose at level 1")

	if(level_1) answer <- FAiR_get_answer(text, check_args = list(items = items, 
				  		checked = checked))
	else        answer <- as.character(NULL)

	foo <- function(x) x[[1]]
	if(length(answer)) criteria <- sapply(constraints[names(answer)], FUN = foo)
	else               criteria <- list()

	if(!level_2) return(criteria)

	## Level 2
	constraints <- FAiR_constraints_2nd(onecol = factors[2] == 1)
	items <- sapply(constraints, FUN = function(x) x[[2]])
	names(items) <- names(constraints)
	if(!SEFA) constraints$block_2nd <- NULL
	checked <- rep(FALSE, length(items))
	text <- paste( "Please indicate which, if any, restrictions you would like",
			"to impose at level 2")

	answer <- FAiR_get_answer(text, check_args = list(items = items, 
							checked = checked))

	if(length(answer)) criteria2 <- sapply(constraints[names(answer)], FUN = foo)
	else               criteria2 <- list()

	criteria <- c(criteria2, criteria)
	return(criteria)
}

## get criteria to use in SEFA and CFA models
FAiR_criterionator_extraction <-
function(criteria, methodArgs, discrepancy, factors, manifest) {
	if(!is.list(criteria)) stop("'criteria' must be a list")
	else if(length(criteria) > 0) {
		CRITERIA <- c(FAiR_constraints_2nd(), FAiR_constraints_1st())
		for(i in 1:length(criteria)) {
			if(is.character(criterionname <- criteria[[i]])) {
				criteria[[i]] <- CRITERIA[[criteria[[i]]]][[1]]
				names(criteria)[i] <- criterionname
			}
			else if(!is.function(criteria[[i]])) {
				stop("criteria must be a list of character",
					" strings or functions, it is usually best",
					" to leave criteria unspecified")
			}
		}
	}
	if(!(discrepancy %in% c("MLE", "YWLS")) && is(manifest, "manifest.data") &&
		is(manifest@acov, "diagonalMatrix")) discrepancy <- "DWLS"

	criteria[[length(criteria) + 1]] <- FAiR_discrepancy(discrepancy)
	if(discrepancy == "SHK") {
		criteria2 <- try(FAiR_constantanator(manifest, "SHK", criteria))
		if(!is.list(criteria2)) {
			discrepancy <- "ELLIPTICAL"
			l <- length(criteria)
			criteria[[l]] <- FAiR_discrepancy(discrepancy)
			criteria <- FAiR_constantanator(manifest, discrepancy, criteria)
		}
		else criteria <- criteria2
	}
	else criteria <- FAiR_constantanator(manifest, discrepancy, criteria)
	names(criteria)[length(criteria)] <- discrepancy
	criteria <- FAiR_fill_methodArgs_extraction(criteria, factors, manifest,
							methodArgs)
	return(criteria)
}

## this does the first menu that pops up
# FAiR_menunator <-
# function(factors, auto) {
# 	if(factors[2] >= 2) {
# 		items2 <- c(paste("Change bounds on correlations among primary factors",
# 				  "at level 2 (currently [-1,1] for each)"), 
# 				paste("Change bounds on coefficients at level 2",
# 				"(currently [-1.5,1.5] for each)"),
# 				"Impose functional inequality restrictions at level 2")
# 		names(items2) <- c("2nd_bounds_cormat", "2nd_bounds_coef",
# 					"2nd_inequalities")
# 		items_list <- list()
# 	}
# 	else if(factors[2] == 1) {
# 		items2 <- c(paste("Change bounds on coefficients at level 2",
# 				"(currently [-1.0,1.0] for each)"),
# 				"Impose functional inequality restrictions at level 2")
# 		names(items2) <- c("2nd_bounds_coef", "2nd_inequalities")
# 		items_list <- list()
# 	}
# 	else { # one-level model
# 		items2 <- paste("Change bounds on correlations among primary factors at",
# 				"level 1 (currently [-1,1] for each)")
# 		items_list <- list(TRUE)
# 		names(items_list[[1]]) <- names(items2) <- "1st_bounds_cormat"
# 	}
# 
# 	bound <- ifelse(factors[1] > 1, 1.5, 1.0)
# 	items1 <- c(paste("Change bounds on coefficients at level 1",
# 				"(currently [", -bound, ",", bound, "] for each)"),
# 				"Impose functional inequality restrictions at level 1")
# 	names(items1) <- c("1st_bounds_coef", "1st_inequalities")
# 
# 	items <- c(items2, items1)
# 	if(factors[2] == 0) {
# 		items_list[[1]] <- c(items_list[[1]], rep(TRUE, length(items1)))
# 		names(items_list[[1]]) <- c(names(items_list[[1]][1]), names(items1))
# 	}
# 	else {
# 		items_list <- list(rep(TRUE, length(items1)), rep(TRUE, length(items2)))
# 		names(items_list[[1]]) <- names(items1)
# 		names(items_list[[2]]) <- names(items2)
# 	}
# 
# 	if(auto) return(items_list)
# 	text <- paste("Would you like to impose inequality restrictions on the model?\n",
# 		"Anything you specified at the command prompt will be implemented",
# 		"automatically,\nbut you can mark the appropriate boxes to verify it")
# 	checked <- rep(FALSE, length(items))
# 	names(checked) <- names(items)
# 	check_args <- list(items = items, checked = checked)
# 	answer <- FAiR_get_answer(text, check_args = check_args)
# 	if(any(names(answer) %in% names(items1))) {
# 		temp <- intersect(names(answer), names(items1))
# 		items_list[[1]][temp] <- FALSE
# 	}
# 	if(factors[2] > 0 && any(names(answer) %in% names(items2))) {
# 		temp <- intersect(names(answer), names(items2))
# 		items_list[[2]][temp] <- FALSE
# 	}
# 	else if(factors[2] == 0 && "1st_bounds_cormat" %in% names(answer)) {
# 		items_list[[1]]["1st_bounds_cormat"] <- FALSE
# 	}
# 	return(items_list)
# }

## this does the first menu that pops up
FAiR_menunator <-
function(factors, SEFA) {
	if(factors[2] >= 2) {
		items2 <- c(paste("Change bounds on correlations among primary factors",
				  "at level 2 (currently [-1,1] for each)"), 
			paste("Fix specific coefficients at level 2 to particular values",
				"(or designate them unfree for other reasons)"),
			paste("Change bounds on coefficients at level 2",
				"(currently [-1.5,1.5] for each)"),
			"Constrain coefficients at level 2 to be equal to each other",
			"Use a nondefault mapping rule at level 2",
			paste("Change the number of zero coefficients required",
			"for each factor at level 2 (currently", factors[2], "for each)"),
			"Impose functional inequality restrictions at level 2")
		names(items2) <- c("bounds_cormat", "peg_coef",
				"bounds_coef", "equalities",
				 "mapping_rule", "zeros", "inequalities")
		if(!SEFA) items2 <- items2[c("bounds_cormat", "equalities",
						"bounds_coef", "inequalities")]
		items_list <- list()
	}
	else if(factors[2] == 1) {
		items2 <- c(paste("Fix specific coefficients at level 2 to particular",
				"values (or designate them unfree for other reasons)"),
			paste("Change bounds on coefficients at level 2",
				"(currently [-1.0,1.0] for each)"),
			"Constrain coefficients at level 2 to be equal to each other",
			"Impose functional inequality restrictions at level 2")
		names(items2) <- c("peg_coef", "bounds_coef",
				"equalities", "inequalities")
		items_list <- list()
	}
	else { # one-level model
		items2 <- paste("Change bounds on correlations among primary factors at",
				"level 1 (currently [-1,1] for each)")
		items_list <- list(FALSE)
		names(items_list[[1]]) <- names(items2) <- "bounds_cormat"
	}

	bound <- ifelse(factors[1] > 1, 1.5, 1.0)
	SEFA2 <- SEFA && factors[1] > 1
	items1 <- c(paste("Fix specific coefficients at level 1 to particular values",
				"(or designate them unfree for other reasons)"),
			paste("Change bounds on coefficients at level 1",
				"(currently [", -bound, ",", bound, "] for each)"),
			"Constrain coefficients at level 1 to be equal to each other",
			"Use a nondefault mapping rule at level 1",
			paste("Change the number of zero coefficients required",
			"for each factor at level 1 (currently", factors[1], "for each)"),
			"Impose functional inequality restrictions at level 1")
	names(items1) <- c("peg_coef", "bounds_coef", "equalities",
			"mapping_rule", "zeros", "inequalities")
	if(!SEFA) items1 <- items1[c("bounds_coef", "equalities", "inequalities")]

	items <- c(items2, items1)
	if(factors[2] == 0) { 
		items_list[[1]] <- c(items_list[[1]], rep(FALSE, length(items1)))
		names(items_list[[1]]) <- c(names(items_list[[1]][1]), names(items1))
	}
	else {
		items_list <- list(rep(FALSE, length(items1)), rep(FALSE, length(items2)))
		names(items_list[[1]]) <- names(items1)
		names(items_list[[2]]) <- names(items2)
	}

	text <- paste("Would you like to change the restrictions imposed on the model?\n",
		"Anything you specified at the command prompt will be implemented",
		"automatically,\nbut you can mark the appropriate boxes to verify")
	checked <- rep(FALSE, length(items))
	check_args <- list(items = items, checked = checked)
	answer <- FAiR_get_answer(text, check_args = check_args)
	if(any(answer %in% items1)) {
		temp <- names(items1)[items1 %in% answer]
		items_list[[1]][temp] <- TRUE
	}
	if(factors[2] > 0 &&  any(answer %in% items2)) {
		temp <- names(items2)[items2 %in% answer]
		items_list[[2]][temp] <- TRUE
	}
	else if(factors[2] == 0 && "bounds_cormat" %in% names(answer)) {
		items_list[[1]]["bounds_cormat"] <- TRUE
	}
	return(items_list)
}

FAiR_GPAnator <-
function(Lambda, method = NULL, methodArgs) {

	items <- c("geomin", "quartimin", "target", "partially-specified target", 
			"oblimax", "simplimax", "factor simplicity (Bentler)",
			"Crawford-Ferguson", "infomax", "entropy (McCammon)", "oblimin")
	names(items) <- c("geomin", "quartimin", "target", "pst", "oblimax",
			"simplimax", "bentler", "cf", "infomax", "mccammon", "oblimin")
	text <- "Which oblique criterion would you like to use?"
	if(is.null(method) | method == "GPA") {
		answer <- FAiR_get_answer(text, radio_items = items)
		method <- names(items)[answer == items]
	}
	else method <- match.arg(method, names(items))

	if(is.null(matrix <- methodArgs$matrix)) {
		text <- "Which matrix should be used for the optimization?"
		items <- c("Primary Pattern", "Reference Structure", 
				"Factor Contribution")
		names(items) <- c("PP", "RS", "FC")
		answer <- FAiR_get_answer(text, radio_items = items)
		matrix <- names(items)[items == answer]
	}
	else {
		matrix <- match.arg(matrix, c("PP", "RS", "FC"))
		methodArgs$matrix <- NULL
	}

	if(method == "geomin") {
		if(is.null(methodArgs$delta)) {
			text <- paste("What small constant would you like to use for the",
				"geomin criterion?", if(ncol(Lambda) > 3)
				"\nThe default of 0.01 might be too small")
			methodArgs$delta <- FAiR_get_number(text, from = 0, to = 0.1,
								by = .005, value = .01)
		}
	}
	else if(method == "target") {
		if(is.null(methodArgs$Target)) {
			methodArgs$Target <- FAiR_make_target(Lambda)
		}
		else valid <- FAiR_check_target(Lambda, methodArgs$Target)
		methodArgs <- list(Target = methodArgs$Target)
	}
	else if(method == "pst") {
		if(is.null(methodArgs$Target)) {
			methodArgs$Target <- FAiR_make_pst(Lambda)
		}
		else valid <- FAiR_check_pst(Lambda, methodArgs$Target)

		if(!is.null(methodArgs$W)) {
			W <- methodArgs$W
			methodArgs$Target[W == 0] <- NA_real_
			warning("an unnecessary 'W' matrix was passed and some cells of",
				" the target matrix were converted to NA")
		}
		methodArgs <- list(Target = methodArgs$Target)
	}
	else if(method == "simplimax") {
		if(is.null(methodArgs$k)) {
			text <- paste("How many near zeros would you like the simplimax",
					"criterion to seek?")
			rows <- nrow(Lambda)
			methodArgs$k <- FAiR_get_number(text, from = 1, value = rows,
							to = length(Lambda), by = 1) 
		}
		else if(methodArgs$k < 0) {
			stop("'k' must be positive")
		}
		else if(methodArgs$k > length(Lambda)) {
			stop("'k' must be smaller than the total number of loadings")
		}
		methodArgs <- list(k = methodArgs$k)
	}
	else if(method == "cf") {
		if(is.null(methodArgs$kappa)) {
			text <- paste("What parameter would you like to use for the",
					"Crawford-Ferguson criterion?")
			methodArgs$kappa <- FAiR_get_number(text, from = 0, to = 1, 
								by = .01, value = 0)
		}
		else if(methodArgs$kappa < 0) {
			stop("'kappa' must be non-negative")
		}
		else if(methodArgs$kappa > 1) {
			stop("'kappa' must be less than or equal to 1.0")
		}
		methodArgs <- list(kappa = methodArgs$kappa)
	}
	else if(method == "oblimin") {
		if(is.null(methodArgs$gam)) {
			text <- paste("What parameter you would like to use for the",
					"oblimin criterion?")
			methodArgs$gam <- FAiR_get_number(text, from = 0, to = 1, 
								by = .01, value = 0)
		}
		else if(methodArgs$gam < 0) {
			stop("'gam' must be non-negative")
		}
		else if(methodArgs$gam > 1) {
			stop("'gam' must be less than or equal to 1.0")
		}
		methodArgs <- list(gam = methodArgs$gam)
	}

	criterion <- FAiR_analytic("GPA")[[1]]
	formals(criterion)$matrix <- matrix
	formals(criterion)$method <- method
	formals(criterion)$methodArgs <- methodArgs
	return (criterion)
}

## make a target matrix
FAiR_make_target <-
function(Lambda) {
	target <- matrix(0, nrow = nrow(Lambda), ncol = ncol(Lambda))
	rownames(target) <- rownames(Lambda)
	colnames(target) <- colnames(Lambda)
	text <- paste("Please edit this matrix to specify the *target* values.\n",
			"Press OK when finished")
	target <- FAiR_edit(text, target)
	while(any(is.na(target))) {
		message <- paste("Target must be *fully* specified, without any NAs",
				"Please try again")
		gmessage(message)
		target <- FAiR_edit(text, target)
	}
	return(target[1:nrow(Lambda),1:ncol(Lambda),drop = FALSE])
}

## make a partially-specified target matrix
FAiR_make_pst <-
function(Lambda) {
	target <- matrix(NA_real_, nrow = nrow(Lambda), ncol = ncol(Lambda))
	rownames(target) <- rownames(Lambda)
	colnames(target) <- colnames(Lambda)
	text <- paste("Please edit this matrix to specify the *target* values.\n",
			"Leave any untargeted cells as NA.\n",
			"Press OK when finished")
	target <- FAiR_edit(text, target)
	return(target[1:nrow(Lambda),1:ncol(Lambda),drop = FALSE])
}

## fills in methodArgs for make_restrictions, possibly with GUI
FAiR_fill_methodArgs_extraction <-
function(criteria, factors, manifest, methodArgs) { 
	fixed  <- matrix(NA_real_, nrow = nrow(cormat(manifest)), ncol = factors[1])
	fixed2 <- matrix(NA_real_, nrow = factors[1], ncol = factors[2])

	if(length(mark <- which(names(criteria) == "ranks_rows_1st"))) {
		if(is.null(row_ranks <- methodArgs$row_ranks_1st) &&
		   is.null(row_ranks <- methodArgs$row_ranks)) {
			orders <- matrix(factors[1], nrow(fixed), factors[1])
			rownames(orders) <- rownames(fixed)
			colnames(orders) <- colnames(fixed)
			row_ranks <- FAiR_order_FC(orders, level = 1, MARGIN = 1)
		}
# 		else methodArgs$row_ranks <- NULL

		if(!identical(dim(row_ranks), dim(fixed))) {
			stop("'row_ranks' must have the same dimension",
				"as the loadings matrix at level 1")
		}

		vals <- unique(c(row_ranks))
		if(any(! (vals %in% c(1:factors[1])))) {
			stop("'row_ranks' must only contain positive integers up to ",
				factors[1])
		}

		if(all(row_ranks == factors[1])) {
			criteria <- criteria[-mark]
			warning("row ranks constraint omitted because no ",
				"constraints were specified")
		}
		else {
			formals(criteria[[mark]])$user_ranks <- row_ranks
			formals(criteria[[mark]])$MARGIN <- 1
		}
	}

	if(length(mark <- which(names(criteria) == "ranks_rows_2nd"))) {
		if(is.null(row_ranks <- methodArgs$row_ranks_2nd) &&
		   is.null(row_ranks <- methodArgs$row_ranks)) {
			orders <- matrix(factors[2], factors[1], factors[2])
			rownames(orders) <- rownames(fixed2)
			colnames(orders) <- colnames(fixed2)
			row_ranks <- FAiR_order_FC(orders, level = 2, MARGIN = 1)
		}
# 		else methodArgs$row_ranks <- NULL

		if(!identical(dim(row_ranks), dim(fixed2))) {
			stop("'row_ranks' must have the same dimension",
				"as the loadings matrix at level 2")
		}

		vals <- unique(c(row_ranks))
		if(any(! (vals %in% c(1:factors[2])))) {
			stop("'row_ranks' must only contain positive integers up to ",
				factors[1])
		}

		if(all(row_ranks == factors[1])) {
			criteria <- criteria[-mark]
			warning("row ranks constraint omitted because no ",
				"constraints were specified")
		}
		else {
			formals(criteria[[mark]])$user_ranks <- row_ranks
			formals(criteria[[mark]])$MARGIN <- 1
		}
	}

	if(length(mark <- which(names(criteria) == "ranks_cols_1st"))) {
		if(is.null(col_ranks <- methodArgs$col_ranks_1st) &&
		   is.null(col_ranks <- methodArgs$col_ranks)) {
			orders <- matrix(nrow(fixed), nrow(fixed), factors[1])
			rownames(orders) <- rownames(fixed)
			colnames(orders) <- colnames(fixed)
			col_ranks <- FAiR_order_FC(orders, level = 1, MARGIN = 2)
		}
# 		else methodArgs$col_ranks <- NULL

		if(!identical(dim(col_ranks), dim(fixed))) {
			stop("'col_ranks' must have the same dimension",
				"as the loadings matrix at level 1")
		}
		vals <- unique(c(col_ranks))
		if(any(! (vals %in% c(1:nrow(fixed))))) {
			stop("'col_ranks' must only contain positive integers up to ",
				nrow(fixed))
		}

		if(all(col_ranks == nrow(fixed))) {
			criteria <- criteria[-mark]
			warning("column ranks constraint omitted because no ",
				"constraints were specified")
		}
		else {
			formals(criteria[[mark]])$user_ranks <- col_ranks
			formals(criteria[[mark]])$MARGIN <- 2
		}
	}

	if(length(mark <- which(names(criteria) == "ranks_cols_2nd"))) {
		if(is.null(col_ranks <- methodArgs$col_ranks_2nd) &&
		   is.null(col_ranks <- methodArgs$col_ranks)) {
			orders <- matrix(nrow(fixed2), nrow(fixed2), factors[2])
			rownames(orders) <- rownames(fixed2)
			colnames(orders) <- colnames(fixed2)
			col_ranks <- FAiR_order_FC(orders, level = 2, MARGIN = 2)
		}
# 		else methodArgs$col_ranks <- NULL

		if(!identical(dim(col_ranks), dim(fixed2))) {
			stop("'col_ranks' must have the same dimension",
				"as the loadings matrix at level 2")
		}
		vals <- unique(c(col_ranks))
		if(any(! (vals %in% c(1:nrow(fixed2))))) {
			stop("'col_ranks' must only contain positive integers up to ",
				nrow(fixed2))
		}

		if(all(col_ranks == nrow(fixed2))) {
			criteria <- criteria[-mark]
			warning("column ranks constraint omitted because no ",
				"constraints were specified")
		}
		else {
			formals(criteria[[mark]])$user_ranks <- col_ranks
			formals(criteria[[mark]])$MARGIN <- 2
		}
	}

	if(length(mark <- which(names(criteria) == "indicators_1st"))) {
		if(is.null(indicators <- methodArgs$indicators_1st) &&
		   is.null(indicators <- methodArgs$indicators)) {
			indicators <- matrix(FALSE, nrow(fixed), ncol(fixed))
			dimnames(indicators) <- dimnames(fixed)
			indicators <- FAiR_indicators_FC(indicators, level = 1)
		}

		if(any(!(indicators %in% c(1:nrow(fixed), NA)))) {
			stop("all entries in 'indicators' must be NA or positive",
				" integers less than ", nrow(fixed))
		}
		formals(criteria[[mark]])$indicators <- indicators
	}

	if(length(mark <- which(names(criteria) == "indicators_2nd"))) {
		if(is.null(indicators <- methodArgs$indicators_2nd) &&
		   is.null(indicators <- methodArgs$indicators)) {
			indicators <- matrix(FALSE, nrow(fixed2), ncol(fixed2))
			dimnames(indicators) <- dimnames(fixed2)
			indicators <- FAiR_indicators_FC(indicators, level = 2)
		}

		if(any(!(indicators %in% c(1:nrow(fixed2), NA)))) {
			stop("all entries in 'indicators' must be NA or positive",
				" integers less than ", nrow(fixed2))
		}
		formals(criteria[[mark]])$indicators_2nd <- indicators
	}

	if(length(mark <- which(names(criteria) == "no_neg_suppressors_1st"))) {
		text <- paste("Select the minimum acceptable factor contribution",
			"coefficient\nsuch that lower values indicate a",
			"suppressor variable at level 1.")
		if(is.null(FC_threshold <- methodArgs$FC_threshold_1st) &&
		   is.null(FC_threshold <- methodArgs$FC_threshold)) {
			FC_threshold <- FAiR_get_number(text, from = -1, to = 0,
							by = .001, value = -.01)
		}
# 		else methodArgs$FC_threshold <- NULL

		if(FC_threshold <= -1) stop("'FC_threshold' must be > -1")
		if(FC_threshold > 0)   stop("'FC_threshold' must be < 0")
		formals(criteria[[mark]])$FC_threshold <- FC_threshold
	}

	if(length(mark <- which(names(criteria) == "no_neg_suppressors_2nd"))) {
		text <- paste("Select the minimum acceptable factor contribution",
			"coefficient\nsuch that lower values indicate a",
			"suppressor variable at level 2.")
		if(is.null(FC_threshold <- methodArgs$FC_threshold_2nd) &&
		   is.null(FC_threshold <- methodArgs$FC_threshold)) {
			FC_threshold <- FAiR_get_number(text, from = -1, to = 0,
							by = .001, value = -.01)
		}
# 		else methodArgs$FC_threshold <- NULL

		if(FC_threshold <= -1) stop("'FC_threshold' must be > -1")
		if(FC_threshold > 0)   stop("'FC_threshold' must be < 0")
		formals(criteria[[mark]])$FC_threshold_2nd <- FC_threshold
	}

	if(length(mark <- which(names(criteria) == "dist_cols_1st"))) {
		text <- paste("Select the minimum acceptable binary distance between",
				"columns of the logically negated loadings matrix",
				"at level 1.")
		if(is.null(cutpoint <- methodArgs$cutpoint_1st) &&
		   is.null(cutpoint <- methodArgs$cutpoint)) {
			cutpoint <- FAiR_get_number(text, from = 0, to = 1,
							by = .01, value = 0.5)
		}

		if(cutpoint <= 0) stop("'cutpoint' must be > 0")
		if(cutpoint >  1) stop("'cutpoint' must be <= 1")
		formals(criteria[[mark]])$cutpoint <- cutpoint
	}

	if(length(mark <- which(names(criteria) == "dist_cols_2nd"))) {
		text <- paste("Select the minimum acceptable binary distance between",
				"columns of the logically negated loadings matrix",
				"at level 2.")
		if(is.null(cutpoint <- methodArgs$cutpoint_2nd) &&
		   is.null(cutpoint <- methodArgs$cutpoint)) {
			cutpoint <- FAiR_get_number(text, from = 0, to = 1,
							by = .01, value = 0.5)
		}

		if(cutpoint <= 0) stop("'cutpoint' must be > 0")
		if(cutpoint >  1) stop("'cutpoint' must be <= 1")
		formals(criteria[[mark]])$cutpoint_2nd <- cutpoint
	}

	if(length(mark <- which(names(criteria) == "volume_1st"))) {
		S <- cormat(manifest)
		formals(criteria[[mark]])$det_S <- det(S)
	}

	if(length(mark <- which(names(criteria) == "block_2nd"))) {
		if(is.null(blockers <- methodArgs$blockers_2nd) &&
		   is.null(blockers <- methodArgs$blockers)) {
			blockers <- FAiR_blockanator(fixed2, level = 2)
		}

		if(!identical(dim(blockers), dim(fixed2))) {
			stop("the dimensions of blockers must be the same as the ",
				"loading matrix at level 1")
		}
		if(!is.logical(blockers)) {
			stop("all cells of 'blockers' must be 'NA', 'TRUE', or 'FALSE'")
		}
		if(!any(blockers, na.rm = TRUE)) {
			criteria <- criteria[-mark]
			warning("blockers constraint omitted at level 2 because no ",
				"coefficients were blocked")
		}
		else formals(criteria[[mark]])$blockers_2nd <- blockers
	}

	if(length(mark <- which(names(criteria) == "block_1st"))) {
		if(is.null(blockers <- methodArgs$blockers_1st) &&
		   is.null(blockers <- methodArgs$blockers)) {
			blockers <- FAiR_blockanator(fixed, level = 1)
		}

		if(!identical(dim(blockers), dim(fixed))) {
			stop("the dimensions of blockers must be the same as the ",
				"loading matrix at level 1")
		}
		if(!is.logical(blockers)) {
			stop("all cells of 'blockers' must be 'NA', 'TRUE', or 'FALSE'")
		}
		if(!any(blockers, na.rm = TRUE)) {
			criteria <- criteria[-mark]
			warning("blockers constraint omitted at level 1 because no ",
				"coefficients were blocked")
		}
		else formals(criteria[[mark]])$blockers <- blockers
	}
	return(criteria)
}

## fills in methodArgs for Rotate, possibly with GUI
FAiR_fill_methodArgs_Rotate <-
function(FAobject, criteria, methodArgs) {
	Lambda  <- loadings(FAobject)
	factors <- ncol(Lambda)

	if(length(mark <- which(names(criteria) == "no_factor_collapse"))) {
		text <- paste("Select the minimum acceptable effective variance",
				"among primary factors.",
				"\nSetting this threshold too close to zero risks",
				"factor collapse.")
		if(is.null(nfc_threshold <- methodArgs$nfc_threshol)) {
			nfc_threshold <- FAiR_get_number(text, from = 0, to = 1, 
							by = .01, value = .25)
		}
		else methodArgs$nfc_threshold <- NULL
		if(nfc_threshold <  0) stop("'nfc_threshold' must be >= 0")
		if(nfc_threshold >= 1) stop("'nfc_threshold' must be < 1")
		formals(criteria[[mark]])$nfc_threshold <- nfc_threshold
	}

	if(length(mark <- which(names(criteria) == "limit_correlations"))) {
		text  <- paste("Select the minimum acceptable correlation",
				"among primary factors")
		if(is.null(lower <- methodArgs$lower)) { 
			lower <- FAiR_get_number(text, from = -1, to = 1,
							by = .01, value = -1)
		}
		else methodArgs$lower <- NULL

		if(lower < -1) stop("'lower' must be >= -1")
		if(lower >= 1) stop("'lower' must be < 1")
		formals(criteria[[mark]])$lower <- lower

		text  <- paste("Select the maximum acceptable correlation",
				"among primary factors")
		if(is.null(upper <- methodArgs$upper)) { 
			upper <- FAiR_get_number(text, from = lower, to = 1,
							by = .01, value = 1)
		}
		else methodArgs$upper <- NULL

		if(lower > upper) stop("'lower' must be <= 'upper'")
		if(upper > 1)     stop("'upper' must be <= 1")
		formals(criteria[[mark]])$upper <- upper

		if(lower == 0 && upper == 0) {
			warning("Rotate() does not do orthogonal rotation.",
				"Ignoring bounds on inter-factor correlations.")
			criteria <- criteria[-mark]
		}
	}

	if(length(mark <- which(names(criteria) == "positive_manifold"))) {
		text <- paste("Select the minimum acceptable reference structure",
				"correlation.\nA slightly negative number is",
				"recommended rather than zero")
		if(is.null(pm_threshold <- methodArgs$pm_threshold)) {
			pm_threshold <- FAiR_get_number(text, from = -1, to = 0, 
							by = .01, value = -.1)
		}
		else methodArgs$pm_threshold <- NULL

		if(pm_threshold <= -1) stop("'pm_threshold' must be > -1")
		if(pm_threshold > 0)   stop("'pm_threshold' must be < 0")
		formals(criteria[[mark]])$pm_threshold <- pm_threshold
	}

	if(length(mark <- which(names(criteria) == "ranks_rows_1st"))) {
		if(is.null(row_ranks <- methodArgs$row_ranks)) {
			orders <- matrix(factors, nrow(Lambda), factors[1])
			rownames(orders) <- rownames(Lambda)
			colnames(orders) <- colnames(Lambda)
			row_ranks <- FAiR_order_FC(orders, level = 1, MARGIN = 1)
		}
		else methodArgs$row_ranks <- NULL

		if(!identical(dim(row_ranks), dim(Lambda))) {
			stop("'row_ranks' must have the same dimension",
				"as the loading matrix")
		}

		vals <- unique(c(row_ranks))
		if(any(! (vals %in% c(1:factors[1])))) {
			stop("'row_ranks' must only contain positive integers up to ",
				factors[1])
		}

		if(all(row_ranks == factors)) {
			criteria <- criteria[-mark]
			warning("row ranks constraint omitted because no ",
				"constraints were specified")
		}
		else {
			formals(criteria[[mark]])$user_ranks <- row_ranks
			formals(criteria[[mark]])$MARGIN <- 1
		}
	}

	if(length(mark <- which(names(criteria) == "ranks_cols_1st"))) {
		if(is.null(col_ranks <- methodArgs$col_ranks)) {
			orders <- matrix(nrow(Lambda), nrow(Lambda), factors[1])
			rownames(orders) <- rownames(Lambda)
			colnames(orders) <- colnames(Lambda)
			col_ranks <- FAiR_order_FC(orders, level = 1, MARGIN = 2)
		}
		else methodArgs$col_ranks <- NULL

		if(!identical(dim(col_ranks), dim(Lambda))) {
			stop("'col_ranks' must have the same dimension",
				"as the loading matrix")
		}
		vals <- unique(c(col_ranks))
		if(any(! (vals %in% c(1:nrow(Lambda))))) {
			stop("'col_ranks' must only contain positive integers up to ",
				nrow(Lambda))
		}

		if(all(col_ranks == nrow(Lambda))) {
			criteria <- criteria[-mark]
			warning("column ranks constraint omitted because no ",
				"constraints were specified")
		}
		else {
			formals(criteria[[mark]])$user_ranks <- col_ranks
			formals(criteria[[mark]])$MARGIN <- 2
		}
	}

	if(length(mark <- which(names(criteria) == "indicators_1st"))) {
		if(is.null(indicators <- methodArgs$indicators_1st) &&
		   is.null(indicators <- methodArgs$indicators)) {
			indicators <- matrix(FALSE, nrow(Lambda), ncol(Lambda))
			dimnames(indicators) <- dimnames(Lambda)
			indicators <- FAiR_indicators_FC(indicators, level = 1)
		}

		if(any(!(indicators %in% c(1:nrow(Lambda), NA)))) {
			stop("all entries in 'indicators' must be NA or positive",
				" integers less than ", nrow(Lambda))
		}
		formals(criteria[[mark]])$indicators <- indicators
	}

	if(length(mark <- which(names(criteria) == "evRF_1st"))) {
		C <- fitted(FAobject, reduced = FALSE, standardized = TRUE)
		formals(criteria[[mark]])$threshold <- det(C)^(1/ncol(C))
	}

	if(length(mark <- which(names(criteria) == "evPF_1st"))) {
		C <- fitted(FAobject, reduced = FALSE, standardized = TRUE)
		formals(criteria[[mark]])$threshold <- det(C)^(1/ncol(C))
	}

	if(length(mark <- which(names(criteria) == "no_neg_suppressors_1st"))) {
		text <- paste("Select the minimum acceptable factor contribution",
			"coefficient\nsuch that lower values indicate a",
			"suppressor variable.")
		if(is.null(FC_threshold <- methodArgs$FC_threshold)) {
			FC_threshold <- FAiR_get_number(text, from = -1, to = 0,
							by = .001, value = -.01)
		}
		else methodArgs$FC_threshold <- NULL


		if(FC_threshold <= -1) stop("'FC_threshold' must be > -1")
		if(FC_threshold > 0)   stop("'FC_threshold' must be < 0")
		formals(criteria[[mark]])$FC_threshold <- FC_threshold
	}

	mark <- length(criteria)
	if(names(criteria)[mark] == "phi") {
		text <- paste("In calculating phi, the reference structure correlations",
				"will be raised to the power 2/c.\nPlease choose c,",
				"which need not be an integer but is usually taken to be",
				"1.0 by Thurstone.")
		if(is.null(c <- methodArgs$c)) {
			c <- FAiR_get_number(text, from = 1, to = 20, by = .1, value = 1)
		}
		else methodArgs$c <- NULL

		if(c < 1) stop("'c' must be >= 1.0")
		formals(criteria[[mark]])$c <- c
	}
	else if(names(criteria)[mark] == "varphi") {
		if(is.null(weights <- methodArgs$weights)) {
			text  <- paste("What kind of weights would you like to use?")
			items <- c("Dynamic weights", "User-specified weights")
			weights <- FAiR_get_answer(text, radio_items = items)
			if(weights == items[[2]]) {
				weight  <- 1
				weights <- rep(0, factors - 1)
				for(i in 1:length(weights)) {
					weight <- FAiR_get_number(paste("Select weight",
							i), from = 0, to = weight, 
							by = .01, value = 0)
					weights[i] <- weight
				}
			}
			else weights <- -1
		}
		else {
			if(length(weights) != factors - 1) {
				stop("the length of 'weights' must be ", factors - 1)
			}
			else if(any(weights < 0)) {
				stop("all 'weights' must be >= 0")
			}
			else if(any(weights > 1)) {
				stop("all 'weights' must be <= 1")
			}
			else if(!identical(weights, sort(weights, decreasing = TRUE))) {
				stop("'weights' must be in weakly decreasing order")
			}
		}
		formals(criteria[[mark]])$weights <- weights
	}
	else if(names(criteria)[mark] == "LS") {
		if(!is.null(eps <- methodArgs$eps)) {
			formals(criteria[[mark]])$eps <- eps
		}
		else methodArgs$eps <- NULL

		if(!is.null(scale <- methodArgs$scale)) {
			formals(criteria[[mark]])$scale <- scale
		}
		else methodArgs$scale <- NULL

		if(!is.null(E <- methodArgs$E)) {
			formals(critiera[[mark]])$E <- E
		}
		else methodArgs$E <- NULL
	}
	else if(names(criteria)[mark] == "minimaximin") {
		## Do nothing
	}
	else criteria[[mark]] <- FAiR_GPAnator(Lambda, names(criteria)[mark], methodArgs)

	return(criteria)
}

## make criteria with GUI
FAiR_make_criteria_GUI <-
function() {
	CRITERIA <- c(FAiR_constraints(), FAiR_constraints_1st())
	CRITERIA$cohyperplanarity_1st <- CRITERIA$dist_cols_1st <-
		CRITERIA$volume_1st <- CRITERIA$block_1st <- NULL

	criteria <- list()

	# No factor collapse criterion (mandatory)
	criteria[[1]] <- CRITERIA[["no_factor_collapse"]][[1]]
	CRITERIA <- CRITERIA[2:length(CRITERIA)]

	# Other restriction criteria
	items <- sapply(CRITERIA, FUN = function(x) x[[2]])
	checked <- rep(FALSE, length(items))
	text <-  paste( "Please indicate which, if any, additional constraints you",
			"would like to impose.")
	answer <- FAiR_get_answer(text, check_args = list(items = items, 
					checked = checked))
	if(length(answer)) {
		CRITERIA <- sapply(CRITERIA, FUN = function(x) x[[1]])[names(answer)]
		criteria <- c(criteria, CRITERIA)
		names(criteria) <- c("no_factor_collapse", names(answer))
	}
	else    names(criteria) <- "no_factor_collapse"

	# Ultimate analytic criterion
	text   <- paste("Please indicate which analytic criterion should be the ultimate",
			"criterion during the lexical minimization.")
	analytics <- FAiR_analytic()
	items  <- sapply(analytics, FUN = function(x) x[[2]])
	answer <- FAiR_get_answer(text, radio_items = items, select = 1)
	answer <- names(analytics)[items == answer]
	criteria[[length(criteria) + 1]]  <- analytics[[answer]][[1]]
	names(criteria)[length(criteria)] <- answer
	return(criteria)
}
