#############################################################################
#
#  This file is a part of the R package "RoughSets".
#
#  Author: Lala Septem Riza and Andrzej Janusz
#  Supervisors: Chris Cornelis, Francisco Herrera, Dominik Slezak and Jose Manuel Benitez
#  Copyright (c):
#       DiCITS Lab, Sci2s group, DECSAI, University of Granada and
#       Institute of Mathematics, University of Warsaw
#
#  This package is free software: you can redistribute it and/or modify it under
#  the terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 2 of the License, or (at your option) any later version.
#
#  This package is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#############################################################################
# This is a function for implementing the algorithms of quick reduct for feature selection
# based on rough set theory and fuzzy rough set theory.
#
# There are several variety of algorithms based on QuickReduct algorithm (\code{\link{FS.quickreduct.RST}}) and
# (\code{\link{FS.quickreduct.FRST}})
# It should be noted that new decision table is produced based on one chosen reduct.
#
# @title QuickReduct algorithm based on rough set theory and fuzzy rough set theory
#
# @param decision.table a data frame representing decision table. See \code{\link{BC.IND.relation.FRST}}.
# @param type.method a type of the methods
# @param type.QR a type of quick reduct
# @param control a list of other parameters
# @seealso \code{\link{FS.brute.force.RST}}, \code{\link{FS.heuristic.filter.RST}}
# @return reduct.
quickreduct.alg <- function(decision.table, type.method = "fuzzy.dependency", type.QR = "fuzzy.QR", control = list()){
	options(warn=-1)
	## set default values of all parameters
	control <- setDefaultParametersIfMissing(control, list(type.aggregation = c("t.tnorm", "lukasiewicz"), alpha = 0.95,
                                     	t.implicator = "lukasiewicz", type.relation = c("tolerance", "eq.3"), q.some = c(0.1, 0.6), q.most = c(0.2, 1),
										alpha.precision = 0.05, penalty.fact = 0.8, m.owa = round(0.5*nrow(decision.table)),
										k.rfrs = round(0.5*nrow(decision.table)), beta.quasi = 0.05, type.rfrs = "k.trimmed.min", randomize = FALSE))
	## get data
	objects <- decision.table
	desc.attrs <- attr(decision.table, "desc.attrs")
	nominal.att <- attr(decision.table, "nominal.attrs")
	decision.attr <- attr(decision.table, "decision.attr")

	if (is.null(attr(decision.table, "decision.attr"))){
		decision.attr <- ncol(objects)
	} else {
		decision.attr <- attr(decision.table, "decision.attr")
		## check if index of decision attribute is not in last index, then change their position
		if (decision.attr != ncol(objects)){
			objects <- cbind(objects[, -decision.attr, drop = FALSE], objects[, decision.attr, drop = FALSE])
			nominal.att <- c(nominal.att[-decision.attr], nominal.att[decision.attr])
		}
	}

	num.att <- ncol(objects)
	t.implicator <- control$t.implicator
	decision.attr <- control$decision.attr
	type.relation <- control$type.relation
	t.similarity <- type.relation[2]
	type.aggregation <- control$type.aggregation
	type.aggregation <- ch.typeAggregation(type.aggregation)
	t.tnorm <- type.aggregation[[2]]
	t.aggregation <- type.aggregation[[1]]
	alpha <- control$alpha
	q.some <- control$q.some
	q.most <- control$q.most
	m.owa <- control$m.owa
	alpha.precision <- control$alpha.precision
	penalty.fact <- control$penalty.fact
	type.rfrs <- control$type.rfrs
	k.rfrs <- control$k.rfrs
	beta.quasi <- control$beta.quasi
	randomize <- control$randomize

	## get all conditional attributes
    P <- matrix(c(seq(1, (num.att - 1), by = 1)), nrow = 1)
	init.IND <- list()

	## generate indiscernibility relation for each attributes as initial values
	 if (any(type.method == c("fuzzy.dependency", "fuzzy.boundary.reg", "min.positive.reg",
							    "hybrid.rule.induction", "fvprs", "sfrs", "vqrs", "owa", "rfrs", "beta.pfrs"))){
		  control.ind <- list(type.aggregation = type.aggregation, type.relation = type.relation)
		  for (i in 1 : ncol(P)){
			 ## perform indiscernibility relation
			  temp.IND <- BC.IND.relation.FRST(decision.table = decision.table, attributes = c(P[, i]),
								     control = control.ind)

			## Save indiscernibility relation of each attributes
			  init.IND[[i]] <- temp.IND$IND.relation
		  }
		  obj.IND.decAttr <- BC.IND.relation.FRST(decision.table = decision.table, attributes = c(num.att), control = control.ind)
	  }


	## get indiscernibility relation for ALL attribute for stopping criteria
	if (type.method == "quickreduct.RST"){
		## make equivalence class
		IND <- BC.IND.relation.RST(decision.table, feature.set = c(P))

		## get rough set
		roughset <- BC.LU.approximation.RST(decision.table, IND)

		## using def.region.RST to get degree of dependency
		region <- BC.positive.reg.RST(decision.table, roughset)

		## get degree of dependency for all attributes
		degree.all <- region$degree.dependency
	}

	else if (any(type.method == c("vqrs", "hybrid.rule.induction"))){
		## perform fuzzy indiscernibility relation for all attributes by reducing init.IND
		 IND <- Reduce.IND(init.IND, t.tnorm)
		 list.IND <- list(IND.relation = IND, type.relation = type.relation, type.aggregation = type.aggregation, type.model = "FRST")
		 obj.IND <- ObjectFactory(list.IND, classname = "IndiscernibilityRelation")

		 # compute LU approximations
		 if (type.method == "vqrs"){
			control.vqrs <- list(q.some = q.some, q.most = q.most, t.tnorm = t.tnorm)
			FRST <- BC.LU.approximation.FRST(decision.table, obj.IND, obj.IND.decAttr,
								    type.LU = "vqrs", control = control.vqrs)
		 } else {
			control.LU <- list(t.implicator = t.implicator, t.tnorm = t.tnorm)
			FRST <- BC.LU.approximation.FRST(decision.table = decision.table, obj.IND, obj.IND.decAttr,
										 type.LU = "implicator.tnorm", control = control.LU)
		 }
		 ## define fuzzy regions
		 res.temp <- BC.def.region.FRST(decision.table = decision.table, FRST)

		 ## get positive regions and degree of dependency
		 pos.region.all <- res.temp$positive.freg
	}

	## quickreduct based on fuzzy discernibility matrix
	else if (type.method == c("fuzzy.discernibility")){
		## Calculate discernibility matrix based on fuzzy rough set
		IND.cond <- build.discMatrix.FRST(decision.table, type.discernibility = "fuzzy.disc", type.relation = type.relation,
		                               t.implicator = t.implicator, type.LU = "implicator.tnorm", show.discernibilityMatrix = FALSE,
									   type.aggregation = type.aggregation)$IND.conditionAttr
		IND.dec <- build.discMatrix.FRST(decision.table, type.discernibility = "fuzzy.disc", type.relation = type.relation,
		                               t.implicator = t.implicator, type.LU = "implicator.tnorm", show.discernibilityMatrix = FALSE,
									   type.aggregation = type.aggregation)$IND.decisionAttr

		## calculate degree of satisfaction for all condition attributes
		SAT.all <- sum(calc.SAT(IND.cond, IND.dec, attributes = c(P)))
	}

	## initialization
	exit <- FALSE
	cov.rules <- NULL
	cover.indx <- c()
	all.indx.obj <- seq(1, nrow(objects))
	gamma.best <- -1
	gamma.best.bound <- 100000
	gamma.best.level <- 0
	gamma.best.bound.level <- 100000

	while (exit == FALSE){
		gamma.prev <- gamma.best
		gamma.prev.bound <- gamma.best.bound
		gamma.prev.level <- gamma.best.level
		gamma.prev.bound <- gamma.best.bound.level
		reject.attr <- c()

		## calculate degree of dependency
		for (i in 1 : ncol(P)){
			if (type.method == "quickreduct.RST"){
				## make equivalence class
				IND <- BC.IND.relation.RST(decision.table, feature.set = c(P[, i]))

				## get rough set
				roughset <- BC.LU.approximation.RST(decision.table, IND)

				## get degree of dependency
				region <- BC.positive.reg.RST(decision.table, roughset)

				## calculate degree of dependency
				degree <- region$degree.dependency
			}
			## check the type of concept
			else if (any(type.method == c("fuzzy.dependency", "fuzzy.boundary.reg", "min.positive.reg",
			                              "hybrid.rule.induction", "fvprs", "sfrs", "vqrs", "owa", "rfrs", "beta.pfrs"))){

				## calculate and construct the class of indiscernibility relation of each attribute in collection P
				IND <- Reduce.IND(init.IND[c(P[, i])], t.tnorm)
				list.IND <- list(IND.relation = IND, type.relation = type.relation, type.aggregation = type.aggregation, type.model = "FRST")
				obj.IND <- ObjectFactory(list.IND, classname = "IndiscernibilityRelation")

				if (any(type.method == c("fuzzy.dependency", "hybrid.rule.induction", "fuzzy.boundary.reg", "min.positive.reg"))){
					## get lower and upper approximation
					control.LU <- list(t.implicator = t.implicator, t.tnorm = t.tnorm)
					FRST <- BC.LU.approximation.FRST(decision.table = decision.table, obj.IND, obj.IND.decAttr,
													type.LU = "implicator.tnorm", control = control.LU)

					## define fuzzy regions
					fuzzy.region <- BC.def.region.FRST(decision.table = decision.table, FRST)

					if (any(type.method == c("fuzzy.dependency", "hybrid.rule.induction"))){
						## save the degree
						degree <- fuzzy.region$degree.dependency

						## perform rule induction: check positive region of each iteration
						## Based on R. Jensen, C. Cornelis, Q. Shen, "Hybrid fuzzy-rough rule induction and feature selection"
						if (type.method == "hybrid.rule.induction"){
							pos.i <- fuzzy.region$positive.freg
							IND <- obj.IND$IND.relation

							## all.indx.obj is indexes of objects which are considered,
							## initially, it is all indexes of objects
							## cover.indx is indexes which are convered
							cons.indx <- all.indx.obj[all.indx.obj %in% cover.indx == FALSE]

							## loop until all objects are covered
							ii <- 1
							while (ii <= length(cons.indx)){
								## condition when positive region of subset attributes is the same as positive region on all attributes
							   if (pos.i[cons.indx[ii]] == pos.region.all[cons.indx[ii]]){

									## the following step refers to the "check" procedure in the paper
									add <- TRUE

									## it is for first time/initialization
									if (is.null(cov.rules)){

										 ## collect rule
										 rules <- list(attributes = list(colnames(objects[c(P[, i])])), objects = list(cons.indx[ii]))
										 IND.rules <- IND[cons.indx[ii], ,drop = FALSE]
										 cov.rules <- IND[cons.indx[ii], ,drop = FALSE]
										 num.rules <- 1


										 ## set add == FALSE to avoid redudant adding
										 add <- FALSE
									}	else {
										## repeat for each existing rules
										for (j in 1 : num.rules){
											## if existing rule is subset of new rule
											if (all(rules$attributes[[j]] %in% colnames(objects[c(P[, i])])) == TRUE){
												## compare degree of coverage between existing rule and new one

												if (all(IND[cons.indx[ii], ,drop = FALSE] <= cov.rules) == TRUE){
													add <- FALSE
													j <- nrow(IND.rules)
												}	else if (all(IND.rules[j, ,drop = FALSE] < IND[cons.indx[ii], ,drop = FALSE]) == TRUE){
													IND.rules[j, , drop = FALSE] <- NA
													rules$attributes[[j]] <- NULL
													rules$objects[[j]] <- NULL
												}

											}

										}
									}
									## add a new rule and update coverage
									if (add == TRUE){
										rules$attributes <- append(rules$attributes, list(colnames(objects[c(P[, i])])))
										rules$objects <- append(rules$objects, list(cons.indx[ii]))
										IND.rules <- rbind(IND.rules, IND[cons.indx[ii], ,drop = FALSE])
										cov.rules <- apply(rbind(cov.rules, IND[cons.indx[ii], ,drop = FALSE]),
													 2, function(x) max(x))

										cover.indx <- which(cov.rules == 1)

										cons.indx <- all.indx.obj[all.indx.obj %in% cover.indx == FALSE]

										num.rules <- num.rules + 1
										ii <- 0
									}
								}
								ii <- ii + 1
							}
						}
					}
					else if (type.method == c("fuzzy.boundary.reg")){
						## get boundary fuzzy region
						boundary.reg <- fuzzy.region$boundary.freg

						## degree of uncertainty
						degree <- sum(unlist(boundary.reg))/length(boundary.reg)
					}	else if (type.method == c("min.positive.reg")){
						positive.reg <- fuzzy.region$positive.freg
						degree <- min(positive.reg) #/min(pos.region.all)
					}
				}	else if (any(type.method == c("fvprs", "sfrs", "owa", "rfrs", "beta.pfrs"))){
					if (type.method == "fvprs"){
						## get lower and upper approximation
						control.fvprs <- list(t.implicator = t.implicator, t.tnorm = t.tnorm, alpha = alpha.precision)
						FRST <- BC.LU.approximation.FRST(decision.table = decision.table, obj.IND, obj.IND.decAttr,
													type.LU = "fvprs", control = control.fvprs)
					}	else if (type.method == "sfrs"){
						## get lower and upper approximation
						control.sfrs <- list(penalty.fact = penalty.fact)
						FRST <- BC.LU.approximation.FRST(decision.table = decision.table, obj.IND, obj.IND.decAttr,
													type.LU = "sfrs", control = control.sfrs)
					}	else if (type.method == "owa"){
						control.owa <- list(t.implicator = t.implicator, t.tnorm = t.tnorm, m.owa = m.owa)
						FRST <- BC.LU.approximation.FRST(decision.table, obj.IND, obj.IND.decAttr,
							  type.LU = "owa", control = control.owa)
					}	else if (type.method == "rfrs"){
						control.rfrs <- list(t.implicator = t.implicator, t.tnorm = t.tnorm, type.rfrs = type.rfrs, k.rfrs = k.rfrs)
						FRST <- BC.LU.approximation.FRST(decision.table, obj.IND, obj.IND.decAttr,
								  type.LU = "rfrs", control = control.rfrs)
					}	else if (type.method == "beta.pfrs"){
						control.beta.pfrs <- list(t.implicator = t.implicator, t.tnorm = t.tnorm, beta.quasi = beta.quasi)
						FRST <- BC.LU.approximation.FRST(decision.table = decision.table, obj.IND, obj.IND.decAttr,
											 type.LU = "beta.pfrs", control = control.beta.pfrs)
					}
					## define fuzzy regions
					fuzzy.region <- BC.positive.reg.FRST(decision.table = decision.table, FRST)

					## calculate degree of dependency
					degree <- fuzzy.region$degree.dependency
				}	else if (type.method == "vqrs"){
					## perform LU approximation
					control.vqrs <- list(q.some = q.some, q.most = q.most, t.tnorm = t.tnorm)
					FRST.VQRS <- BC.LU.approximation.FRST(decision.table, obj.IND, obj.IND.decAttr,
											  type.LU = "vqrs", control = control.vqrs)

					## determine regions
					pos.region <- BC.positive.reg.FRST(decision.table, FRST.VQRS)$positive.freg
					degree <- min(1, sum(pos.region)/sum(pos.region.all))
				}
			}	else if (type.method == c("fuzzy.discernibility")){
				SAT <- sum(calc.SAT(IND.cond, IND.dec, attributes = c(P[, i])))
				degree <- SAT/SAT.all
			}

			if (type.QR == "modified.QR"){
				if (type.method == "fuzzy.boundary.reg"){
					if ((degree < gamma.best.bound)  || (degree < gamma.prev.bound)){
						reject.attr <- append(reject.attr, i)
					}
				}	else {
					if ((degree < gamma.best.level)  || (degree < gamma.prev.level)){
						reject.attr <- append(reject.attr, i)
					}
				}
			}

			## check to get reduct
			if (type.method == "fuzzy.boundary.reg"){
				if (degree < gamma.best.bound){
					reduct <- P[, i]
					degree.reduct <- degree
					gamma.best.bound <- degree
					indx.min <- i
				}
			}	else {
				if (degree > gamma.best){
					reduct <- P[, i]
					degree.reduct <- degree
					gamma.best <- degree
					indx.max <- i
				}
			}
		}

		########### Procedure to select the attributes ###############
		## for fuzzy.boundary.reg: get min of degree
		if (type.method == "fuzzy.boundary.reg"){
			## check terminating condition
			if (abs(gamma.best.bound - gamma.prev.bound) < 0.00005){
				exit <- TRUE
			}	else if (nrow(P) == (num.att - 1)){
				exit <- TRUE
			}	else {
				if (type.QR == "modified.QR"){
					gamma.best.level <- degree.reduct

					if (!is.null(indx.min)){

						## add attribute which has max degree into P
						## e.g. P <- [1,2,3,4; 2,2,2,2] when indx attr 2 is max
						P <- rbind(P, P[1,indx.min])

						## delete the column of max degree attribute
						P <- P[, -c(indx.min, reject.attr), drop = FALSE]

						## check whether randomize = TRUE means we select attributes randomly
						if (randomize == TRUE){
							P <- randomize.attrs(attributes = P)
						}

						if (ncol(P) <= 1){
							exit <- TRUE
							reduct <- P
						}
					}	else {
						exit <- TRUE
						reduct <- P
					}
					indx.min <- NULL
				}	else {
					if (!is.null(indx.min)){
						## add attribute which has max degree into P
						P <- rbind(P, P[1,indx.min])

						## delete the column of max degree attribute
						P <- P[, -indx.min, drop = FALSE]

						## check whether randomize = TRUE means we select attributes randomly
						if (randomize == TRUE){
							P <- randomize.attrs(attributes = P)
						}

					}
					else {
						exit <- TRUE
					}

					indx.min <- NULL
				}
			}
		}

		## for others: get max of degree
		else if (any(type.method == c("fuzzy.discernibility", "min.positive.reg", "vqrs"))){
			if (degree.reduct >= alpha){
				exit <- TRUE
			} else if (nrow(P) == (num.att - 1)){
				exit <- TRUE
			}	else {
				if (type.QR == "modified.QR"){
					if (!is.null(indx.max)){
						gamma.best.level <- degree.reduct

						## add attribute which has max degree into P
						## e.g. P <- [1,2,3,4; 2,2,2,2] when indx attr 2 is max
						P <- rbind(P, P[1,indx.max])

						## delete the column of max degree attribute
						P <- P[, -c(indx.max, reject.attr), drop = FALSE]

						## check whether randomize = TRUE means we select attributes randomly
						if (randomize == TRUE){
							P <- randomize.attrs(attributes = P)
						}

						if (ncol(P) <= 1){
							exit <- TRUE
							reduct <- P
						}
					}	else {
						exit <- TRUE
						reduct <- P
					}
					indx.max <- NULL
				}	else {
					if (!is.null(indx.max)){
						## add attribute which has max degree into P
						## e.g. P <- [1,2,3,4; 2,2,2,2] when indx attr 2 is max
						P <- rbind(P, P[1,indx.max])

						## delete the column of max degree attribute
						P <- P[, -indx.max, drop = FALSE]

						## check whether randomize = TRUE means we select attributes randomly
						if (randomize == TRUE){
							P <- randomize.attrs(attributes = P)
						}
					}	else {
						exit <- TRUE
					}
					indx.max <- NULL
				}
			}
		}	else if (any(type.method == c("quickreduct.RST"))){
			if (degree.reduct == degree.all){
				exit <- TRUE
			}	else {
				## add attribute which has max degree into P
				## e.g. P <- [1,2,3,4; 2,2,2,2] when indx attr 2 is max
				P <- rbind(P, P[1,indx.max])

				## delete the column of max degree attribute
				P <- P[, -indx.max, drop = FALSE]

				## check whether randomize = TRUE means we select attributes randomly
				if (randomize == TRUE){
					P <- randomize.attrs(attributes = P)
				}
			}
		}	else if (any(type.method == c("fuzzy.dependency", "owa", "rfrs", "fvprs", "sfrs", "beta.pfrs", "hybrid.rule.induction"))){
			if (abs(gamma.best - gamma.prev) < 0.000005){
				exit <- TRUE
			}	else if (nrow(P) == (num.att - 1)){
				exit <- TRUE
			}	else {
				if (type.QR == "modified.QR"){
					if (!is.null(indx.max)){
						gamma.best.level <- degree.reduct

						## add attribute which has max degree into P
						## e.g. P <- [1,2,3,4; 2,2,2,2] when indx attr 2 is max
						P <- rbind(P, P[1,indx.max])

						## delete the column of max degree attribute
						P <- P[, -c(indx.max, reject.attr), drop = FALSE]

						## check whether randomize = TRUE means we select attributes randomly
						if (randomize == TRUE){
							P <- randomize.attrs(attributes = P)
						}

						if (ncol(P) <= 1){
							exit <- TRUE
							reduct <- P
						}
					}	else {
						exit <- TRUE
					}
					indx.max <- NULL
				}	else {
					if (!is.null(indx.max)){
						## add attribute which has max degree into P
						## e.g. P <- [1,2,3,4; 2,2,2,2] when indx attr 2 is max
						P <- rbind(P, P[1,indx.max])

						## delete the column of max degree attribute
						P <- P[, -indx.max, drop = FALSE]

						## check whether randomize = TRUE means we select attributes randomly
						if (randomize == TRUE){
							P <- randomize.attrs(attributes = P)
						}
					}	else {
						exit <- TRUE
					}
					indx.max <- NULL
				}
			}
		}
	}

	if (type.method == "hybrid.rule.induction"){
		## change into standard form of rules and then produce result
		return(std.rules(rules, type.method = "RI.hybridFS.FRST", decision.table, t.similarity, t.tnorm))
	}	else {
		reduct <- sort(reduct)
		return(reduct)
	}
}

# this function is used to calculate degree of satisfaction
# @param IND.cond indiscernibility matrix for all conditional attributes
# @param IND.dec indiscernibility matrix for decision attribute
# @param attributes a list of considered attributes
calc.SAT <- function(IND.cond, IND.dec, attributes){

	SAT	<- matrix()
	ii <- 1
	for (j in 1 : nrow(IND.cond)){
		for (k in 1 : ncol(IND.cond)){
			if (j < k) {
				temp <- matrix()
				for (i in 1 : length(attributes)){
						temp[i] <- unlist(IND.cond[j, k])[attributes[i]]
				}
				if (length(attributes) > 1){
					val <- calc.conorm(func.conorm, data = temp[-1], init.val = temp[1], t.conorm = "lukasiewicz")
				}	else {
					val <- temp
				}
				SAT[ii] <- calc.implFunc(IND.dec[j, k], val, type.implication.func = "lukasiewicz")
				ii <- ii + 1
			}
		}
	}

	return(SAT)
}

# This function is used to calculate conorm all relations
#
# @title conorm calculation
# @param func.conorm a function to calculate 2 variable using a certain conorm
# @param data a vector of data
# @param init.val a value of attribute
# @param t.conorm a type of t-conorm
calc.conorm <- function(func.conorm, data, init.val, t.conorm = "lukasiewicz", ...){
	n <- length(data)
	conorm.val <- func.conorm(data[1], init.val, t.conorm)
	if (n > 1) {
		for (i in rev(data[-1])) {
			conorm.val <- func.conorm(i, conorm.val, t.conorm)
		}
	}
	return(conorm.val)

}

# This function is used to calculate conorm of two variables
#
# @title conorm calculation on two variables
# @param right.val a value of one attribute on right side
# @param init.val a value of attribute
# @param cotnorm a type of t-norm
func.conorm <- function(right.val, init.val, t.conorm){
	if (t.conorm == "lukasiewicz"){
		return (min(right.val + init.val, 1))
	}	else if (t.conorm == "max"){
		return (min(right.val, init.val))
	}
}

#' An auxiliary function for the \code{qualityF} parameter in the \code{\link{FS.greedy.heuristic.reduct.RST}}, \code{\link{FS.DAAR.heuristic.RST}} and \code{\link{FS.greedy.heuristic.superreduct.RST}} functions.
#' It is based on the \emph{gini} index as a measure of information (Stoffel and Raileanu, 2000).
#'
#' @title The gini-index measure
#' @author Andrzej Janusz
#'
#' @param decisionDistrib an integer vector corresponding to a distribution of attribute values.
#'
#' @return a numeric value indicating the gini index of an attribute.
#' @references
#' K. Stoffel and L. E. Raileanu, "Selecting Optimal Split-Functions with Linear Threshold Unit Trees and Madaline-Style Networks",
#' in: Research and Development in Intelligent Systems XVII, BCS Conference Series (2000).
#'
#' @export
X.gini <- function(decisionDistrib)  {
  return(1 - sum((decisionDistrib/sum(decisionDistrib))^2))
}

#' An auxiliary function for the \code{qualityF} parameter in the \code{\link{FS.greedy.heuristic.reduct.RST}}, \code{\link{FS.DAAR.heuristic.RST}} and \code{\link{FS.greedy.heuristic.superreduct.RST}} functions.
#' It is based on \emph{entropy} as a measure of information(Shannon, 1948).
#'
#' @title The entropy measure
#' @author Andrzej Janusz
#'
#' @param decisionDistrib an integer vector corresponding to a distribution of attribute values
#'
#' @return a numeric value indicating entropy of an attribute.
#' @references
#' C. E. Shannon, "A Mathematical Theory of Communication", Bell System Technical Journal, vol. 27, p. 379 - 423, 623 - 656 (1948).
#'
#' @export
X.entropy <- function(decisionDistrib)  {
  probabilities = decisionDistrib/sum(decisionDistrib)
  return(-sum(probabilities*log2(probabilities), na.rm=T))
}

#' It is an auxiliary function for the \code{qualityF} parameter in the \code{\link{FS.greedy.heuristic.reduct.RST}} and \code{\link{FS.greedy.heuristic.superreduct.RST}} functions.
#'
#' @title  The discernibility measure
#' @author Andrzej Janusz
#'
#' @param decisionDistrib an integer vector corresponding to a distribution of decision attribute values
#' @return a numeric value indicating a number of conflicts in a decision attribute
#'
#' @export
X.nOfConflicts <- function(decisionDistrib)  {
  return(as.numeric(sum(as.numeric(sum(decisionDistrib) - decisionDistrib) * as.numeric(decisionDistrib))))
}

# It is an auxiliary function for the \code{qualityF} parameter in the \code{\link{FS.greedy.heuristic.reduct.RST}} and \code{\link{FS.greedy.heuristic.superreduct.RST}} functions.
#
# @title The discernibility measure function based on \code{log2}
# @param decisionDistrib a distribution of values on the decision attribute.
# @return \code{qualityF}.
X.nOfConflictsLog <- function(decisionDistrib)  {
  return(log2(1+as.numeric(sum(as.numeric(sum(decisionDistrib) - decisionDistrib) * decisionDistrib))))
}

# It is an auxiliary function for the \code{qualityF} parameter in the \code{\link{FS.greedy.heuristic.reduct.RST}} and \code{\link{FS.greedy.heuristic.superreduct.RST}} functions.
#
# @title The discernibility measure function based on \code{sqrt}
# @param decisionDistrib a distribution of values on the decision attribute.
# @return \code{qualityF}.
X.nOfConflictsSqrt <- function(decisionDistrib)  {
  return(sqrt(as.numeric(sum(as.numeric(sum(decisionDistrib) - decisionDistrib) * as.numeric(decisionDistrib)))))
}

# This is an auxiliary function for computing decision reducts
# qualityGain2 <- function(vec, uniqueValues, decisionVec, uniqueDecisions,
#                         INDclasses, INDclassesSizes, baseChaos, chaosFunction = X.gini)  {
#
#   classCountsList = compute_indiscernibility_and_chaos2(INDclasses, vec, uniqueValues,
#                                                         decisionVec, uniqueDecisions)
#   INDclassesLengths = classCountsList[[1]]
#   tmpInd = which(INDclassesLengths > 1)
#   if(length(tmpInd) > 0) {
#     remainingChaos = sapply(classCountsList[tmpInd + 1], chaosFunction)
#     remainingChaos = sum(remainingChaos * INDclassesLengths[tmpInd]) / length(decisionVec)
#   } else remainingChaos = 0
#
#   return(as.numeric(baseChaos - remainingChaos))
# }

qualityGain <- function(vec, uniqueValues, decisionVec, uniqueDecisions,
                        INDclassesList, INDclassesSizes, baseChaos, chaosFunction = X.gini)  {

  classCounts = .C("computeIndiscernibilityAndChaos",
                   INDclasses = as.integer(unlist(INDclassesList)),
                   INDsizes = as.integer(INDclassesSizes),
                   NOfINDClasses = as.integer(length(INDclassesSizes)),
                   attrValues = as.integer(vec),
                   NOfAttrValues = as.integer(length(uniqueValues)),
                   decValues = as.integer(decisionVec),
                   NOfDecs = as.integer(length(uniqueDecisions)),
                   output = as.integer(rep(0,length(INDclassesList)*length(uniqueValues)*length(uniqueDecisions))))

  classCounts = matrix(classCounts$output,
                       nrow = length(INDclassesList)*length(uniqueValues),
                       ncol = length(uniqueDecisions), byrow=TRUE)
  newINDsizes = rowSums(classCounts)
  tmpInd = newINDsizes > 1
  if(sum(tmpInd) > 0) {
    remainingChaos = apply(classCounts[tmpInd, ,drop=FALSE], 1, chaosFunction)
    remainingChaos = sum(remainingChaos * newINDsizes[tmpInd]) / length(decisionVec)
  } else remainingChaos = 0

  return(as.numeric(baseChaos - remainingChaos))
}


# It is used to randomize attributes
# @param attributes a matrix of attributes
randomize.attrs <- function(attributes){
	## check whether randomize = TRUE means we select attributes randomly
	rand.num <- sample(0:1, ncol(attributes), replace = TRUE)
	if (sum(rand.num) == 0){
		attributes <- attributes[, 1, drop = FALSE]
	}	else {
		attributes <- attributes[, which(rand.num == 1), drop = FALSE]
	}
	return (attributes)
}

# This function is used to calculate a discernibility function used to get all reducts
# @param discernibilityMatrix an object of a class "DiscernibilityMatrix
calc.all.reducts <- function(discernibilityMatrix) {

  disc.list = discernibilityMatrix$disc.list
  disc.mat = discernibilityMatrix$disc.mat
  type.discernibility = discernibilityMatrix$type.discernibility
  type.model = discernibilityMatrix$type.model
  names.attr = discernibilityMatrix$names.attr

  ## delete redundant of discernibility matrix
  disc.list <- unique(disc.list)

  ## delete blank characters
  disc.list.temp <- c()
  j <- 1
  for (i in 1 : length(disc.list)){
    if (!(is.character(disc.list[[i]]) & length(disc.list[[i]]) == 0)){
      disc.list.temp[j] <- list(disc.list[[i]])
      j <- j + 1
    }
  }
  disc.list <- disc.list.temp

  ## perform boolean function
  res <- boolean.func(disc.list)

  ## results
  if (length(res$dec.red) == 0){
    ## construct DecisionReduct class
    mod <- list(discernibility.matrix = disc.mat, decision.reduct = NULL, core = NULL, discernibility.type = type.discernibility,
                type.task = "computation of all reducts", type.model = type.model)
    class.mod <- ObjectFactory(mod, classname = "ReductSet")

    return(class.mod)
  }
  else if (length(res$dec.red) == 1){

	## construct DecisionReduct class
    reduct = which(as.character(names.attr) %in% res$dec.red)
    names(reduct) = names.attr[reduct]

    mod <- list(discernibility.matrix = disc.mat, decision.reduct = list(reduct), core = reduct, discernibility.type = type.discernibility,
                type.task = "computation of all reducts", type.model = type.model)
    class.mod <- ObjectFactory(mod, classname = "ReductSet")

    return(class.mod)
  }
  else {
    ## transform CNF into DNF
    CNFclauses = res$dec.red
    CNFlengths = sapply(CNFclauses, length)
    CNFclauses = CNFclauses[order(CNFlengths)]
    tmpDNF = CNFclauses[[1]]

	for(i in 2:length(CNFclauses))  {
      tmpDNF = expand.grid(tmpDNF, CNFclauses[[i]], KEEP.OUT.ATTRS = F, stringsAsFactors = F)
      tmpDNF = split(tmpDNF, 1:nrow(tmpDNF))
      tmpDNF = lapply(tmpDNF, function(x) unique(unlist(x)))
      tmpLengths = sapply(tmpDNF, length)
      tmpDNF = tmpDNF[order(tmpLengths)]
      tmpDNFlength = length(tmpDNF)
      j = 2
      while(j <= tmpDNFlength)  {
        tmpIdx = sapply(tmpDNF[j:tmpDNFlength], function(x,y) all(y %in% x), tmpDNF[[j-1]])
        tmpDNF = tmpDNF[c(rep(T, j-1), !tmpIdx)]
        j = j + 1
        tmpDNFlength = length(tmpDNF)
      }
    }

    DNFclauses = lapply(tmpDNF, function(x) x[order(x)])

    if (length(DNFclauses) == 1){
      core <- DNFclauses
    }
	else {
      core <- res$core
    }

    if(length(core) > 0) {
      core = which(as.character(names.attr) %in% core)
      names(core) = names.attr[core]
    }

    DNFclauses = lapply(DNFclauses, function(x,namesVec) {x = which(namesVec %in% x);
                                                          names(x) = namesVec[x];
                                                          return(x)}, as.character(names.attr))

    mod <- list(decision.reduct = DNFclauses, core = core, discernibility.type = type.discernibility,
                type.task = "computation of all reducts", type.model = type.model)
    class.mod <- ObjectFactory(mod, classname = "ReductSet")

    return(class.mod)
  }
}

# it is used to perform boolean function
# @param disc.list a list of attributes resulted by build.discMatrix.FRST
boolean.func <- function(disc.list){
#	req.suc <- require("sets", quietly=TRUE)
#	if(!req.suc) stop("In order to use this function, you need to install the package sets.")

	## check subset among objects
	if (length(disc.list) > 1){
		## create a function to be vectorize
		func.ch.duplicate <- function(i, j, disc.list){
			if (j > i){
#				if (set_is_subset(as.set(disc.list[[i]]), as.set(disc.list[[j]]))) {
			  if (all(disc.list[[i]] %in% disc.list[[j]])) {
					disc.list[[j]] <<- NA
				}
#				else if (set_is_subset(as.set(disc.list[[j]]), as.set(disc.list[[i]]))) {
        else if (all(disc.list[[j]] %in% disc.list[[i]])) {
					disc.list[[i]] <<- NA
				}
			}
		}

		## vectorize func.ch.duplicate
		VecFun <- Vectorize(func.ch.duplicate, vectorize.args=list("i","j"))
		outer(1 : length(disc.list), 1 : length(disc.list), VecFun, disc.list)
	}

	## filter the list
	## collect decision reduct and core
	dec.red <- list()
	core <- c()
	for (i in 1 : length(disc.list)){
		if (!is.na(unlist(disc.list[[i]][1]))){
			dec.red <- append(dec.red, list(disc.list[[i]]))
			if (length(disc.list[[i]]) == 1){
				core <- union(core, disc.list[[i]])
			}
		}
	}

	return(list(dec.red = dec.red, core = core))
}

# a function for converting formulas in a CNF form to a DNF form
convertCNFtoDNF <- function(CNFclauses){

  if(length(CNFclauses) > 1) {
    ## sort the clauses by their length
    CNFlengths = sapply(CNFclauses, length)
    if(any(CNFlengths == 0)) {
      CNFclauses = CNFclauses[CNFlengths > 0]
      CNFlengths = CNFlengths[CNFlengths > 0]
    }
    CNFclauses = CNFclauses[order(CNFlengths)]

    ## eliminate unnecessary clauses for efficiency
    tmpCNFlength = length(CNFclauses)
    j = 2
    while(j <= tmpCNFlength){
      tmpIdx = sapply(CNFclauses[j:tmpCNFlength], function(x,y) all(y %in% x), CNFclauses[[j-1]])
      CNFclauses = CNFclauses[c(rep(T, j-1), !tmpIdx)]
      j = j + 1
      tmpCNFlength = length(CNFclauses)
    }

    ## convert to DNF form - start from the first CNF clause
    tmpDNF = CNFclauses[[1]]
    for(i in 2:length(CNFclauses))  {
      ## expand two clauses into possible DNFs
      tmpDNF = expand.grid(tmpDNF, CNFclauses[[i]],
                           KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      tmpDNF = split(tmpDNF, 1:nrow(tmpDNF))
      ## take only those which are unique
      tmpDNF = lapply(tmpDNF, function(x) unique(unlist(x)))
      tmpLengths = sapply(tmpDNF, length)
      tmpDNF = tmpDNF[order(tmpLengths)]
      ## eliminate unnecessary clauses for efficiency
      tmpDNFlength = length(tmpDNF)
      j = 2
      while(j <= tmpDNFlength)  {
        tmpIdx = sapply(tmpDNF[j:tmpDNFlength], function(x,y) all(y %in% x), tmpDNF[[j-1]])
        tmpDNF = tmpDNF[c(rep(T, j-1), !tmpIdx)]
        j = j + 1
        tmpDNFlength = length(tmpDNF)
      }
    }
  } else {
    tmpDNF = CNFclauses
  }

  DNFclauses = lapply(tmpDNF, function(x) x[order(x)])
  names(DNFclauses) = paste("reduct", 1:length(DNFclauses), sep = "")
  return(DNFclauses)
}

# a function for computing a core from a list of all reducts of a data set
computeCore = function(reductList) {

  if(length(reductList) > 0) {
    stopFlag = FALSE
    i = 2
    core = reductList[[1]]
    while(!stopFlag && i < length(reductList)) {
      core = core[core %in% reductList[[i]]]
      i = i + 1
      if(length(core) < 1) stopFlag = TRUE
    }
  } else stop("empty reduct list")

  return(core)
}

# An auxiliary function for computing attribute relevance using random probes.
# It is used by FS.DAAR.heuristic.RST function.
# computeRelevanceProb2 = function(INDclasses, INDclassesSizes, attributeVec, uniqueValues,
#                                 attrScore, decisionVec, uniqueDecisions, baseChaos,
#                                 qualityF = X.gini, nOfProbes = 100, withinINDclasses = FALSE)
# {
#   if(withinINDclasses) {
#     attributeList = lapply(INDclasses, function(x, y) y[x], attributeVec)
#     tmpAttribute = attributeVec
#     flattenINDclasses = unlist(INDclasses, use.names = FALSE)
#     probeScores = replicate(nOfProbes,
#                             {tmpAttributePermutation = unlist(lapply(attributeList, sample),
#                                                               use.names = FALSE);
#                              tmpAttribute[flattenINDclasses] = tmpAttributePermutation;
#                              qualityGain(tmpAttribute, uniqueValues, decisionVec, uniqueDecisions,
#                                          INDclasses, INDclassesSizes, baseChaos, chaosFunction = qualityF)})
#   } else {
#     probeScores = replicate(nOfProbes,
#                             {tmpAttribute = sample(attributeVec);
#                              qualityGain(tmpAttribute, uniqueValues, decisionVec, uniqueDecisions,
#                                          INDclasses, INDclassesSizes, baseChaos, chaosFunction = qualityF)})
#   }
#
# return(mean(attrScore > probeScores))
# }

# An auxiliary function for computing attribute relevance using random probes.
# It is used by FS.DAAR.heuristic.RST function.
computeRelevanceProb = function(INDclasses, INDclassesSizes, attributeVec, uniqueValues,
                                attrScore, decisionVec, uniqueDecisions, baseChaos,
                                qualityF = X.gini, nOfProbes = 100, withinINDclasses = FALSE)
{
  if(withinINDclasses) {
    attributeList = lapply(INDclasses, function(x, y) y[x], attributeVec)
    flattenINDclasses = unlist(INDclasses, use.names = FALSE)
    probesList = replicate(nOfProbes,
                           {tmpAttribute = attributeVec;
                            tmpAttributePermutation = unlist(lapply(attributeList, sample),
                                                             use.names = FALSE);
                            tmpAttribute[flattenINDclasses] = tmpAttributePermutation;
                            tmpAttribute},
                           simplify = FALSE)

    probeScores = sapply(probesList, qualityGain,
                         uniqueValues = uniqueValues,
                         decisionVec = decisionVec,
                         uniqueDecisions = uniqueDecisions,
                         INDclassesList = INDclasses,
                         INDclassesSizes = INDclassesSizes,
                         baseChaos = baseChaos,
                         chaosFunction = qualityF,
                         USE.NAMES = FALSE)
    rm(probesList, flattenINDclasses, attributeList)
  } else {
    probeScores = replicate(nOfProbes,
                            {tmpAttribute = sample(attributeVec);
                            qualityGain(tmpAttribute, uniqueValues, decisionVec, uniqueDecisions,
                                        INDclasses, INDclassesSizes, baseChaos, chaosFunction = qualityF)})
  }

  return(mean(attrScore > probeScores))
}
