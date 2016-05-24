#############################################################################
#
#  This file is a part of the R package "frbs".
#
#  Author: Lala Septem Riza
#  Co-author: Christoph Bergmeir
#  Supervisors: Francisco Herrera Triguero and Jose Manuel Benitez
#  Copyright (c) DiCITS Lab, Sci2s group, DECSAI, University of Granada.
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
#' This is the internal function that implements genetic fuzzy systems for fuzzy rule learning 
#' based on the MOGUL methodology (GFS.FR.MOGUL). It is used to solve regression tasks. 
#' Users do not need to call it directly, but just use \code{\link{frbs.learn}} and \code{\link{predict}}.
#' 
#' This method was proposed by Herrera et al.   
#' GFS.FR.MOGUL implements a genetic algorithm determining the structure 
#' of the fuzzy IF-THEN rules and the membership function parameters. 
#' There are two general types of fuzzy IF-THEN rules, namely the descriptive 
#' and the approximative/free semantic approaches. A descriptive approach means that
#' the linguistic labels represent a real-world semantic; the linguistic labels are 
#' uniformly defined for all rules. In contrast, in the approximative approach 
#' there isn't any associated linguistic label. This method is based on the latter one. 
#' We model a fuzzy IF-THEN rule on a chromosome which consists of the parameter values of 
#' the membership function. So, every rule has its own membership function values. 
#' A population contains many such generated chromosomes, based on the iterative rule learning 
#' approach (IRL). IRL means that the chromosomes will be generated one by one, taking into
#' account the fitness value and covering factor, until there are sufficient chromosomes 
#' in the population. After having obtained the population, the genetic algorithm is started,
#' using the genetic operators selection, mutation, and crossover. 
#'  
#' @title GFS.FR.MOGUL model building 
#'
#' @param data.train a matrix (\eqn{m \times n}) of normalized data for the training process, where \eqn{m} is the number of instances and 
#'        \eqn{n} is the number of variables; the last column is the output variable. Note the data must be normalized between 0 and 1. 
#' @param persen_cross a real number between 0 and 1 determining the probability of crossover.
#' @param persen_mutant a real number between 0 and 1 determining the probability of mutation.
#' @param max.iter the maximal number of iterations.
#' @param max.gen the maximal number of generations of the genetic algorithm.
#' @param max.tune the maximal number of tuning iterations.
#' @param range.data.ori a matrix containing the ranges of the original data. 
#' @param epsilon a real number between 0 and 1 determining the boundary of covering factor. 
#' @seealso \code{\link{GFS.FR.MOGUL.test}}, \code{\link{frbs.learn}}, and \code{\link{predict}}
# @return list of the model data. please have a look at \code{\link{frbs.learn}} for looking complete components.
#' @references
#' F. Herrera, M. Lozano, and J.L. Verdegay, 
#' "A learning process for fuzzy control rules using genetic algorithms", 
#' Fuzzy Sets and Systems, vol. 100, pp. 143 - 158 (1998).
#'
#' O. Cordon, M.J. del Jesus, F. Herrera, and M. Lozano, 
#' "MOGUL: A methodology to obtain genetic fuzzy rule-based systems 
#' under the iterative rule learning approach", International Journal of Intelligent Systems, 
#' vol. 14, pp. 1123 - 1153 (1999).
# @export
GFS.FR.MOGUL <- function(data.train, persen_cross = 0.6, persen_mutant = 0.3, max.iter = 10, max.gen = 10, max.tune = 10,
                         range.data.ori, epsilon = 0.4){
	## A genetic generation stage
	data.train.ori <- data.train
	data.sample <- data.train
	num.var <- ncol(data.train)
	num.inputvar <- num.var - 1
	best.rules <- matrix(NA, nrow = 1, ncol = (3 * ncol(data.sample)))
	iter <- 1
	stop.cr <- FALSE
	num.rules.popu <- 10
	
	## create progress bar
	progressbar <- txtProgressBar(min = 0, max = max.iter, style = 3)	
	while (stop.cr == FALSE){	
		## get rule based on data.train 
		temp.rule <- matrix(NA, nrow = 1, ncol = (3 * ncol(data.sample)))
		
		## check number of rules
		if (nrow(data.sample) < num.rules.popu){
			num.rules.popu <- nrow(data.sample)
		}
		
		## generate initial population
		for (i in 1 : num.rules.popu){
			rule <- generate.popu(method.type = "GFS.FR.MOGUL", data.train = data.sample[i, ,drop = FALSE])
			temp.rule <- rbind(temp.rule, rule)
		}
		popu <- na.omit(temp.rule)		
		popu.size <- nrow(popu)
		old.best.fit <- -1
		gen <- 1
		
		while (gen < max.gen){
			##calculate individual fitness
			ind.fit <- eval.indv.fitness(data.train = data.sample, rule = popu, method.type = "GFS.FR.MOGUL", epsilon = epsilon)

			## combine rule.gen and its fitness
			rule.withFit <- cbind(popu, ind.fit)

			## get the best individual
			new.best.fit <- max(rule.withFit[, ncol(rule.withFit)])
			if (new.best.fit > old.best.fit){
				indx <- which.max(rule.withFit[, ncol(rule.withFit), drop = FALSE])
				new.best.rule <- rule.withFit[indx, -ncol(rule.withFit), drop = FALSE]
				old.best.fit <- new.best.fit
			}
			
			popu.rule <- rule.withFit[, -ncol(rule.withFit), drop = FALSE]
			if (nrow(popu.rule) > 2){
				## GA.Crossover
				popu <- GA.crossover(popu.rule, persen_cross, type = "MOGUL", num.inputvar)
			}
			## GA.Mutation
			popu <- GA.mutation(popu, persen_mutant, type = "MOGUL", iter = gen, 
			                 max.iter = max.gen)			
						
			gen = gen + 1
		}		
		best.rules <- rbind(best.rules, new.best.rule)
				
		data.sample <- red.data(data.sample, new.best.rule, epsilon, type = "MOGUL")
		iter = iter + 1	
		
		## update the number of data after reducing	
		num.sample.fin <- nrow(data.sample)	
		if (num.sample.fin <= 3 || iter >= max.iter){
			stop.cr <- TRUE
			iter <- max.iter
		}
			
		## progress bar
		setTxtProgressBar(progressbar, iter)
	}
	close(progressbar)
	best.rules <- na.omit(best.rules)

	mod <- list(rule = best.rules)
	
	params <- list(mod = mod, data.train = data.train.ori, max.iter = max.tune)
	new.mod <- tune.MF(method.type = "GFS.FR.MOGUL", params)
	best.rules <- new.mod$rule
	
	mod <- list(rule = best.rules, range.data.ori = range.data.ori, type.tnorm = "MIN", 
	   type.snorm = "SUM", type.model = "APPROXIMATE", type.implication.func = "MIN")	
	return(mod)
}

#' This is the internal function that implements the Thrift's technique based 
#' on a genetic algorithm. It is used to tackle regression tasks. Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}.
#'
#' This method was developed by Thrift using Mamdani's model as fuzzy IF-THEN rules. 
#' In this method, we consider a table as a genotype with alleles that are fuzzy set indicators 
#' over the output domain. The phenotype is produced by the behavior produced 
#' by the fuzzification, max-* composition, and defuzzification operations. 
#' A chromosome (genotype) is formed from the decision table by going rowwise 
#' and producing a string of numbers from the code set. Standard crossover 
#' and mutation operators can act on these string. 
#' 
#' @title GFS.Thrift model building 
#' 
#' @param data.train a matrix (\eqn{m \times n}) of normalized data for the training process, 
#'        where \eqn{m} is the number of instances and 
#'        \eqn{n} is the number of variables; the last column is the output variable.
#'        Note the data must be normalized between 0 and 1. 
#' @param popu.size the size of the population which is generated in each generation.
#' @param num.labels a matrix describing the number of linguistic terms.
#' @param persen_cross a real number between 0 and 1 representing the probability of crossover.
#' @param persen_mutant a real number between 0 and 1 representing the probability of mutation.
#' @param max.gen the maximal number of generations for the genetic algorithm.
#' @param range.data.ori a matrix containing the ranges of the original data. 
#' @param type.defuz the type of the defuzzification method. For more detail, see \code{\link{defuzzifier}}.
#'        The default value is \code{WAM}.
#' @param type.tnorm the type of t-norm. For more detail, please have a look at \code{\link{inference}}. 
#' @param type.snorm the type of s-norm. For more detail, please have a look at \code{\link{inference}}. 
#' @param type.mf the type of shape of membership function. See \code{\link{fuzzifier}}. 
#' @param type.implication.func the type of implication function. See \code{\link{WM}}.
#' @seealso \code{\link{GFS.Thrift.test}}, \code{\link{frbs.learn}}, and \code{\link{predict}}
# @return list of the model data. please have a look at \code{\link{frbs.learn}} for looking complete components.
#' @references
#' P. Thrift, "Fuzzy logic synthesis with genetic algorithms", In Proceedings of the Fourth
#' International Conference on Genetic Algorithms (ICGA91), San Diego (United States of
#' America), pp. 509 - 513 (1991).
# @export
GFS.Thrift <- function(data.train, popu.size = 10, num.labels, persen_cross = 0.6, persen_mutant = 0.3,
 max.gen = 10, range.data.ori, type.defuz = "WAM", type.tnorm = "MIN", type.snorm = "MAX", type.mf = "TRIANGLE", type.implication.func = "ZADEH"){
		
	## create progress bar
	progressbar <- txtProgressBar(min = 0, max = max.gen, style = 3)		
	range.data <- matrix(nrow = 2, ncol = ncol(data.train))
	range.data[1, ] <- 0
	range.data[2, ] <- 1
	
	## initialize mod
	mod <- NULL

	#### Initialize population
	num.var <- ncol(data.train)
	num.inputvar <- num.var - 1
	if (popu.size > num.labels[1,1]^num.var){
		popu.size = num.labels[1,1]^num.var
	}
	
	## generate initial population
	complete.popu <- matrix(nrow = (num.labels[1,1]^num.var), ncol = num.var)
	x <- list(1:num.labels[1,1])
	tmp = expand.grid(rep(x, num.var))
	complete.popu <- as.matrix(tmp)
	rule.data.num  <- complete.popu[sample(nrow(complete.popu)),]	
	
	## change scale terms
	rule.data.num <- scale.StrToRule(rule.data.num, num.labels)
	rule.tmp <- generate.rule(rule.data.num, num.labels)	
	var.mf <- partition.MF(range.data, num.labels, type.mf)
	
	## Construct parameters
	mod$range.input <- range.data[, 1 : (ncol(range.data) - 1), drop = FALSE]
	mod$range.output <- range.data[, ncol(range.data), drop = FALSE]
	mod$num.varinput <- (ncol(data.train) - 1)
	mod$num.fvalinput <- num.labels[1, 1 : (ncol(num.labels) - 1), drop = FALSE]
	mod$names.varinput <- rule.tmp$names.varinput
	mod$varinp.mf <- var.mf[, 1 : (num.labels[1,1] * (num.var - 1)), drop = FALSE]
	colnames(mod$varinp.mf) <- mod$names.varinput
	mod$names.varoutput <- rule.tmp$names.varoutput
	mod$rule <- rule.tmp$rule
	mod$rule.data.num = rule.data.num
	mod$varout.mf <- var.mf[, (num.labels[1,1] * (num.var - 1) + 1) : ncol(var.mf), drop = FALSE]
	colnames(mod$varout.mf) <- mod$names.varoutput
	mod$type.defuz <- type.defuz
	mod$type.tnorm <- type.tnorm
	mod$type.snorm <- type.snorm
	mod$type.model <- "MAMDANI"
	mod$method.type <- "GFS.THRIFT"
	
	## calculate cover factor
	heu.val <- calc.heu(mod, data.train, num.labels, type.implication.func)
	tmp <- cbind(rule.data.num, heu.val)
	
	## sort based on cover factor
	rule.data.num.sorted <- tmp[order(tmp[, ncol(tmp)], decreasing = TRUE), ]
	
	## get rule.data.num 
	rule.data.num.sorted <- rule.data.num.sorted[1 : popu.size, ,drop = FALSE]
	rule.data.num <- rule.data.num.sorted[, 1 : (ncol(rule.data.num.sorted) - 1), drop = FALSE]
	
	## make testing data from training data
	data.test <- data.train[, 1 : (ncol(data.train) - 1), drop = FALSE]
	real.val <- data.train[, ncol(data.train), drop = FALSE]

	## initialize for looping
	stop <- 0	
	best.fit <- 1000000000
	temp.fit <- matrix()
	bench <- matrix(ncol = 2)
	
	while (stop <= max.gen){	

		## delete or change unnecessary rules 
		res.rule.data.num <- prune.rule(rule.data.num = rule.data.num, method = "THRIFT", num.labels = num.labels)
		mod$rule.data.num <- res.rule.data.num
		pred.val <-	GFS.Thrift.test(mod, data.test)
		
		bench <- cbind(real.val, pred.val)
		
		## calculate RMSE as fitness function
		residuals <- matrix(bench[, 1] - bench[, 2], nrow = 1)
		RMSE <- sqrt(mean(residuals^2))		
		fitness	<- RMSE	
		if (fitness < best.fit){
			best.fit <- fitness
			best.rule.data.num <- mod$rule.data.num	
		}		
				
		### crossover
		rule.data.num <- GA.crossover(rule.data.num, persen_cross, type = "1Pt")

		### mutation
		rule.data.num <- GA.mutation(rule.data.num, persen_mutant, num.labels, num.inputvar = NULL, type = "MICHIGAN")
		
		stop = stop + 1
		## progress bar
		setTxtProgressBar(progressbar, stop)
	}
	
	close(progressbar)
	## get the best rules
	mod$rule.data.num <- best.rule.data.num
	tmp <- generate.rule(best.rule.data.num, num.labels)
	mod$rule <- tmp$rule	
	mod$range.data.ori <- range.data.ori
	mod$type.defuz <- type.defuz
	mod$type.tnorm <- type.tnorm
	mod$type.snorm <- type.snorm
	mod$type.mf <- type.mf
	mod$num.labels <- num.labels
	mod$type.implication.func <- type.implication.func
    return(mod)
}

#' This is the internal function that implements the Ishibuchi's method based on 
#' genetic cooperative-competitive learning (GFS.GCCL). It is used to handle classification tasks. 
#' Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}.
#' 
#' This method is based on Ishibuchi's method. In this method, a chromosome describes 
#' each linguistic IF-THEN rule using integer as its representation of the antecedent part. 
#' In the consequent part of the fuzzy rules, the heuristic method is applied to 
#' automatically generate the class. The evaluation is calculated for each rule 
#' which means that the performance is not based on the entire rule set. 
#' The outline of the method is as follows. 
#' \itemize{
#' \item Step 1: Generate an initial population of fuzzy IF-THEN rules.
#' \item Step 2: Evaluate each fuzzy IF-THEN rule in the current population.
#' \item Step 3: Generate new fuzzy IF-THEN rules by genetic operators.
#' \item Step 4: Replace a part of the current population with the newly generated rules.
#' \item Step 5: Terminate the algorithm if a stopping condition is satisfied, 
#' otherwise return to Step 2. 
#' } 
#' Additionally, to handle high dimensional data, this method uses "don't care"
#' attributes on the antecedent fuzzy set.
#' 
#' @title GFS.GCCL model building 
#'
#' @param data.train a matrix (\eqn{m \times n}) of normalized data for the training process, 
#'         where \eqn{m} is the number of instances and 
#'         \eqn{n} is the number of variables; the last column is the output variable.
#'         Note the data must be normalized between 0 and 1. 
#' @param popu.size the size of the population which is generated in each generation.
#' @param range.data.input a matrix containing the ranges of the normalized input data.
#' @param num.labels a matrix describing the number of linguistic terms.
#' @param persen_cross a real number between 0 and 1 representing the probability of crossover.
#' @param persen_mutant a real number between 0 and 1 representing the probability of mutation.
#' @param max.gen the maximal number of generations for the genetic algorithm.
#' @param range.data.ori a matrix containing the ranges of the input data.
# @return list of the model data. please have a look at \code{\link{frbs.learn}} for looking complete components.
#' @references
#' H. Ishibuchi, T. Nakashima, and T. Murata, "Performance evaluation of fuzzy classifier systems for 
#' multidimensional pattern classification problems",
#' IEEE trans. on Systems, Man, and Cybernetics - Part B: Sybernetics, vol. 29. no. 5, pp. 601 - 618 (1999).
# @export
GFS.GCCL <- function(data.train, popu.size = 10, range.data.input, num.labels, persen_cross = 0.6, persen_mutant = 0.3, 
                    max.gen = 10, range.data.ori){

	num.labels.inp <- num.labels[1, 1:(ncol(num.labels) - 1), drop = FALSE]
	num.inputvar <- (ncol(data.train) - 1)
	
	## make range of data according to class on data training
	range.data.out <- matrix(c(min(data.train[, ncol(data.train)], na.rm = TRUE) - 0.4999, 
	                           max(data.train[, ncol(data.train)], na.rm = TRUE) + 0.4999), nrow = 2)
	range.data <- cbind(range.data.input, range.data.out)
	type.mf = "TRIANGLE"
	
	## generate initial model using Wang Mendel/Chi technique as heuristics
	if (popu.size < nrow(data.train)){
		num.data <- popu.size
	} else {
		num.data <- nrow(data.train)
	}
	
	mod <- WM(data.train[1 : num.data, ,drop = FALSE], num.labels, type.mf, type.tnorm = "PRODUCT", type.implication.func = "ZADEH", classification = TRUE)
	rule.data.num <- mod$rule.data.num
	num.r <- as.integer(runif(1, min = 1, max = nrow(rule.data.num) + 0.9))
	row.r <- as.integer(runif(num.r, min = 1, max = nrow(rule.data.num) + 0.9))	
	p.dcare	 <- 0.8
	for (j in row.r){
		col.r <- as.integer(runif(1, min = 1, max = ((ncol(rule.data.num) - 1) + 0.9)))
		rand <- runif(1, min = 0, max = 1)
		if (rand <= p.dcare)
			rule.data.num[j, col.r] <- 0
	}

	ant.rules <- rule.data.num[, -ncol(rule.data.num), drop = FALSE]
	num.rules <- nrow(ant.rules)
	if (popu.size > num.rules){
		ant.rules.a <- generate.popu(method.type = "GFS.GCCL", num.var = num.inputvar, 
		                             num.labels = num.labels.inp, popu.size = (popu.size - num.rules))
		ant.rules <- rbind(ant.rules, ant.rules.a)
	}
		
	varinp.mf <- mod$varinp.mf
	
	mod <- GFS.Michigan(data.train, ant.rules, varinp.mf, num.labels, persen_cross, persen_mutant, max.gen, only.GCCL = TRUE)
	
	best.rule <- mod$rule
	best.grade.cert <- mod$grade.cert

	red.mod <- prune.rule(rule.data.num = best.rule, method = "GCCL", indv.fit = best.grade.cert)
		
	mod$rule.data.num <- red.mod$rule
	if (nrow(mod$rule.data.num) < 1){
		stop("frbs cannot generate the fuzzy rules, please set the higher value of popu.size") 
	}
	
	mod$grade.cert <- red.mod$grade.cert
	mod$rule <- generate.rule(red.mod$rule, num.labels)$rule
	mod$rule[, ncol(mod$rule)] <- mod$rule.data.num[, ncol(mod$rule.data.num)]
	mod$num.labels <- num.labels
	mod$range.data.ori <- range.data.ori
	mod$type.mf <- "TRIANGLE"
	mod$type.tnorm <- "PRODUCT"
	mod$type.snorm <- "MAX"
	mod$type.model <- "FRBCS"
	mod$type.implication.func <- "ZADEH"
	return(mod)
}

#' This is the internal function that implements the Ishibuchi's method based on 
#' hybridization of genetic cooperative-competitive learning (GCCL) and Pittsburgh (FH.GBML). It is used to solve classification tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}.
#'
#' This method is based on Ishibuchi's method using the hybridization 
#' of GCCL and the Pittsburgh approach for genetic fuzzy systems. 
#' The algorithm of this method is as follows:
#' \itemize{
#' \item Step 1: Generate population where each individual in the population is a fuzzy rule set. 
#' \item Step 2: Calculate the fitness value of each rule set in the current population.
#' \item Step 3: Generate new rule sets by the selection, crossover, and mutation in 
#'             the same manner as the Pittsburgh-style algorithm. Then, apply iterations of 
#'             the GCCL to each of the generated rule sets with a probability. 
#' \item Step 4: Add the best rule set in the current population to newly generated rule sets 
#'             to form the next population.
#' \item Step 5: Return to Step 2 if the prespecified stopping condition is not satisfied. 
#' }
#'
#' @title FH.GBML model building 
#' 
#' @param data.train a matrix (\eqn{m \times n}) of normalized data for the training process, 
#'         where \eqn{m} is the number of instances and 
#'         \eqn{n} is the number of variables; the last column is the output variable.
#'         Note the data must be normalized between 0 and 1. 
#' @param popu.size the size of the population which is generated in each generation.
#' @param max.num.rule the maximum number of rules.
#' @param persen_cross a real number between 0 and 1 determining the probability of crossover.
#' @param persen_mutant a real number between 0 and 1 determining the probability of mutation.
#' @param max.gen the maximal number of generations for the genetic algorithms.
#' @param num.class a number of the classes.
#' @param range.data.input a matrix containing the ranges of the normalized input data.
#' @param p.dcare a probability of "don't care" attributes occurred.
#' @param p.gccl a probability of GCCL process occurred.
# @return list of the model data. please have a look at \code{\link{frbs.learn}} for looking complete components.
#' @references
#' H. Ishibuchi, T. Yamamoto, and T. Nakashima, "Hybridization of fuzzy GBML approaches for pattern classification 
#' problems," IEEE Trans. on Systems, Man, and Cybernetics-Part B: Cybernetics, 
#' vol. 35, no. 2, pp. 359 - 365 (2005).
# @export
FH.GBML <- function(data.train, popu.size = 10, max.num.rule = 5, persen_cross = 0.6, persen_mutant = 0.3, max.gen = 10, num.class, range.data.input, p.dcare = 0.5, p.gccl = 0.5){
	
	## create progress bar
	progressbar.gbml <- txtProgressBar(min = 0, max = max.gen, style = 3)	
	
	##Inisialitation of Data
	range.data <- range.data.input
	range.data[1, ] <- 0
	range.data[2, ] <- 1	
	type.mf = 1
	
	num.inputvar <- ncol(data.train) - 1 
	num.var <- ncol(data.train)
	num.labels.inp <- matrix(rep(14, num.inputvar), nrow = 1)
	num.labels <- cbind(num.labels.inp, num.class)
	
	##define parameter values of MF on input variables
	varinp.mf <- create.MF(num.inputvar)

	##Inisialitation of population
	popu.rule <- matrix(NA, nrow = 1, ncol = (max.num.rule * num.inputvar))
	popu.comp.rule <- matrix(NA, nrow = 1, ncol = (max.num.rule * (num.inputvar + 1)))
	popu.grade.cert <- matrix(NA, nrow = 1, ncol = max.num.rule)
	popu.fit.rule <- matrix(NA, nrow = 1, ncol = max.num.rule)
	
	for (i in 1 : popu.size){
		num.inputvar <- (ncol(data.train) - 1)
		
		## generate antecent of rules
		rule.data.num <- generate.popu(method.type = "FH.GBML", num.var = num.inputvar, 
		                       num.labels = num.labels.inp, popu.size = max.num.rule)
			
		num.r <- as.integer(runif(1, min = 1, max = nrow(rule.data.num) + 0.9))
		row.r <- as.integer(runif(num.r, min = 1, max = nrow(rule.data.num) + 0.9))		
		for (j in row.r){
			col.r <- as.integer(runif(1, min = 1, max = (ncol(rule.data.num) + 0.9)))
			rand <- runif(1, min = 0, max = 1)
			if (rand <= p.dcare)
				rule.data.num[j, col.r] <- 0
		}
	
		## determine Class-q and grade of certainty
		buffer.rule <- det.fit.rule(data.train, rule.data.num, varinp.mf, num.labels, type.grade.cert = 2)
		
		## get components on buffer.rule
		comp.rule <- buffer.rule$rule
		grade.cert <- buffer.rule$grade.cert
		fit.rule <- buffer.rule$fit.rule
		
		## get grade of certainty and fit.rule as vector by transposing
		temp.grade.cert <- t(grade.cert)
		temp.fit.rule <- t(fit.rule)
		
		## make antecendent of rule and complete rule as vector (to be individual on population)
		ant.rule.tmp <- convert.MatToVec(comp.rule[, -ncol(comp.rule), drop = FALSE])
		comp.rule.tmp <- convert.MatToVec(comp.rule)
	
		## arange all components as population
		popu.rule <- rbind(popu.rule, ant.rule.tmp)
		popu.comp.rule <- rbind(popu.comp.rule, comp.rule.tmp)
		popu.grade.cert <- rbind(popu.grade.cert, temp.grade.cert)
		popu.fit.rule <- rbind(popu.fit.rule, temp.fit.rule)
	}
	
	## delete element NA on initialization step
	popu.rule <- na.omit(popu.rule)
	popu.comp.rule <- na.omit(popu.comp.rule)
	popu.grade.cert <- na.omit(popu.grade.cert)
	popu.fit.rule <- na.omit(popu.fit.rule)

	## Initialization
	iter = 1
	best.err.percent = 100000
	data.test <- data.train[, -ncol(data.train), drop = FALSE]
	
	while (iter <= max.gen){
		##Evaluate fitness function for each individual (a set of fuzzy rules)
		err.percent <- matrix(nrow = nrow(popu.rule), ncol = 1)
		for (i in 1 : nrow(popu.rule)){
			
			## get one row of population
			rule.vector <- matrix(popu.comp.rule[i, ], nrow = 1)
			grade.cert <- matrix(popu.grade.cert[i, ], ncol = 1)
			
			## convert vector of rule into rule.data.num form
			rule <- convert.VecToMat(rule.vector, num.var)			
	
			## Get the best population by calculating fitness
			classifier.rule <- list(rule = rule, grade.cert = grade.cert)
			resClass.train <- GFS.GCCL.test(data.test, classifier.rule, varinp.mf)
			bench <- cbind(data.train[, ncol(data.train), drop = FALSE], resClass.train)
			counter = 0
			for (j in 1 : nrow(bench)){
				if (bench[j, 1] != bench[j, 2]){
					counter = counter + 1
				}					
			}	
			err.percent[i, 1] <- counter / nrow(bench) * 100
		}	
		
		min.err.percent <- min(err.percent)	
		if (min.err.percent < best.err.percent){
			indx.best.indv <- which.min(err.percent)
			best.err.percent <- min.err.percent
			
			## save the best
			best.popu.rule <- popu.rule[indx.best.indv, ,drop = FALSE]
			best.popu.comp.rule <- popu.comp.rule[indx.best.indv, ,drop = FALSE]
			best.rule <- convert.VecToMat(best.popu.comp.rule, num.var)			
			best.grade.cert <- matrix(popu.grade.cert[indx.best.indv, ], ncol = 1)			
		}
		##Crossover
		new.popu.rule <- GA.crossover(popu.rule, persen_cross, type = "1Pt")
		
		##mutation
		new.popu.rule <- GA.mutation(new.popu.rule, persen_mutant, type = "PITTSBURGH", num.labels.inp, num.inputvar)
		
		## GFS-GCCL
		popu.rule <- matrix(NA, nrow = 1, ncol = (max.num.rule * num.inputvar))
		popu.comp.rule <- matrix(NA, nrow = 1, ncol = (max.num.rule * (num.inputvar + 1)))
		popu.grade.cert <- matrix(NA, nrow = 1, ncol = max.num.rule)
		for (i in 1 : nrow(new.popu.rule)){
			r.m <- runif(1, min=0, max =1)
			
			## convert rule in vector into matrix
			ant.rule.vector <- new.popu.rule[i, ,drop = FALSE]
			ant.rules <- convert.VecToMat(ant.rule.vector, num.inputvar)
			if (r.m <= p.gccl){
				### convert rule as vector into ant.rules
				persen_cross.michigan <- persen_cross
				persen_mutant.michigan <- persen_mutant
				max.gen.michigan <- 1
				mod.michigan <- GFS.Michigan(data.train, ant.rules, varinp.mf, num.labels, persen_cross.michigan, persen_mutant.michigan, max.gen.michigan)
						
				## convert matrix into vector
				comp.rule.tmp <- mod.michigan$rule
				
				grade.cert.tmp <- mod.michigan$grade.cert
				ant.rule.tmp <- comp.rule.tmp[, -ncol(comp.rule.tmp), drop = FALSE]
			} else {
				## determine grade.cert, fitness
				buffer.rule <- det.fit.rule(data.train, ant.rules, varinp.mf, num.labels, type.grade.cert = 2)
				comp.rule.tmp <- buffer.rule$rule
				grade.cert.tmp <- buffer.rule$grade.cert
				ant.rule.tmp <- comp.rule.tmp[, -ncol(comp.rule.tmp), drop = FALSE]
			}
			comp.rule.tmp <- convert.MatToVec(comp.rule.tmp)
			grade.cert.tmp <- t(grade.cert.tmp)
			rule.tmp <- convert.MatToVec(ant.rule.tmp)
			
			## make antecendent of rule and complete rule as vector (to be individual on population)
			popu.rule <- rbind(popu.rule, rule.tmp)
			popu.comp.rule <- rbind(popu.comp.rule, comp.rule.tmp)
			popu.grade.cert <- rbind(popu.grade.cert, grade.cert.tmp)
				
		}
		popu.rule <- na.omit(popu.rule)
		popu.comp.rule <- na.omit(popu.comp.rule)
		popu.grade.cert <- na.omit(popu.grade.cert)
		
		iter = iter + 1
		## progress bar
		setTxtProgressBar(progressbar.gbml, iter)
	}
	
	close(progressbar.gbml)
	red.mod <- prune.rule(rule.data.num = best.rule, method = "GCCL", indv.fit = best.grade.cert)		   
	rule.data.num <- red.mod$rule
	if (nrow(rule.data.num) < 1){
		stop("frbs cannot generate the fuzzy rules, please set the higher value of popu.size and max.num.rule") 
	}
	
	rule <- generate.rule(red.mod$rule, num.labels)$rule
	rule[, ncol(rule)] <- rule.data.num[, ncol(rule.data.num)]
	best.grade.cert <- red.mod$grade.cert
	mod <- list(rule = rule, rule.data.num = rule.data.num, grade.cert = best.grade.cert, varinp.mf = varinp.mf, range.data.ori = range.data.input, 
	         num.labels = num.labels, type.mf = "TRIANGLE", type.tnorm = "PRODUCT", type.model = "FRBCS", type.snorm = "MAX", type.implication.func = "ZADEH")
	return(mod)
}

#' This is the internal function that implements 
#' the structural learning algorithm on vague environment (SLAVE). It is used to handle classification tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}.
#'
#' This method is adopted from A. Gonzalez and R. Perez's paper which is applied for 
#' classification problems. SLAVE is based on the iterative rule learning approach 
#' which means that we get only one fuzzy rule in each execution of the genetic algorithm. 
#' In order to eliminate the irrelevant variables in a rule, 
#' SLAVE has a structure composed of two parts: the first part is to represent 
#' the relevance of variables and the second one is to define values of the parameters. 
#' The following steps are conducted in order to obtain fuzzy rules:
#' \itemize{
#' \item Step 1: Use the genetic algorithm process to obtain ONE RULE for the system.
#' \item Step 2: Collect the rule into the final set of rules.
#' \item Step 3: Check and penalize this rule.
#' \item Step 4: If the stopping criteria is satisfied, 
#' the system returns the set of rules as solution. Otherwise, back to Step 1.
#' }
#' This method uses binary codes as representation of the population and 
#' applies the basic genetic operators, i.e., selection, crossover, and mutation on it.
#' And, the best rule is obtained by calculating the degree of consistency and completeness. 
#'
#' @title SLAVE model building 
#' 
#' @param data.train a matrix (\eqn{m \times n}) of normalized data for the training process, 
#'        where \eqn{m} is the number of instances and 
#'        \eqn{n} is the number of variables; The last column is the output variable.
#'        Note the data must be normalized between 0 and 1. 
#' @param persen_cross a real number between 0 and 1 representing the probability of crossover.
#' @param persen_mutant a real number between 0 and 1 representing the probability of mutation.
#' @param max.iter the maximal number of iterations.
#' @param max.gen the maximal number of generations for the genetic algorithm.
#' @param num.labels a number of the linguistic terms.
#' @param range.data.input a matrix containing the ranges of the normalized input data.
#' @param k.lower a lower bound of the noise threshold.
#' @param k.upper an upper bound of the noise threshold.
#' @param epsilon a value between 0 and 1 representing the covering factor.
# @return list of the model data. please have a look at \code{\link{frbs.learn}} for looking complete components.
#' @references
#' A. Gonzalez and R. Perez, "Selection of relevant features in a fuzzy genetic learning algorithm",  
#' IEEE Transactions on Systems, Man, and Cybernetics, Part B: Cybernetics, vol. 31, no. 3, 
#' pp.  417 - 425 (2001).
# @export
SLAVE <- function(data.train, persen_cross = 0.6, persen_mutant = 0.3, max.iter = 10, max.gen = 10, num.labels, range.data.input, k.lower = 0.25, k.upper = 0.75, epsilon = 0.1){
	if (max.iter < 5){
		stop("please set maximum iteration greater than 5")
	}
	
	## initialize 
	type.mf = 1
	best.comp.rule <- matrix(ncol = ncol(data.train))
	num.inputvar <- (ncol(data.train) - 1)
	rule <- matrix(NA, nrow = 1, ncol = (num.inputvar + ncol(data.train)))

	## split data.train based on class
	num.labels.inp <- num.labels[, -ncol(num.labels), drop = FALSE]
	num.class <- num.labels[1, ncol(num.labels)]
		
	## generate varinp.mf
	range.data.input.norm <- range.data.input
	range.data.input.norm[1, ] <- 0
	range.data.input.norm[2, ] <- 1

	## create progress bar
	progressbar <- txtProgressBar(min = 0, max = (num.class * max.iter), style = 3)	
	iter.pb <- 1
			
	## iterate for each class on data.train
	for (i in 1 : num.class){
		max.fit <- matrix(0, nrow = 1, ncol = 1)	
		popu <- generate.popu(method.type = "SLAVE", data.train = data.train, num.var = num.inputvar, 
						num.labels = num.labels.inp, range.data = range.data.input.norm, 
						type.mf = 1, name.class = i)		
		name.class <- i
		data.sample <- data.train[which(data.train[, ncol(data.train)] == name.class), ]
		
		## calculate fitness for each rule
		popu.var <- popu$popu.var
		comp.rule <- popu$popu.val
		varinp.mf <- popu$varinp.mf
		res <- eval.indv.fitness(data.train = data.train, rule = comp.rule, method.type = "SLAVE", data.sample = data.sample, varinp.mf = varinp.mf, 
								 num.labels = num.labels, popu.var = popu.var, k.lower = k.lower, k.upper = k.upper)
		
		### get index which has max value
		best.indv <- which.max(res)		

		## get the best rule(VAL) and VAR
		best.comp.rule <- comp.rule[best.indv, ,drop = FALSE]
		best.var <- popu.var[best.indv, ,drop = FALSE]
		rule.new <- cbind(best.var, best.comp.rule)		
		best.rule <- rule.new
		max.fit <- res[best.indv, 1]

		## check degree considering with rule.new
		VAR.rule.new <- rule.new[, 1 : num.inputvar, drop = FALSE]
		VAL.rule.new <- rule.new[, (num.inputvar + 1) : (ncol(rule.new) - 1), drop = FALSE]

		## get actual rule (with don't care attributes)
		ant.best.rule <- VAR.rule.new * VAL.rule.new
		data.sample <- red.data(data.sample, ant.best.rule, epsilon, type = "SLAVE" , varinp.mf)
		old.num.data <- nrow(data.sample)

		stop.cr <- FALSE
		iter <- 1
		num.sample.fin <- -1	
		while (stop.cr == FALSE){
			## split class
			data.input <- data.sample[, -ncol(data.sample), drop = FALSE]			
			if (old.num.data != num.sample.fin){
				## update the number of data
				old.num.data <- num.sample.fin
				
				popu <- generate.popu(method.type = "SLAVE", data.train = data.train, num.var = num.inputvar,
				                    num.labels = num.labels.inp, range.data = range.data.input.norm,
									type.mf = 1, name.class = i)		
				popu.var <- popu$popu.var
				comp.rule <- popu$popu.val
				ant.rule <- comp.rule[, -ncol(comp.rule), drop = FALSE]
							
				## Genetic Algorithms Procedure
				## 1. Representation Individual in binary form
				res.popu <- convert.IntToBin(ant.rule, num.labels.inp)
				popu.var <- res.popu$popu.var
				popu.val <- res.popu$popu.val
				popu <- cbind(popu.var, popu.val)		
				num.popu <- nrow(popu)
			}
			
			gen <- 1			
			init.indv.fit <- -1			
			while (gen < max.gen && num.popu > 3){
				popu.var <- popu[, 1 : num.inputvar, drop = FALSE]
				popu.val <- popu[, (num.inputvar + 1) : ncol(popu), drop = FALSE]
					
				## crossover
				if (nrow(popu.val) > 3) {
					new.popu.val <- GA.crossover(popu.val, persen_cross, type = "VAL", num.inputvar, num.labels)
					new.popu.var <- GA.crossover(popu.var, persen_cross, type = "2pt", num.inputvar, num.labels)
					new.popu <- cbind(new.popu.var, new.popu.val)
				}
		
				## mutation
				new.popu <- GA.mutation(new.popu, persen_mutant, type = "IRL-BINARY", num.labels, num.inputvar)
		
				# convert binary into int
				new.popu.var <- new.popu[, 1 : num.inputvar ,drop = FALSE]
				new.popu.val.bin <- new.popu[, (num.inputvar + 1) : ncol(new.popu) ,drop = FALSE]				
				new.popu.val.int <- convert.BinToInt(new.popu.val.bin, num.labels.inp)
		
				## convert string to rule form
				new.popu.val.int <- scale.StrToRule(new.popu.val.int, num.labels.inp)
				
				## build complete rules
				new.popu.val <- cbind(new.popu.val.int, name.class)
				
				## calculate individual fitness		
				res <- eval.indv.fitness(data.train = data.train, rule = new.popu.val, method.type = "SLAVE", data.sample = data.sample, varinp.mf = varinp.mf, 
										 num.labels = num.labels, popu.var = new.popu.var, k.lower = k.lower, k.upper = k.upper)
				
				## choose the best indv and add into best population	
				best.indv <- which.max(res)
				new.indv.fit <- max(res)				
				if (new.indv.fit > init.indv.fit) {
					## get the best rule(VAL) and VAR
					best.comp.rule <- new.popu.val[best.indv, ,drop = FALSE]
					best.var <- new.popu.var[best.indv, ,drop = FALSE]
					rule.new <- cbind(best.var, best.comp.rule)
					init.indv.fit <- new.indv.fit
				}

				popu <- cbind(new.popu.var, new.popu.val.bin)
				gen = gen + 1
			}	
			
			## adding the rule
			best.rule <- rbind(best.rule, rule.new)
			
			## check degree considering with rule.new
			VAR.rule.new <- rule.new[, 1 : num.inputvar, drop = FALSE]
			VAL.rule.new <- rule.new[, (num.inputvar + 1) : (ncol(rule.new) - 1), drop = FALSE]
	
			## get actual rule (with don't care attributes)
			ant.best.rule <- VAR.rule.new * VAL.rule.new
			data.sample <- red.data(data.sample, ant.best.rule, epsilon, type = "SLAVE", varinp.mf)
			
			iter = iter + 1	
			iter.pb = iter.pb + 1
			## update the number of data after reducing	
			num.sample.fin <- nrow(data.sample)
			if (num.sample.fin <= 3 || iter >= max.iter){
				stop.cr <- TRUE
				iter.pb <- i * max.iter
			}
			## progress bar
			setTxtProgressBar(progressbar, iter.pb)
		}
		
		## check for fuzzy rule sets for each class
		rule <- rbind(rule, best.rule)
	}	
	close(progressbar)
	
	rule <- na.omit(rule)
	VAR <- rule[, 1 : num.inputvar, drop = FALSE]
	VAL <- rule[, (num.inputvar + 1) : ncol(rule), drop = FALSE]
	rule <- VAL * cbind(VAR, 1)
	rule.data.num <- unique(rule)
	gen.rule.temp <- generate.rule(rule.data.num, num.labels)
	rule <- gen.rule.temp$rule
	colnames(varinp.mf) <- gen.rule.temp$names.varinput
	rule[, ncol(rule)] <- rule.data.num[, ncol(rule.data.num)]
	mod <- list(rule = rule, rule.data.num = rule.data.num, varinp.mf = varinp.mf, num.labels = num.labels, range.data.ori = range.data.input, type.model = "FRBCS",
	             type.mf = "TRIANGLE", type.tnorm = "PRODUCT", type.snorm = "MAX", type.implication.func = "ZADEH")	
	return(mod)
}


#' This is the internal function that implements genetic lateral tuning and rule selection of linguistic fuzzy systems (GFS.LT.RS). 
#' It is used to solve regression tasks. 
#' Users do not need to call it directly, but just use \code{\link{frbs.learn}} and \code{\link{predict}}.
#' 
#' This method was proposed by R. Alcala et al.   
#' GFS.LT.RS implements a evolutionary algorithm for postprocessing in constructing FRBS model. 
#' It uses a new rule representation model based on the linguistic 2-tupples representation that allows
#' the lateral displacement of the labels. This function allows two different tuning which are global and local tuning.
#' 
#' Regarding with evolutionary algorithms, the following are main components:
#' \itemize{
#' \item coding scheme and initial gene pool;
#' \item chromosome evalution;
#' \item crossover operator;
#' \item restarting approach;
#' \item evolutionary model;
#' }
#' In first time, population is constructed by Wang & Mendel's technique. Mean square error (MSE) is used to 
#' calculate chromosome evaluation. This method performs BLX-a in crossover process. 
#' Additionally, rule selection method is performed in order to minimize the number of rules.
#'
#' @title GFS.LT.RS model building 
#'
#' @param data.train a matrix (\eqn{m \times n}) of normalized data for the training process, where \eqn{m} is the number of instances and 
#'        \eqn{n} is the number of variables; the last column is the output variable. Note the data must be normalized between 0 and 1. 
#' @param popu.size the size of the population which is generated in each generation.
#' @param range.data a matrix representing interval of data.
#' @param num.labels a matrix representing the number of linguistic terms in each variables.
#' @param persen_mutant a real number between 0 and 1 determining the probability of mutation.
#' @param max.gen the maximal number of generations of the genetic algorithm.
#' @param mode.tuning a type of tuning which are \code{"LOCAL"} or \code{"GLOBAL"}. 
#' @param type.tnorm a type of t-norm. See \code{\link{inference}}.
#' @param type.snorm a type of s-norm. See \code{\link{inference}}.
#' @param type.implication.func a type of implication function. See \code{\link{WM}}.
#' @param type.defuz a type of defuzzification methods. See \code{\link{defuzzifier}}.
#' @param rule.selection a boolean value representing whether performs rule selection or not.
#' @param range.data.ori a matrix containing the ranges of the original data. 
#' @seealso \code{\link{GFS.LT.RS.test}}, \code{\link{frbs.learn}}, and \code{\link{predict}}
# @return list of the model data. please have a look at \code{\link{frbs.learn}} for looking complete components.
#' @references
#' R. Alcala, J. Alcala-Fdez, and F. Herrera, "A proposal for the genetic lateral tuning of linguistic fuzzy systems and its interaction with
#' rule selection", IEEE Trans. on Fuzzy Systems, Vol. 15, No. 4, pp. 616 - 635 (2007). 
# @export
GFS.LT.RS <- function(data.train, popu.size = 10, range.data, num.labels, persen_mutant, max.gen = 10, mode.tuning = "GLOBAL", 
                      type.tnorm = "MIN", type.snorm = "MAX", type.implication.func = "ZADEH", type.defuz = "WAM", rule.selection = FALSE, range.data.ori) {
	## create progress bar
	progressbar <- txtProgressBar(min = 0, max = max.gen, style = 3)	

	## make testing data from training data
	data.test <- data.train[, 1 : (ncol(data.train) - 1), drop = FALSE]
	real.val <- data.train[, ncol(data.train), drop = FALSE]
	
	## generate initial rules and membership function from  Wang & Mendel's method
	params <- list(mode.tuning = mode.tuning, rule.selection = rule.selection, type.tnorm = type.tnorm, 
	              type.implication.func = type.implication.func)
	res <- generate.popu(method.type = "GFS.LT.RS", data.train = data.train, range.data = range.data,
           	            num.labels = num.labels, popu.size = popu.size, params = params)
	var.mf <- res$var.mf
	popu <- res$popu
	num.rule <- res$num.rule
	rule.data.num <- res$rule.data.num
	init.popu <- popu
	
	## iteration of GA
	gen = 1
	init.err <- 100000
	err <- matrix(nrow = popu.size + 1, ncol = 1)
	L <- 0
		
	while (gen <= max.gen){
		
		## Calculate fitness
		params <- list(mode.tuning = mode.tuning, rule.selection = rule.selection)
		indv.fit <- eval.indv.fitness(data.train = data.train, rule = rule.data.num, method.type = "GFS.LT.RS", varinp.mf = var.mf, 
			num.labels = num.labels, popu.var = popu, params = params)
		
		min.err <- min(indv.fit[, 1])
		if (min.err < init.err) {
			indx.best.indv <- which.min(indv.fit[, 1])
			best.indv <- popu[indx.best.indv, , drop = FALSE]
			init.err <- min.err
			L <- 1
		}
		
		## sort population according to individual fitness
		temp.popu <- cbind(popu, indv.fit)
		temp.popu <- temp.popu[order(temp.popu[, ncol(temp.popu)]), ]
		popu.parent <- temp.popu[1:2, -ncol(temp.popu), drop = FALSE]
	
		## perform crossover
		params <- list(mode.tuning = mode.tuning, rule.selection = rule.selection, num.rule = num.rule)
		new.popu <- GA.crossover(population = popu.parent, persen_cross = 1, type = "2tupple", params = params)
		
		## perform mutation for rule selection only
		popu.rs <- new.popu[, (ncol(new.popu) - num.rule + 1) : ncol(new.popu), drop = FALSE]
		new.popu.rs <- GA.mutation(population = popu.rs, persen_mutant = persen_mutant, type = "BINARY")
		new.popu[, (ncol(new.popu) - num.rule + 1) : ncol(new.popu)] <- new.popu.rs
		
		## calculate individual fitness
		params$mode.tuning <- mode.tuning
		params$rule.selection <- rule.selection
		new.indv.fit <- eval.indv.fitness(data.train = data.train, rule = rule.data.num, method.type = "GFS.LT.RS", varinp.mf = var.mf, 
			num.labels = num.labels, popu.var = new.popu, params = params)
			
		## sort new individual
		temp.new.popu <- cbind(new.popu, new.indv.fit)
		temp.new.popu <- temp.new.popu[order(temp.new.popu[, ncol(temp.new.popu)]), ]
		
		## add new individual into population
		if (temp.popu[(nrow(temp.popu) - 1), ncol(temp.popu)] >= temp.new.popu[2, ncol(temp.new.popu)]){
			temp.popu[(nrow(temp.popu) - 1), ] <- temp.new.popu[1, ]
			temp.popu[nrow(temp.popu), ] <- temp.new.popu[2, ]
		}
		
		## update population
		popu <- temp.popu[, -ncol(temp.popu)]	
		gen = gen + 1
		L = L + 1
		
		if (L > (1/3 * max.gen)){
			popu <- init.popu
		}
		## progress bar
		setTxtProgressBar(progressbar, gen)
	}
	close(progressbar)
	
	if (rule.selection == TRUE){
		rule.data.num <- check.active.rule(rule.data.num, num.rule, best.indv)				
	} else {
		rule.data.num <- rule.data.num				
	}
	params$popu <- best.indv
	params$var.mf <- var.mf
	params$num.rule <- num.rule
	params$mode.tuning <- mode.tuning
	params$rule.selection <- rule.selection
	res.tune <- tune.MF(method.type = "GFS.LT.RS", params = params)
	var.mf <- res.tune$var.mf
	var.mf.tune <- res.tune$var.mf.tune
	rule <- generate.rule(rule.data.num, num.labels)$rule
	varinp.mf <- var.mf[, 1 : (ncol(var.mf) - num.labels[1, ncol(num.labels)]), drop = FALSE]
	varout.mf <- var.mf[, (ncol(varinp.mf) + 1) : ncol(var.mf), drop = FALSE]
	mod <- list(rule = rule, rule.data.num = rule.data.num, var.mf = var.mf, varinp.mf = varinp.mf, varout.mf = varout.mf, var.mf.tune = var.mf.tune, num.labels = num.labels, 
	                mode.tuning = mode.tuning, rule.selection = rule.selection, range.data.ori = range.data.ori, 
					type.defuz = type.defuz, type.tnorm = type.tnorm, type.snorm = type.snorm, type.mf = "TRIANGLE",
					type.implication.func = type.implication.func, type.model = "2TUPPLE")
	return(mod)
}

