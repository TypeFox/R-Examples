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
#
# This function applies genetic algorithm based on genetic cooperative-competitive learning approach. 
#
# @title Genetic algorithm based on genetic cooperative-competitive learning approach 
#
# @param data.train a matrix (\eqn{m \times n}) of normalized data for the training process, where \eqn{m} is the number of instances and 
# \eqn{n} is the number of variables; the last column is the output variable. Note the data must be normalized between 0 and 1. 
# @param ant.rules a matrix of antecedent parts of rules.
# @param varinp.mf a matrix which is parameter values of membership function.
# @param num.labels a matrix defining the number of linguistic terms.
# @param persen_cross a real number between 0 and 1 representing the probability of crossover.
# @param persen_mutant  a real number between 0 and 1 representing the probability of mutation.
# @param max.gen the maximal number of generation on genetic algorithms.
# @export
GFS.Michigan <- function(data.train, ant.rules, varinp.mf, num.labels, persen_cross, 
				persen_mutant, max.gen, only.GCCL = FALSE){
	## create progress bar
	if (only.GCCL){
		progressbar <- txtProgressBar(min = 0, max = max.gen, style = 3)	
	}
	
	## Initialize population
	mod <- NULL
	err.percent <- matrix()
	num.inputvar <- (ncol(data.train) - 1)
	
	##### Determine Class C and Certanty Factor grade.cert
	buffer.rule <- det.fit.rule(data.train, ant.rules, varinp.mf, num.labels, type.grade.cert = 1)
	
	#### Process of GA
	iter = 1
	best.err.percent = 100000
	data.test <- data.train[, -ncol(data.train), drop = FALSE]
	
	while (iter <= max.gen){
		
		old.buffer.rule <- buffer.rule
		rule <- old.buffer.rule$rule
		grade.cert <- old.buffer.rule$grade.cert
		fit.rule <- old.buffer.rule$fit.rule
	
		## Get the best population by calculating fitness
		classifier.rule <- old.buffer.rule
		resClass.train <- GFS.GCCL.test(data.test, classifier.rule, varinp.mf)
		err.percent = 100*sum(data.train[, ncol(data.train), drop = FALSE] != resClass.train)/nrow(resClass.train)
		
		if (err.percent < best.err.percent){
			best.classifier.rule <- classifier.rule
			best.err.percent <- err.percent
		}
		
		## crossover	
		rule.ant <- rule[, -ncol(rule), drop = FALSE]
		new.popu <- GA.crossover(rule.ant, persen_cross, type = "1Pt")
	
		## mutation
		new.popu <- GA.mutation(new.popu, persen_mutant, num.labels, num.inputvar = NULL, type = "MICHIGAN")	
		
		## determine class of new population
		temp.buffer.rule <- det.fit.rule(data.train, new.popu, varinp.mf, num.labels, type.grade.cert = 1)

		#### Update rule	
		temp <- update.rule(old.buffer.rule, temp.buffer.rule)

		buffer.rule <- list(rule = temp[, 1 : (num.inputvar + 1), drop = FALSE], 
		                    grade.cert = temp[, (ncol(temp) - 1), drop = FALSE],
							fit.rule = temp[, ncol(temp), drop = FALSE])

		iter = iter + 1
		if (only.GCCL){
			## progress bar
			setTxtProgressBar(progressbar, iter)
		}
	}	
	
	if (only.GCCL){
		close(progressbar)
	}
	best.rule <- best.classifier.rule$rule
	best.grade.cert <- best.classifier.rule$grade.cert
	
	mod <- list(rule = best.rule, grade.cert = best.grade.cert, varinp.mf = varinp.mf)
	
	return(mod)
}

# This function is to add new rule into a set of rules
#
# @param old.rule a matrix of a set of rules
# @param new.rule a matrix of a new rule
update.rule <- function(old.rule, new.rule){
	old.data <- cbind(old.rule$rule, old.rule$grade.cert, old.rule$fit.rule)
	new.data <- cbind(new.rule$rule, new.rule$grade.cert, new.rule$fit.rule)
	res <- old.data
	res[which.min(old.data[, ncol(old.data)]), ] <- new.data[which.max(new.data[, ncol(new.data)]), ]
	return(res)
}

# This function calculates a fitness value of each chromosome.
#
# @param data.train a matrix of training data
# @param ant.rules a matrix of the antecedent part of rules
# @param num.labels a matrix of the number of linguistic terms
# @param type.grade.cert a type of a method used to calculate certainty factor on consequent part of rules.
det.fit.rule <- function(data.train, ant.rules, varinp.mf, num.labels, type.grade.cert){
	
	res.rule <- det.grade.cert(data.train, ant.rules, varinp.mf, num.labels, type.grade.cert)
	
	##### Classify all the given training pattern and calculate the fitness value on each rule
	rule <- res.rule$comp.rule.data.num
	grade.cert <- res.rule$grade.cert
	fit.rule = eval.indv.fitness(data.train = data.train, rule = rule, method.type = "GFS.GCCL", grade.cert = grade.cert, varinp.mf = varinp.mf)	
	buffer.rule <- list(rule = rule, grade.cert = grade.cert, fit.rule = fit.rule)
	
	return(buffer.rule)
}

# This function performs mutation operator on genetic algorithm
#
# @param population a matrix of population
# @param persen_mutant a value of mutation probability will be occurred.
# @param num.labels a matrix of the number of linguistic terms.
# @param num.inputvar a number of input variables
# @param iter a interation counter
# @param max.iter a maximal number of iteration
# @param type a type of genetic algorithm approach
GA.mutation <- function(population, persen_mutant, type = "MICHIGAN", num.labels = NULL, 
                        num.inputvar = NULL, iter = NULL, max.iter = NULL){
	if (type == "MICHIGAN"){
		new.popu <- population
		
		## iterate along number of population
		for (i in 1 : nrow(population)){
			
			## generate a real number randomly
			r = runif(1, min = 0, max = 1)
			if (r <= persen_mutant){
				
				## generate a location of mutation
				point.mutation <- as.integer(runif(1, min = 1, max = ncol(population) + 0.99))
				
				## min & max values
				min.labels.cons <- (num.labels[1,1] * (point.mutation - 1) + 1)
				max.labels.cons <- min.labels.cons  + num.labels[1,1] - 1 
				
				## mutation change one up and down of label
				run.updown <- sample(c(-1,1), 1, replace = TRUE)
				new.label = population[i, point.mutation] + run.updown
								
				if (new.label >= (max.labels.cons + 1)){
					new.popu[i, point.mutation] <- 0
				}
				else if (new.label <= (min.labels.cons - 1)){
					new.popu[i, point.mutation] <- 0
				} else {
					new.popu[i, point.mutation] <- new.label
				}
			}
		}
	}
	else if (type == "PITTSBURGH"){
		new.popu <- population
		if (is.null(num.inputvar)) {
			stop("please insert the number of input variables")
		}
		for (i in 1 : nrow(population)){
			r = runif(1, min = 0, max = 1)		
			if (r <= persen_mutant){
				point.mutation <- as.integer(runif(1, min = 1, max = ncol(population) + 0.99))
				
				indx <- point.mutation %% num.inputvar
				if (indx == 0){
					point.term <- num.inputvar
				} else {
					point.term <- indx
				}
				
				min.labels.cons <- (num.labels[1,1] * (point.term - 1) + 1)
				max.labels.cons <- min.labels.cons  + num.labels[1,1] - 1 
				
				run.updown <- sample(c(-1,1), 1, replace = TRUE)
				new.label = population[i, point.mutation] + run.updown
								
				if (new.label >= (max.labels.cons + 1)){
					new.popu[i, point.mutation] <- 0
				}
				else if (new.label <= (min.labels.cons - 1)){
					new.popu[i, point.mutation] <- 0
				} else {
					new.popu[i, point.mutation] <- new.label
				}
			}
		}
	}	
	else if (type == "IRL-BINARY"){
		new.popu <- population
		new.popu.var <- new.popu[, 1 : num.inputvar, drop = FALSE]
		new.popu.val <- new.popu[, (num.inputvar + 1) : ncol(new.popu), drop = FALSE]
		for (i in 1 : nrow(population)){
			r = runif(1, min = 0, max = 1)		
			if (r <= persen_mutant){
				## mutation on var
				temp <- matrix(0, nrow = 1, ncol = num.inputvar)
				loc.var.r <- as.integer(runif(1, min = 1, max = num.inputvar + 0.99))
				temp[1, loc.var.r] = 1
				new.popu.var[i, ] <- matrix(as.integer(xor(new.popu.var[i, ], temp)), nrow = 1)
				
				## mutation on val
				# get variable which is changed
				temp.r <- as.integer(runif(1, min = 0, max = (num.inputvar - 1) + 0.99))
				# get index exactly 
				start <- temp.r * num.labels[1,1] + 1
				end <- start + num.labels[1,1] - 1
				# location of mutation
				loc.val.r <- as.integer(runif(1, min = start, max = end + 0.99))
				
				# perform mutation
				new.popu.val[i, start : end] <- 0
				new.popu.val[i, loc.val.r] <- 1
			}
		}
		new.popu <- cbind(new.popu.var, new.popu.val)
	}
	else if (type == "MOGUL"){
		##mutation		
		new.popu <- population
		for (i in 1 : nrow(new.popu)){
			loc.mutation <- runif(ncol(new.popu), min = 0, max = 1)		
			for (j in 1 : ncol(new.popu)){
				if (loc.mutation[j] < persen_mutant){
					cond <- as.integer(runif(1, min = 0, max = 1.99))
					r <- runif(1, min = 0, max = 1)
					b = 1
					if (cond == 0){
						y <- (1 - new.popu[i, j])
						delta <- y * (1 -  r ^ ((1 - iter/max.iter) ^ b))
						new.popu[i, j] <- (new.popu[i, j] + delta)					
					} else {
						y <- new.popu[i, j]
					 	delta <- y * (1 -  r ^ ((1 - iter/max.iter) ^ b))
						new.popu[i, j] <- new.popu[i, j] - delta
					}
				}
				if (j %% 3 == 1){
					new.popu[i, j] <- min(new.popu[i, j], new.popu[i, j + 1])
				}
				else if (j %% 3 == 2){
					new.popu[i, j] <- max(new.popu[i, j - 1], new.popu[i, j])
					new.popu[i, j] <- min(new.popu[i, j], new.popu[i, j + 1])
				} else{
					new.popu[i, j] <- max(new.popu[i, j], new.popu[i, j - 1])
				}
			}
		}
	}
	else if (type == "BINARY"){
		new.popu <- population
		for (i in 1 : nrow(population)){
			r = runif(1, min = 0, max = 1)
			if (r <= persen_mutant) {
				op.rand <- matrix(round(runif(ncol(population), min = 0, max = 1)), nrow = 1)
				new.popu[i, ] <- matrix(as.integer(xor(op.rand, population[i, , drop = FALSE])), nrow = 1)
			}
		}
	}
	return(new.popu)
}

# This function performs crossover operator on genetic algorithms
#
# @param population a matrix of population
# @param persen_cross a value of the crossover probability
# @param type a type of crossover methods
# @param num.inputvar a number of input variables
# @param num.labels a matrix of the number of linguistic terms
# @param params a list of other parameters
GA.crossover <- function(population, persen_cross, type = "1Pt", num.inputvar = NULL, 
                num.labels = NULL, params = list()){
	if (type == "1Pt"){
		median.popu <- floor(nrow(population)/2)	
		num.parent <- runif(1, min = 1, max = median.popu)
		indx.parent <- sample(seq(1, median.popu), num.parent, replace = FALSE)
		new.popu <- population
		
		### crossover
		for (i in 1 : length(indx.parent)){			
			r = runif(1, min = 0, max = 1)
			if (r <= persen_cross){
				point.cross <- as.integer(runif(1, min = 2, max = ncol(population) + 0.99))
				new.popu[indx.parent[i], point.cross : ncol(new.popu)] <- population[indx.parent[i] + median.popu, point.cross : ncol(population)]
				new.popu[indx.parent[i] + median.popu, point.cross : ncol(new.popu)] <- population[indx.parent[i], point.cross : ncol(population)]
			}
		}
	}
	else if (type == "VAL"){
		new.popu <- population
		median.popu <- floor(nrow(population)/2)	
		num.parent <- runif(1, min = 1, max = median.popu)
		indx.parent <- sample(seq(1, median.popu), num.parent, replace = FALSE)
		
		for (i in 1 : length(indx.parent)){
			r = runif(1, min = 0, max = 1)
			if (r <- persen_cross){
				## we use this technique in order to move variable code instead of binary itself
				temp.r.1 <- as.integer(runif(1, min = 0, max = (num.inputvar - 1) + 0.99))
				temp.r.2 <- as.integer(runif(1, min = (temp.r.1 + 1), max = num.inputvar + 0.99))
				# get index exactly 
				point.cross.1 <- temp.r.1 * num.labels[1,1] + 1
				point.cross.2 <- temp.r.2 * num.labels[1,1]
				new.popu[indx.parent[i], point.cross.1 : point.cross.2] <- population[indx.parent[i] + median.popu, point.cross.1 : point.cross.2]
				new.popu[indx.parent[i] + median.popu, point.cross.1 : point.cross.2] <- population[indx.parent[i], point.cross.1 : point.cross.2]
			}
		}
	}
	else if (type == "2pt"){
		median.popu <- floor(nrow(population)/2)	
		num.parent <- runif(1, min = 1, max = median.popu)
		indx.parent <- sample(seq(1, median.popu), num.parent, replace = FALSE)
		new.popu <- population
		
		### crossover
		for (i in 1 : length(indx.parent)){			
			r = runif(1, min = 0, max = 1)
			if (r <= persen_cross){
				point.cross.1 <- as.integer(runif(1, min = 1, max = (ncol(population) - 1) + 0.99))
				point.cross.2 <- as.integer(runif(1, min = point.cross.1, max = ncol(population) + 0.99))
				new.popu[indx.parent[i], point.cross.1 : point.cross.2] <- population[indx.parent[i] + median.popu, point.cross.1 : point.cross.2]
				new.popu[indx.parent[i] + median.popu, point.cross.1 : point.cross.2] <- population[indx.parent[i], point.cross.1 : point.cross.2]
			}
		}
	}
	else if (type == "MOGUL"){
		### crossover 
		## determine the number of parents
		new.popu <- population
		median.popu <- floor(nrow(population)/2)	
		num.parent <- runif(1, min = 1, max = median.popu)
		indx.parent <- sample(seq(1, median.popu), num.parent, replace = FALSE)
	
		for (i in 1 : length(indx.parent)){
			r = runif(1, min = 0, max = 1)
			if (r <- persen_cross){
				## we use this technique in order to move variable code instead of binary itself
				temp.r.1 <- as.integer(runif(1, min = 0, max = (num.inputvar - 1) + 0.99))
				temp.r.2 <- as.integer(runif(1, min = (temp.r.1 + 1), max = num.inputvar + 0.99))
				# get index exactly 
				point.cross.1 <- temp.r.1 * 3 + 1
				point.cross.2 <- temp.r.2 * 3
				new.popu[indx.parent[i], point.cross.1 : point.cross.2] <- population[indx.parent[i] + median.popu, point.cross.1 : point.cross.2]
				new.popu[indx.parent[i] + median.popu, point.cross.1 : point.cross.2] <- population[indx.parent[i], point.cross.1 : point.cross.2]
			}
		}	
	}
	else if (type == "2tupple"){
		
		mode.tuning <- params$mode.tuning
		rule.selection <- params$rule.selection
		num.rule <- params$num.rule
		popu.parent <- population
		
		if (rule.selection == TRUE){
			popu.lateral <- popu.parent[, 1 : (ncol(popu.parent) - num.rule), drop = FALSE]
			popu.rule <- popu.parent[, (ncol(popu.lateral) + 1) : ncol(popu.parent), drop = FALSE]
			new.popu.rule <- GA.crossover(population = popu.rule, persen_cross, type = "2pt")
			new.popu.lateral <- BLX.alpha(population = popu.lateral, persen_cross = persen_cross, alpha = 0.3)
			new.popu <- cbind(new.popu.lateral, new.popu.rule)	
		} else {
			popu.lateral <- popu.parent
			new.popu.lateral <- BLX.alpha(population = popu.lateral, persen_cross = persen_cross, alpha = 0.3)
			new.popu <- new.popu.lateral
		}
		
	}
	return(new.popu)
}

# This function is used to perform crossover based on BLX.alpha technique.
#
# @param population a matrix of population
# @param persen_cross a percentage of crossover occurred
# @param alpha a multiplication value
BLX.alpha <- function(population, persen_cross, alpha){
	new.popu <- population
	Cmax <- matrix(apply(new.popu, 2, max), nrow = 1)
	Cmin <- matrix(apply(new.popu, 2, min), nrow = 1)
	I <- matrix(Cmax - Cmin, nrow = 1)

	for (j in 1 : ncol(new.popu)){				
		new.popu[1, j] <- runif(1, min = Cmin[1, j] - I[1, j] * alpha, max = Cmax[1, j] + I[1, j] * alpha)
		new.popu[2, j] <- runif(1, min = Cmin[1, j] - I[1, j] * alpha, max = Cmax[1, j] + I[1, j] * alpha)		
	}

	return(new.popu)
}

# This function is to determine class as consequent part value
#
# @param data.test a matrix of testing data
# @param rule a matrix of a set of rules
# @param grade.cert a matrix of certainty factor
# @param varinp.mf a matrix of parameter values of membership function
# @param type a type of consequent part whether it is weighted or not.
det.class <- function(data.test, rule, grade.cert = NULL, varinp.mf, type = "WEIGHTED"){
	data.input <- data.test
	rule.ant <- matrix(rule[, -ncol(rule)], nrow=nrow(rule))
	degree <- matrix(nrow=nrow(rule.ant), ncol=ncol(rule.ant))
	
	f.deg <- function(j, k, rule.ant){	
		if (rule.ant[j, k] == 0){
			degree[j, k] <<- 1
		} else {
			#check for triangular
			#get term on rule
			term.act <- rule.ant[j,k]
			if (varinp.mf[1, term.act] == 1){
				if (data.input[1,k] <= varinp.mf[2, term.act]){
					degree[j, k] <<- 0
				}
				else if (data.input[1,k] <= varinp.mf[3, term.act]) {
					degree[j, k] <<- (data.input[1,k] - varinp.mf[2, term.act]) / (varinp.mf[3, term.act] - varinp.mf[2, term.act])
				}
				else if (data.input[1,k] <= varinp.mf[4, term.act]) {
					degree[j, k] <<- (varinp.mf[4, term.act] - data.input[1,k]) / (varinp.mf[4, term.act] - varinp.mf[3, term.act])
				} else {
					degree[j, k] <<- 0
				}
			}
		}		
	}	
	vec.f.deg <- Vectorize(f.deg, vectorize.args=list("j","k"))
	outer(1 : nrow(rule.ant), 1 : ncol(rule.ant), vec.f.deg, rule.ant)
	
	prod.degree = matrix(apply(degree, 1, prod), ncol = 1) 
	
	if (type == "WEIGHTED"){
		temp = prod.degree * grade.cert
	}
	else if (type == "NON-WEIGHTED"){
		temp = prod.degree
	}
	
	ind.rule.act = which.max(temp)
	
	res.class = rule[ind.rule.act, ncol(rule)]
	
	res <- list(res.class = res.class, ind.rule.act = ind.rule.act)
	return(res)
}

# This function determines class and calculates certainty factor on consequent part
#
# @param data.train a matrix of training data
# @param ant.rules a matrix of the antecedent part of rules
# @param num.labels a matrix of the number of linguistic terms
# @type.grade.cert a type of certainty factor
det.grade.cert <- function(data.train, ant.rules, varinp.mf, num.labels, type.grade.cert = 1){
	data.input <- data.train[, 1 : (ncol(data.train) - 1), drop = FALSE]
	data.class <- data.train[, ncol(data.train), drop = FALSE]

	### Step 1: Calculate compatibility grade
	degree <- matrix(nrow=nrow(data.input), ncol=ncol(data.input))
	miu.rule <- matrix(nrow=nrow(data.input), ncol=1)
	grade.cert <- matrix(nrow=nrow(ant.rules), ncol =1)
	
	comp.rule.data.num <- cbind(ant.rules, 0)
	
	if (type.grade.cert == 1){
		##iterate for each rule
		for (i in 1 : nrow(ant.rules)){
			##apply all data.input for each rule
			ant.rules.i <- ant.rules[i, ,drop = FALSE]
			
			## get product of degree on each attribute
			miu.rule <- calc.degree.ant(data.train, ant.rules.i, varinp.mf)

			## combine miu.rule and class of data
			temp.consq <- cbind(miu.rule, data.class)	
		
			num.class <- num.labels[1, ncol(num.labels)]
			
			## calculate aggregate based on class
			beta.class = as.matrix(aggregate(x=temp.consq[,1], by=list(temp.consq[,2]), FUN="sum"))
		
			## determine class
			ind.max.class <- which.max(beta.class[, 2])

			## check duplication on max of value on certain class
			check.num.max <- which(beta.class[, 2] == beta.class[ind.max.class, 2])
		
			## check condition to determine class
			if (max(beta.class[, 2]) != 0 && length(check.num.max) == 1){
				comp.rule.data.num[i, ncol(comp.rule.data.num)] = ind.max.class
			}

			if (comp.rule.data.num[i, ncol(comp.rule.data.num)] == 0 || sum(beta.class[, 2]) == 0){
				grade.cert[i, 1] = 0
			} else{
				temp = max(beta.class[, 2])
				sum.beta.class = sum(beta.class[, 2])
				grade.cert[i, 1] = (temp - ((sum.beta.class - temp)/(nrow(beta.class) - 1)))/sum.beta.class
			}
		}
	}
	else if (type.grade.cert == 2){
		for (i in 1 : nrow(ant.rules)){
			##apply all data.input for each rule
			ant.rules.i <- ant.rules[i, ,drop = FALSE]
			
			## get product of degree on each attribute
			miu.rule <- calc.degree.ant(data.train, ant.rules.i, varinp.mf)
			
			## combine miu.rule and class of data
			temp.consq <- cbind(miu.rule, data.class)	
			num.class <- num.labels[1, ncol(num.labels)]
			
			## calculate aggregate based on class
			beta.class = as.matrix(aggregate(x=temp.consq[,1], by=list(temp.consq[,2]), FUN="sum"))
			
			if (sum(beta.class[, 2]) != 0){
				Pr.class = matrix(beta.class[, 2] / sum(beta.class[, 2]), ncol = 1)			
				indx.max.class = which.max(Pr.class)
				comp.rule.data.num[i, ncol(comp.rule.data.num)] = indx.max.class
				grade.cert[i, 1] = Pr.class[indx.max.class] - (sum(Pr.class) - Pr.class[indx.max.class])
				if (grade.cert[i, 1] < 0){
					comp.rule.data.num[i, ncol(comp.rule.data.num)] = 0
					grade.cert[i, 1] = 0
				}
			}
			else {
				comp.rule.data.num[i, ncol(comp.rule.data.num)] = 0
				grade.cert[i, 1] = 0
			}			
		}
	}
				
	## delete NA in rule
	comp.rule.data.num = na.omit(comp.rule.data.num)
	grade.cert = na.omit(grade.cert)
	
	res = list(comp.rule.data.num = comp.rule.data.num, grade.cert = grade.cert)
	return (res)
}

# This function is to calculate a degree of antecedent part on rules
#
# @param data.train a matrix of training data
# @param ant.rules.i antecedent part of one rule
# @param varinp.mf a matrix of parameter values of membership function
calc.degree.ant <- function(data.train, ant.rules.i, varinp.mf){
	data.input <- data.train[, 1 : (ncol(data.train) - 1), drop = FALSE]
	data.class <- data.train[, ncol(data.train), drop = FALSE]

	### Step 1: Calculate compatibility grade
	degree <- matrix(nrow=nrow(data.input), ncol=ncol(data.input))
	miu.rule <- matrix(nrow=nrow(data.input), ncol=1)
	
	f.deg <- function(j, k, data.input){	
		if (ant.rules.i[1, k] == 0){
			degree[j, k] <<- 1
		} else {
			#check for triangular
			#get term on rule
			term.act <- ant.rules.i[1,k]
			if (varinp.mf[1, term.act] == 1){
				if (data.input[j,k] < varinp.mf[2, term.act]){
					degree[j, k] <<- 0						
				}
				else if (data.input[j,k] < varinp.mf[3, term.act]) {
					degree[j, k] <<- (data.input[j,k] - varinp.mf[2, term.act]) / (varinp.mf[3, term.act] - varinp.mf[2, term.act])
				}
				else if (data.input[j,k] < varinp.mf[4, term.act]) {
					degree[j, k] <<- (varinp.mf[4, term.act] - data.input[j,k]) / (varinp.mf[4, term.act] - varinp.mf[3, term.act])
				} else {
					degree[j, k] <<- 0
				}
			}				
		}	
	}		
	vec.f.deg <- Vectorize(f.deg, vectorize.args=list("j","k"))
	outer(1 : nrow(data.input), 1 : ncol(data.input), vec.f.deg, data.input)

	## get product of degree on each attribute
	miu.rule = matrix(apply(degree, 1, prod), ncol = 1) 

	return(miu.rule)
}

# This function is to reduce training data as input data
# @param data.train a matrix of training data
# @param ant.rules.i a matrix of single rule
# @param epsilon a parameter representing a threshold value
# @param type a type of method
# @param varinp.mf a matrix of membership function
red.data <- function(data.train, ant.rules.i, epsilon, type = "SLAVE", varinp.mf = NULL){
	if (type == "SLAVE"){
		degree.ant <- calc.degree.ant(data.train, ant.rules.i, varinp.mf)
		data.train[which(degree.ant[, 1] >= epsilon), ] <- NA
	}
	else if (type == "MOGUL"){
		for (j in 1 : nrow(data.train)){
			## check completeness and covering ==> ch.cover
			comp.degree <- ch.cover(data.train[j, ,drop = FALSE], ant.rules.i)	
			if (comp.degree >= epsilon){
				data.train[j, ] <- NA
			}
		}		
	}
	data.train <- na.omit(data.train)
	return(data.train)
}

# This function is to determine class of new instances
#
# @param data.test a matrix of testing data
# @param classifier.rule a complete rule containing rules and certainty factor (grade.cert)
# @param varinp.mf a matrix of parameter values of the membership function.
GFS.GCCL.test <- function(data.test, classifier.rule, varinp.mf){
	rule <- classifier.rule$rule
	grade.cert <- classifier.rule$grade.cert
	res <- matrix(nrow=nrow(data.test), ncol = 1)
	
	for (i in 1 : nrow(data.test)){
		temp <- det.class(matrix(data.test[i, ], nrow = 1), rule, grade.cert, varinp.mf)
		res[i, 1] <- temp$res.class
	}

	return(res)
}

# This function creates a matrix of parameter values of the membership function.
# 
# @param num.inputvar a number of input variables
create.MF <- function(num.inputvar){
	##define parameter values of MF on input variables
	## partition.MF function is on GFS.Thrift
	var.mf.i <- matrix(nrow = 5, ncol = 14)
	
	num.labels.2 <- matrix(c(2), nrow = 1)
	num.labels.3 <- matrix(c(3), nrow = 1)
	num.labels.4 <- matrix(c(4), nrow = 1)
	num.labels.5 <- matrix(c(5), nrow = 1)
	var.mf.2 <- partition.MF(matrix(c(0, 1), nrow = 2), num.labels.2, type.mf = "TRIANGLE")
	var.mf.3 <- partition.MF(matrix(c(0, 1), nrow = 2), num.labels.3, type.mf = "TRIANGLE")
	var.mf.4 <- partition.MF(matrix(c(0, 1), nrow = 2), num.labels.4, type.mf = "TRIANGLE")
	var.mf.5 <- partition.MF(matrix(c(0, 1), nrow = 2), num.labels.5, type.mf = "TRIANGLE")
	var.mf.i <- cbind(var.mf.2, var.mf.3, var.mf.4, var.mf.5)

	var.mf <- matrix(rep(var.mf.i, num.inputvar), nrow = 5)
	names.one.var <- c("s.2", "l.2", "s.3", "m.3", "l.3", "s.4", "ms.4", "ml.4", "l.4", "s.5", "ms.5", "m.5", "ml.5", "l.5")
	names.temp <- list()
	for (i in 1 : num.inputvar){
		n.var <- paste("v", i, sep = ".")
		temp <- paste(n.var, names.one.var, sep = "_")
		names.temp <- append(names.temp, temp)
	}
	names.var <- as.character(names.temp)
	colnames(var.mf) <- names.var
	return(var.mf)
}

# This function is to convert rules using matrix form into vector form
# 
# @param data.mat a matrix of rules
convert.MatToVec <- function(data.mat){
	## define varibles
	tmp <- data.mat[1, ,drop = FALSE]
		
	## convert matrix into vector
	for (i in 2 : nrow(data.mat)){
		tmp <- cbind(tmp, matrix(data.mat[i, ], nrow = 1))
	}	
	
	return(tmp)
}

# This function is to convert rules using vector form into matrix form
#
# @param data.vec a vector of rules
# @param num.var a number of variables
convert.VecToMat <- function(data.vec, num.var){
	data.mat <- matrix(NA, nrow = 1, ncol = num.var)
	max.num.rule <- ncol(data.vec)/num.var
	for (j in 0 : (max.num.rule - 1)){
		start.rule <- j * num.var + 1
		end.rule <- start.rule + num.var - 1
		temp <- data.vec[1, start.rule : end.rule, drop = FALSE]
		data.mat <- rbind(data.mat, temp)
	}
	data.mat <- na.omit(data.mat)
	
	return(data.mat)
}

# This function is used to calculate error
#
# @param mod a frbs object
# @param data.train a matrix of training data
# @param method.type a type of method
calc.MSE <- function(mod, data.train, method.type){
	data.test <- data.train[, -ncol(data.train), drop = FALSE]
	
	if (method.type == "SLAVE"){
		res.test <- SLAVE.test(mod, data.test)
	
		y.pred <- res.test
		y.real <- data.train[, ncol(data.train) ,drop = FALSE]
		bench <- cbind(y.pred, y.real)
		counter <- 0
		for (i in 1 : nrow(bench)){
			if (bench[i, 1] != bench[i, 2]){
				counter = counter + 1
			}
		}
		err <- counter / nrow(bench) * 100
	}
	else if (method.type == "GFS.FR.MOGUL"){
		rule.gen <- mod$rule
		real.val <- data.train[, ncol(data.train), drop = FALSE]
		y.pred <- GFS.FR.MOGUL.test(mod, data.test)
		y.real <- real.val

		#### Measure error for classification
		residuals <- (y.real - y.pred)
		RMSE <- sqrt(mean(residuals^2))
		err <- RMSE
	}
	
	else if (method.type == "GFS.LT.RS"){
		y.pred <-	GFS.LT.RS.test(mod, data.test)
		y.real <- data.train[, ncol(data.train), drop = FALSE]
		residuals <- (y.real - y.pred)
		RMSE <- sqrt(mean(residuals^2))
		err <- RMSE
	}
	return (err)
}

# This function is used to check and define the active rule
# 
# @param rule.data.num rules in matrix format
# @param num.rule number of rules
# @param popu a matrix representing a population
check.active.rule <- function(rule.data.num, num.rule, popu){	
	#NA.rule <- popu.rule

	NA.rule <- popu[, (ncol(popu) - num.rule + 1) : ncol(popu), drop = FALSE]
	
	NA.rule[which(NA.rule == 0)] <- NA
	NA.rule.i <- t(NA.rule[1, ,drop = FALSE])
	rule.data.num.rs <- cbind(rule.data.num, NA.rule.i)
	rule.data.num.rs <- na.omit(rule.data.num.rs)
	rule.data.num.rs <- rule.data.num.rs[, -ncol(rule.data.num.rs), drop = FALSE]
	if (nrow(rule.data.num.rs) < 1){
		rule.data.num.rs <- rule.data.num
	}

	return(rule.data.num.rs)
}


# This function is used to convert rule in binary form into integer form
#
# @param data.bin a matrix of rule in binary form
# @param num.labels.inp a matrix of the number of linguistic terms of input variables
convert.BinToInt <- function(data.bin, num.labels.inp){
	data.int <- matrix(nrow = nrow(data.bin), ncol = ncol(num.labels.inp))
	for (i in 1 : nrow(data.bin)){
		for (j in 1 : ncol(num.labels.inp)){
			start = (j - 1) * num.labels.inp[1, j] + 1
			end = start + num.labels.inp[1, j] - 1
			
			temp <- data.bin[i, start : end, drop = FALSE]
			loc.term <- which.max(temp)
			data.int[i, j] <- loc.term
		}
	}
	
	return(data.int)
}

# This function converts rule in integer form into binary form
#
# @param data.mat a matrix of rules in integer form
# @param num.labels.inp a matrix of the number of linguistic terms of input variables
convert.IntToBin <- function(data.mat, num.labels.inp){
	
	data.mat <- scale.RuleToStr(data.mat, num.labels.inp)	
	n.labels.inp <- num.labels.inp[1,1]		
	popu.var <- matrix(1, nrow = nrow(data.mat), ncol = ncol(data.mat))
	
	seq.k <- seq(0,(ncol(num.labels.inp)-1))
	data.m <- t(apply(data.mat, 1, function(x) x + num.labels.inp * seq.k))	
	popu.val <- matrix(0, nrow = nrow(data.m), ncol = (n.labels.inp * ncol(data.m)))
	f.popu <- function(i,j,data.m){
		popu.val[i, data.m[i,j]] <<- 1
	}
	vec.f.popu <- Vectorize(f.popu, vectorize.args=list("i","j"))
	outer(1 : nrow(data.m), 1 : ncol(data.m), vec.f.popu, data.m)
	
	res <- list(popu.var = popu.var, popu.val = popu.val)	
	return(res)
}

# This function is to change scale of rules
# @param rule.ori a matrix of original rule 
# @param num.labels a matrix of the number of linguistic terms
scale.RuleToStr <- function(rule.ori, num.labels){
	
	seq.k <- seq(0,(ncol(num.labels)-1))
	num <- num.labels * seq.k
	res <- t(apply(rule.ori, 1, function(x) x - sign(x) * num))
	return(res)
}

# This function is to change scale of rule
#
# @param rule.ori a matrix of rules
# @param num.labels a matrix of the number of linguistic terms
scale.StrToRule <- function(rule.ori, num.labels){	
	seq.k <- seq(0,(ncol(num.labels)-1))
	num <- num.labels * seq.k
	res <- t(apply(rule.ori, 1, function(x) x + sign(x) * num))
	
	return(res)
}

# This function is to calculate individual fitness
# @param data.train a matrix of training data
# @param data.sample a matrix of training data which have the same class on consequent part
# @param rule a matrix of rules
# @param grade.cert a matrix of certainty factor
# @param varinp.mf a matrix of parameter values of membership functions
# @param num.labels a matrix of the number of linguistic terms
# @param popu.var a matrix of population describing the selected variables
# @param method.type a kind of method will be used.
# @param k.lower a lower bound of the noise threshold
# @param k.upper a upper bound of the noise threshold
# @param epsilon a value of covering factor
# @param params a list of other parameters
eval.indv.fitness <- function(data.train, rule, method.type = "SLAVE", data.sample = NULL, grade.cert = NULL, varinp.mf = NULL, 
			num.labels = NULL, popu.var = NULL, k.lower = NULL, k.upper = NULL, epsilon = NULL, params = list()){
	if (method.type == "SLAVE") {
		data.input <- data.train[, -ncol(data.train), drop = FALSE]
		res <- matrix()
		
		## degree of completeness
		temp <- degree.completeness(data.sample, rule, varinp.mf, num.labels, popu.var)		
		deg.completeness <- temp$deg.completeness
		n.plus <- temp$n.plus
		
		## soft consistency condition
		soft.cons <- soft.consistency(data.train, rule, varinp.mf, num.labels, n.plus, popu.var, k.lower, k.upper)
		
		indv.fit <- deg.completeness * soft.cons
	}
	
	else if (method.type == "GFS.GCCL"){
		data.input <- data.train[, -ncol(data.train), drop = FALSE]
		data.class <- data.train[, ncol(data.train), drop = FALSE]
		rule.ant <- rule[, -ncol(rule), drop = FALSE]
		degree <- matrix(nrow=nrow(rule.ant), ncol=ncol(rule.ant))
		NCP <- matrix(0, nrow=nrow(rule.ant), ncol = 1)
		
		for (i in 1 : nrow(data.input)){
			res <- det.class(data.input[i, ,drop = FALSE], rule, grade.cert, varinp.mf, type = "WEIGHTED")		
			if(res$res.class == data.train[i, ncol(data.train)]){		
				NCP[res$ind.rule.act, 1] = NCP[res$ind.rule.act, 1] + 1
			}
		}
		
		indv.fit <- NCP
	
	}
	
	else if (method.type == "GFS.FR.MOGUL"){
		data.test <- data.train[, -ncol(data.train), drop = FALSE]
		real.val <- data.train[, ncol(data.train), drop = FALSE]
		ind.fit <- matrix(nrow = nrow(rule), ncol = 1)
		cov.deg <- matrix(nrow = nrow(rule), ncol = 1)
		R <- matrix(nrow = nrow(rule), ncol = nrow(data.train))
		for (i in 1 : nrow(rule)){
			mod <- list(rule = rule[i, ,drop = FALSE])
			res <- GFS.FR.MOGUL.test(mod, data.test)
			residuals <- (real.val - res)
			RMSE <- sqrt(mean(residuals^2))		
			ind.fit[i, 1] <- 1/(1 + RMSE)	
		}
		
		### check high-frequency value
		for(i in 1 : nrow(rule)){
			for (j in 1 : nrow(data.train)){
				R[i, j] <- ch.cover(matrix(data.train[j, ], nrow = 1), rule[i, ,drop = FALSE])
				if (R[i, j] >= epsilon){
					R[i, j] <- 1
				}
				else {
					R[i, j] <- 0
				}
			}
		}
		
		## coverage degree
		cov.deg <- matrix(rowSums(R)/nrow(data.train), ncol = 1)	
		
		## membership function width
		MWR <- MF.width(rule, ncol(data.train))	

		## fitness
		indv.fit <- (ind.fit * cov.deg * MWR)
	}
	
	else if (method.type == "GFS.LT.RS"){
		rule.data.num <- rule
		popu <- popu.var
		var.mf <- varinp.mf
		mode.tuning <- params$mode.tuning
		rule.selection <- params$rule.selection
		num.rule <- nrow(rule.data.num)
		indv.fit <- matrix(nrow = nrow(popu), ncol = 1)
		
		## Calculate fitness
		for (i in 1 : nrow(popu)){					
			if (rule.selection == TRUE){
				rule.data.num.rs <- check.active.rule(rule.data.num, num.rule, popu[i, ,drop = FALSE])	
			} else {
				rule.data.num.rs <- rule.data.num				
			}
			params$popu <- popu[i, ,drop = FALSE]	
			params$var.mf <- var.mf
			params$num.rule <- num.rule
			params$mode.tuning <- mode.tuning
			params$rule.selection <- rule.selection
			res <- tune.MF(method.type = "GFS.LT.RS", params = params)
			var.mf <- res$var.mf
			var.mf.tune <- res$var.mf.tune
			
			## check fitness
			mod <- list(rule.data.num = rule.data.num.rs, var.mf = var.mf, var.mf.tune = var.mf.tune, mode.tuning = mode.tuning, 
			            rule.selection = rule.selection, num.labels = num.labels)
			indv.fit[i, 1] <- calc.MSE(mod, data.train, method.type = "GFS.LT.RS")							
		}
	}
	return(indv.fit)
}

# This function calculates the soft consistency of rules
#
# @param data.train a matrix of training data
# @param rule a matrix of rules
# @param varinp.mf a matrix of parameter values of membership functions
# @param num.labels a matrix of the number of linguistic terms
# @param n.plus a value of the number of positive examples
# @param popu.var a matrix of population describing the selected variables
# @param k.lower a lower bound of the noise threshold
# @param k.upper a upper bound of the noise threshold
soft.consistency <- function(data.train, rule, varinp.mf, num.labels, n.plus, popu.var, k.lower = 0.25, k.upper = 0.75){
	ant.rule <- rule[, -ncol(rule), drop = FALSE]
	
	data.input <- data.train[, -ncol(data.train), drop = FALSE]
	num.labels.inp <- num.labels[, -ncol(num.labels), drop = FALSE]
	
	comp.rule.dt.tra <- create.rule(data.train, varinp.mf, num.labels.inp)
	ant.rule.dt.tra <- comp.rule.dt.tra[, -c(ncol(comp.rule.dt.tra)), drop = FALSE]

	
	n.negative <- matrix(nrow = nrow(rule), ncol = 1)
	
	for (i in 1 : nrow(ant.rule)){
		counter = 0
		ant.rule.i <- ant.rule[i, ,drop = FALSE]
		popu.var.i <- popu.var[i, ,drop = FALSE]
		for (j in 1 : nrow(ant.rule.dt.tra)){
			ant.rule.dt.tra.i <- ant.rule.dt.tra[j, ,drop = FALSE]
			
			similarity <- ch.similarity.rule(ant.rule.i, ant.rule.dt.tra.i, popu.var.i)
			if (similarity == 0){
				if (comp.rule.dt.tra[j, ncol(comp.rule.dt.tra)] != rule[i, ncol(rule)]){
					counter = counter + 1 
				}
			}
		}
		
		n.negative[i, 1] <- counter
	}
	
	k1 = k.lower
	k2 = k.upper
	s.cons <- matrix(nrow = nrow(n.plus), ncol = 1)
	
	for (i in 1 : nrow(n.negative)){
		lower <- k1 * n.plus[i, 1]
		upper <- k2 * n.plus[i, 1]
		if (n.negative[i, 1] <= lower){
			s.cons[i, 1] = 1
		}
		else if (n.negative[i, 1] <= upper){
			if (n.plus[i, 1] == 0){
				s.cons[i, 1] = 0
			} else {
				s.cons[i, 1] = (k2 * n.plus[i, 1] - n.negative[i, 1]) / (n.plus[i, 1] * (k2 - k1))
			}
		} else
			s.cons[i, 1] = 0
	}
	return(s.cons)
}

# This function is to calculate the degree of completeness
#
# @param data.sample a matrix of training data which have the same class on consequent part
# @param rule a matrix of rules
# @param varinp.mf a matrix of the parameter values of membership function
# @param num.labels a matrix of the number of linguistic terms
# @param popu.var a matrix of population describing the selected variable
degree.completeness <- function(data.sample, rule, varinp.mf, num.labels, popu.var){
	ant.rule <- rule[, -ncol(rule), drop = FALSE]
	
	data.input <- data.sample[, -ncol(data.sample), drop = FALSE]
	num.labels.inp <- num.labels[, -ncol(num.labels), drop = FALSE]	
	comp.rule.dt.tra <- create.rule(data.sample, varinp.mf, num.labels.inp)
	
	ant.rule.dt.tra <- comp.rule.dt.tra[, -c(ncol(comp.rule.dt.tra)), drop = FALSE]
	n.plus <- matrix(nrow = nrow(ant.rule), ncol = 1)
	
	for (i in 1 : nrow(ant.rule)){
		counter = 0
		ant.rule.i <- ant.rule[i, ,drop = FALSE]
		popu.var.i <- popu.var[i, ,drop = FALSE]
		for (j in 1 : nrow(ant.rule.dt.tra)){
			ant.rule.dt.tra.i <- ant.rule.dt.tra[j, ,drop = FALSE]
			similarity <- ch.similarity.rule(ant.rule.i, ant.rule.dt.tra.i, popu.var.i)
			if (similarity == 0){
				counter = counter + 1 
			}
		}
		n.plus[i, 1] <- counter
	}
	
	deg.completeness <- n.plus / nrow(data.input)

	res <- list(deg.completeness = deg.completeness, n.plus = n.plus)
	
	return(res)
}

# this function is to count the similarity between two rules
#
# @param rule a matrix of the first rule
# @param rule.data a matrix of the second rule
# @param popu.var a matrix of population describing the selected variables
# @param type a type of similarity
# @param indv.fit a matrix of individual fitness
ch.similarity.rule <- function(rule, rule.data, popu.var = NULL, type = NULL, indv.fit = NULL){
	if (is.null(type)){
		sim = 0
		if (ncol(rule) != ncol(rule.data)){
			stop("the dim is not the same, please check rules")
		}
		
		for (i in 1 : ncol(rule.data)){
			if (popu.var[1, i] == 0){
				count = 0
			} else{
				if (rule[1, i] == rule.data[1, i]){
					count = 0
				} else {
					count = 1
				}
			}
			sim = sim + count
		}
		res <- sim
	} else {
		actual.rule <- rule
		rule.data.num <- rule.data
		nondup.indx <- which(duplicated(actual.rule) == FALSE, arr.ind = TRUE)
		rule.data.num <- rule.data.num[nondup.indx, ,drop = FALSE]
		indv.fit <- indv.fit[nondup.indx, 1, drop = FALSE]
		res <- list(rule = rule.data.num, fit = indv.fit)
	} 
	return(res)
}

# This function is to create rule
# 
# @param data.train a matrix of training data
# @param varinput.mf a matrix of parameter values of membership function of input variables
# @param num.labels a matrix of the number of linguistic terms
create.rule <- function(data.train, varinput.mf, num.labels){
	
	data.input <- data.train[, -ncol(data.train), drop = FALSE]
	num.inputvar <- ncol(data.input)
	
	## fuzzification of training data
	MF <- fuzzifier(data.input, num.inputvar, num.labels, varinput.mf)
	
	MF.max <- matrix(0, nrow = nrow(MF), ncol = ncol(MF))
	k <- 1
	for (i in 1 : length(num.labels)){
		start <- k
		end <- start + num.labels[1, i] - 1
		MF.temp <- MF[, start : end]
		
		for (m in 1 : nrow(MF)){
			max.MFTemp <- max(MF.temp[m, ], na.rm = TRUE)
			max.loc <- which.max(MF.temp[m, ])
			MF.max[m, k + max.loc - 1] <- max.MFTemp		
		}
		k <- k + num.labels[1, i]	
	}
	
	rule.matrix <- MF.max
	rule.matrix[which(rule.matrix > 0)] <- 1
	
	## rule.data.num is the numbers representing the sequence of string of variable names 	
	## delete incomplete rule
	ind.incomplete <- which(rowSums(rule.matrix) != num.inputvar)
	rule.matrix[ind.incomplete, ] <- NA	
	data.class <- data.train[, ncol(data.train), drop = FALSE]
	data.class[ind.incomplete, ] <- NA
	
	## create rule in numeric
	seqq <- seq(1:ncol(rule.matrix))
	rule.data.num <- t(apply(rule.matrix, 1, function(x) x * seqq))	
	rule.data.num <- na.omit(rule.data.num)
	data.class <- na.omit(data.class)
	
	## rule.data.num is the numbers representing the sequence of string of variable names 	
	rule.data.num[which(rule.data.num == 0)] <- NA
	rule.data.num <- t(apply(rule.data.num, 1, na.omit))
	
	comp.rule = cbind(rule.data.num, data.class)
	comp.rule <- comp.rule[which(duplicated(comp.rule) == FALSE), ,drop = FALSE]
	
	return(comp.rule)
}

# This function is to make a partition of membership funciton
#
# @param range.data a matrix of range of data
# @param num.labels a matrix of the number of linguistic terms
# @param type.mf a type of shape of membership function
partition.MF <- function(range.data, num.labels, type.mf = "TRIANGLE"){
	## initialize
	jj <- 0
	ncol.range.data <- ncol(range.data)
	var.mf <- matrix(0, nrow = 5, ncol = sum(num.labels))

	## loop for all variable
	for (i in 1 : ncol.range.data){
		## initialize
		seg <- num.labels[1, i]
		kk <- 1
			
		## Make the depth on each linguistic value, assumed it's similar on all region. 
		delta.point <- (range.data[2, i] - range.data[1, i]) / (seg + 1)
		
			##contruct matrix of parameter of membership function (var.mf) on each variable (max(var) is ncol.data.train
			for (j in 1 : num.labels[1, i]){
			
				##counter to continue index of column var.mf					
				jj <- jj + 1
				
				## type.mf <- 1 means using trapezoid and triangular
				if (type.mf == 1 || type.mf == "TRIANGLE") {
				
					delta.tri <- (range.data[2, i] - range.data[1, i]) / (seg - 1)
					## on the left side
					if (kk %% num.labels[1, i] == 1){   
						var.mf[1, jj] <- 1	
						var.mf[2, jj] <- range.data[1, i]
						var.mf[3, jj] <- var.mf[2, jj]
						var.mf[4, jj] <- var.mf[3, jj] + delta.tri
						var.mf[5, jj] <- NA						
					}
					
					## on the right side
					else if (kk %% num.labels[1, i] == 0){
						var.mf[1, jj] <- 1	
						var.mf[4, jj] <- range.data[2, i]
						var.mf[3, jj] <- var.mf[4, jj] 
						var.mf[2, jj] <- var.mf[3, jj] - delta.tri
						var.mf[5, jj] <- NA						
					}
					
					## on the middle
					else{
						var.mf[1, jj] <- 1	
						var.mf[3, jj] <- range.data[1, i] + (j - 1) * delta.tri
						var.mf[2, jj] <-  var.mf[3, jj] - delta.tri
						var.mf[4, jj] <- var.mf[3, jj] + delta.tri
						var.mf[5, jj] <- NA						
					}	
				}
				## type.mf == 2 means we use trapezoid
				else if (type.mf == 2 || type.mf == "TRAPEZOID") {
					delta.tra <- (range.data[2, i] - range.data[1, i]) / (seg + 2)
					
					## on the left side
					if (kk %% num.labels[1, i] == 1){   
						var.mf[1, jj] <- 2	
						var.mf[2, jj] <- range.data[1, i]
						var.mf[3, jj] <- var.mf[2, jj] + delta.tra
						var.mf[4, jj] <- var.mf[3, jj] + delta.tra
						var.mf[5, jj] <- NA												
					}
					
					## on the right side
					else if (kk %% num.labels[1, i] == 0){
						var.mf[1, jj] <- 3	
						var.mf[4, jj] <- range.data[2, i]
						var.mf[3, jj] <- var.mf[4, jj] - delta.tra
						var.mf[2, jj] <- var.mf[3, jj] - delta.tra
						var.mf[5, jj] <- NA												
					}
					
					## on the middle
					else{
						var.mf[1, jj] <- 4	
						var.mf[2, jj] <- range.data[1, i] + (j - 1) * 1.15 * delta.tra
						var.mf[3, jj] <- var.mf[2, jj] + delta.tra
						var.mf[4, jj] <- var.mf[3, jj] + 0.5 * delta.tra					
						var.mf[5, jj] <- var.mf[4, jj] + delta.tra											
					}
				}
				##Type 5: Gaussian
				##parameter=(mean a, standard deviation b)
				##a=var.mf[2,]
				##b=var.mf[3,]
				else if (type.mf == 3 || type.mf == "GAUSSIAN") {
					delta.gau <- (range.data[2, i] - range.data[1, i]) / (seg - 1)
					##On the left side
					if (kk %% num.labels[1, i] == 1){   
						var.mf[1, jj] <- 5	
						var.mf[2, jj] <- range.data[1, i]
						var.mf[3, jj] <- 0.35 * delta.gau
						var.mf[4, jj] <- NA
						var.mf[5, jj] <- NA
					}
					##On the right side
					else if (kk %% num.labels[1, i] == 0){
						var.mf[1, jj] <- 5	
						var.mf[2, jj] <- range.data[2, i]
						var.mf[3, jj] <- 0.35 * delta.gau
						var.mf[4, jj] <- NA
						var.mf[5, jj] <- NA
					}
					## On the middle
					else {
					var.mf[1, jj] <- 5	
					var.mf[2, jj] <- range.data[1, i] + (j - 1) * delta.gau
					var.mf[3, jj] <- 0.35 * delta.gau
					var.mf[4, jj] <- NA
					var.mf[5, jj] <- NA
					}
				}
				##Type 6: Sigmoid/logistic
				##parameter=(gamma,c)
				##gamma=var.mf[2,]
				##c=var.mf[3,]
				else if (type.mf == 4 || type.mf == "SIGMOID") {
					delta.sig <- (range.data[2, i] - range.data[1, i]) / (seg + 1)
					##On the left side
					if (kk %% num.labels[1, i] == 1){   
						var.mf[1, jj] <- 6	
						var.mf[2, jj] <- range.data[1, i]
						var.mf[3, jj] <- delta.sig
						var.mf[4, jj] <- NA
						var.mf[5, jj] <- NA
					}
					##On the right side
					else if (kk %% num.labels[1, i] == 0){
						var.mf[1, jj] <- 6	
						var.mf[2, jj] <- range.data[2, i]
						var.mf[3, jj] <- delta.sig
						var.mf[4, jj] <- NA
						var.mf[5, jj] <- NA
					}
					## On the middle
					else {
					var.mf[1, jj] <- 6	
					var.mf[2, jj] <- range.data[1, i] + (j - 1) * delta.sig
					var.mf[3, jj] <- delta.sig
					var.mf[4, jj] <- NA
					var.mf[5, jj] <- NA
					}
				}
				##checking for type 7: Generalized Bell
				##parameter=(a,b,c)
				##a=var.mf[2,]
				##b=var.mf[3,]
				##c=var.mf[4,]
				else if (type.mf == 5 || type.mf == "BELL") {
					##On the left side
					if (kk %% num.labels[1, i] == 1){   
						var.mf[1, jj] <- 7	
						var.mf[2, jj] <- 0.6 * delta.point
						var.mf[3, jj] <- 0.6 * delta.point
						var.mf[4, jj] <- range.data[1, i]
						var.mf[5, jj] <- NA
					}
					##On the right side
					else if (kk %% num.labels[1, i] == 0){
						var.mf[1, jj] <- 7	
						var.mf[2, jj] <- 0.6 * delta.point
						var.mf[3, jj] <- 0.6 * delta.point
						var.mf[4, jj] <- range.data[2, i]
						var.mf[5, jj] <- NA
					}
					## On the middle
					else {
					
					var.mf[1, jj] <- 7	
					var.mf[2, jj] <- 0.6 * delta.point
					var.mf[3, jj] <- 0.6 * delta.point
					var.mf[4, jj] <- range.data[1, i] + delta.point * j
					var.mf[5, jj] <- NA
					}
				}						
				kk <- kk + 1
			}
	}

return (var.mf)
}

# This function is used to cut the rules which is redundant with others
# 
# @param rule.data.num a matrix of rules in integer form
# @param num.labels a matrix of the number of linguistic terms
# @param method a kind of method which will be used.
# @param num.inputvar a number of input variables
# @param indv.fit a matrix of individual fitness of rules.
# @param eps a value representing threshold number
prune.rule <- function(rule.data.num, method = "THRIFT", num.labels = NULL, 
                       num.inputvar = NULL, indv.fit = NULL, eps = 0.1){
	if (method == "THRIFT") {
		for (i in 1 : nrow(rule.data.num)){
			min.labels.cons <- (num.labels[1,1] * (ncol(rule.data.num) - 1) + 1)
			max.labels.cons <- min.labels.cons  + num.labels[1,1] - 1 
			
			if ((rule.data.num[i, ncol(rule.data.num)] > max.labels.cons) || (rule.data.num[i, ncol(rule.data.num)] < min.labels.cons)){
				rule.data.num[i, ncol(rule.data.num)] <- as.integer(runif(1, min = min.labels.cons, max = max.labels.cons))
			}	
		}	
		
		res <- unique(rule.data.num)
	}
	
	else if (method == "GCCL"){
		rule.mat <- rule.data.num
		grade.cert <- indv.fit
		
		## check misclass
		indx.mis <- which(rule.mat[, ncol(rule.mat)] == 0)
		for (i in indx.mis){
			rule.mat[i, ncol(rule.mat)] = NA
			grade.cert[i, 1] = NA
		}
		
		rule.mat <- na.omit(rule.mat)
		grade.cert <- na.omit(grade.cert)
		
		if (nrow(rule.mat) > 2){
			temp <- ch.similarity.rule(rule = rule.mat, rule.data = rule.mat, type = "GENERAL", indv.fit = grade.cert)
			rule <- temp$rule
			grade.cert <- temp$fit
		}
		else {
			rule <- rule.mat			
		}
		
		res <- list(rule = rule, grade.cert = grade.cert)	
	}
	return(res)
}

# This function is to calculate heuristic values on consequent part
#
# @param mod a list of parameters values
# @param data.train a matrix of training data
# @param num.labels a matrix of the number of linguistic terms
# @param type.implication.func a type of implication function
calc.heu <- function(mod, data.train, num.labels, type.implication.func){
	## process to evaluate the fitness
	rule <- rulebase(mod$type.model, mod$rule)
	
	## fuzzifier in training phase (for all input and output variables)
	num.varinput <- mod$num.varinput + 1
	num.fvalinput <- cbind(mod$num.fvalinput, num.labels[1,1])
	varinp.mf <- cbind(mod$varinp.mf, mod$varout.mf)
	MF <- fuzzifier(data.train, num.varinput, num.fvalinput, varinp.mf)
	
	ncol.MF <- ncol(MF)
	names.variable <- c(mod$names.varinput, mod$names.varoutput)	
	names.var <- names.variable[1 : ncol.MF]
	colnames(MF) <- c(names.var)
		
	## split fuzzifier result for input variables
	MF.input <- MF[, 1 : (ncol(MF) - num.labels[1,1])]	
	miu.rule <- inference(MF.input, rule, mod$names.varinput, mod$type.tnorm, mod$type.snorm)
		
	## calculate degree	
	f.degree <- function(i, j, miu.rule){	
		degree.cons <- MF[i, mod$rule[j, ncol(mod$rule)]]
		miu.rule[i, j] <<- calc.implFunc(miu.rule[i, j], degree.cons, type.implication.func)
	}	
	vec.f.degree <- Vectorize(f.degree, vectorize.args=list("i","j"))
	outer(1 : nrow(miu.rule), 1 : ncol(miu.rule), vec.f.degree, miu.rule)	
	
	## calculate coverage
	cover.val <- matrix(apply(miu.rule, 2, max))
	
	return(cover.val)	
}

# This function is to calculate compatibility degrees. 
#
# @title The compability degree
# @param data.train.i a matrix(1 x m) of the normalized data.
# @param rule.gen chromosome which represents fuzzy IF-THEN rules.  
# @return The compability degree
# @export
ch.cover <- function(data.train.i, rule.gen){
## calculate degree of MF using chromosome ==> calc.compDegree
comp.degree <- matrix(nrow = nrow(rule.gen))

for (i in 1 : nrow(rule.gen)){
	# do t-norm over all variables
	degree.var <- calc.compDegree(data.train.i, rule.gen[i, ,drop = FALSE], method.type = "GFS.FR.MOGUL")
	comp.degree[i, 1] <- min(degree.var)
}		
return(comp.degree)	
}

# This function is to calculate membership function (MF) degrees. 
#
# @title The MF degree
# @param data.i a matrix(1 x m) of the normalized data.
# @param rule.gen.i a chromosome which represents fuzzy IF-THEN rules. 
# @param method.type a name of method 
# @param var.mf a matrix of membership function parameters
# @param params a list of other parameters
# @return The degree
# @export
calc.compDegree <- function(data.i, rule.gen.i, method.type = "GFS.FR.MOGUL", var.mf = NULL, params = list()){
	## calculate degree of MF
	if (method.type == "GFS.FR.MOGUL"){
		ncol.data.i <- ncol(data.i)
		degree.R <- matrix(nrow = 1, ncol = ncol.data.i)
		for (j in 1 : ncol.data.i){
			xx <- data.i[1, j]
			aa <- rule.gen.i[1, (j * 3 - 2)]
			bb <- rule.gen.i[1, (j * 3 - 1)]
			cc <- rule.gen.i[1, (j * 3)]
			
			if (xx <= aa){
				degree.R[1, j] <- 0
			}
			else if (xx <= bb) {
				degree.R[1, j] <- (xx - aa) / (bb - aa)
			}
			else if (xx <= cc) {
				degree.R[1, j] <- (cc - xx) / (cc - bb)
			} else {
				degree.R[1, j] <- 0
			}
		}
	}
	
	else if (method.type == "GFS.THRIFT"){
		ncol.data.i <- ncol(data.i)
		degree.R <- matrix(nrow = 1, ncol = ncol.data.i)
		for (j in 1 : ncol.data.i){
			if (rule.gen.i[1, j] == 0){
				degree.R[1, j] <- 1
			} else {
				xx <- data.i[1, j]
				aa <- var.mf[2, rule.gen.i[1, j]]
				bb <- var.mf[3, rule.gen.i[1, j]]
				cc <- var.mf[4, rule.gen.i[1, j]]
				if (xx <= aa){
					degree.R[1, j] <- 0
				}
				else if (xx <= bb) {
					degree.R[1, j] <- (xx - aa) / (bb - aa)
				}
				else if (xx <= cc) {
					degree.R[1, j] <- (cc - xx) / (cc - bb)
				} else {
					degree.R[1, j] <- 0
				}
			}
		}
	}
	else if (method.type == "GFS.LT.RS"){
		ncol.data.i <- ncol(data.i)
		degree.R <- matrix(nrow = 1, ncol = ncol.data.i)
		mode.tuning <- params$mode.tuning
		if (mode.tuning == "GLOBAL") {
			for (j in 1 : ncol.data.i){
				if (rule.gen.i[1, j] == 0){
					degree.R[1, j] <- 1
				} else {
					xx <- data.i[1, j]
					aa <- var.mf[2, rule.gen.i[1, j]]
					bb <- var.mf[3, rule.gen.i[1, j]]
					cc <- var.mf[4, rule.gen.i[1, j]]
					if (xx <= aa){
						degree.R[1, j] <- 0
					}
					else if (xx <= bb) {
						degree.R[1, j] <- (xx - aa) / (bb - aa)
					}
					else if (xx <= cc) {
						degree.R[1, j] <- (cc - xx) / (cc - bb)
					} else {
						degree.R[1, j] <- 0
					}
				}
			}
		} else {
			var.mf.tune <- params$var.mf.tune
			num.labels <- params$num.labels
			seqq <- seq(from = num.labels[1, 1], to = ncol(var.mf.tune), by = num.labels[1, 1])
			var.mf.tune <- var.mf.tune[, -seqq, drop = FALSE]
			j.rule <- params$j

			for (j in 1 : ncol.data.i){
				if (rule.gen.i[1, j] == 0){
					degree.R[1, j] <- 1
				} else {
					xx <- data.i[1, j]
					aa <- var.mf[2, rule.gen.i[1, j]] + var.mf.tune[1, (2 * (j.rule - 1) + j)]
					bb <- var.mf[3, rule.gen.i[1, j]] + var.mf.tune[1, (2 * (j.rule - 1) + j)]
					cc <- var.mf[4, rule.gen.i[1, j]] + var.mf.tune[1, (2 * (j.rule - 1) + j)]
					if (xx <= aa){
						degree.R[1, j] <- 0
					}
					else if (xx <= bb) {
						degree.R[1, j] <- (xx - aa) / (bb - aa)
					}
					else if (xx <= cc) {
						degree.R[1, j] <- (cc - xx) / (cc - bb)
					} else {
						degree.R[1, j] <- 0
					}
				}
			}
		}
	}

return(degree.R)
}


# This function is the internal function of the GFS.FR.MOGUL method to compute the predicted values.  
#
# @title GFS.FR.MOGUL: The testing phase
# @param data.test a matrix(m x n) of data for the testing process, where m is the number of instances and 
# n is the number of input variables.
# @param rule.gen fuzzy IF-THEN rules
# @param method.type a type of method used
# @param var.mf a membership function parameters
# @param params a list of other parameters
# @return A matrix of predicted values.
# @export
calc.pred.val <- function(data.test, rule.gen, method.type = "GFS.FR.MOGUL", var.mf = NULL, params = list()){
	output.val.pred <- matrix(nrow = nrow(data.test), ncol = 1)
	if (method.type == "GFS.FR.MOGUL"){
		for (i in 1 : nrow(data.test)){
			deg.R <- matrix()
			deg.R <- ch.cover(data.test[i, ,drop = FALSE], rule.gen)

			##defuzzification
			##get mean value of output variable on rule.gen
			output.val.mean <- rule.gen[, (ncol(rule.gen) - 1), drop = FALSE]
			
			## calculate average values
			if (sum(deg.R) != 0){
				output.val.pred[i, 1] <- (t(deg.R) %*% output.val.mean) / sum(deg.R)		
			} else {
				output.val.pred[i, 1] <- 0.5
			}
		}
	}
	else if (any(method.type == c("GFS.THRIFT", "GFS.LT.RS"))){
		for (i in 1 : nrow(data.test)){
			deg.R <- matrix(0, nrow = nrow(rule.gen), ncol = 1)
			output.val.mean <- matrix(0.5, nrow = nrow(rule.gen), ncol = 1)
			
			for (j in 1 : nrow(rule.gen)){
				# do t-norm over all variables
				if (method.type == "GFS.LT.RS"){
					params$j = j
				}
				degree.var <- calc.compDegree(data.test[i, ,drop = FALSE], rule.gen[j, ,drop = FALSE], method.type, var.mf = var.mf, params = params)
				deg.R[j, 1] <- min(degree.var)
				fired.term <- rule.gen[j, ncol(rule.gen)]
				output.val.mean[j, 1] <- var.mf[3, fired.term]
			}	
		
			## calculate average values
			if (sum(deg.R) != 0){
				output.val.pred[i, 1] <- (t(deg.R) %*% output.val.mean) / sum(deg.R)
			} else {
				output.val.pred[i, 1] <- 0.5
			}
		}
	}
	
return(output.val.pred)
}

# This function is the internal function of the GFS.FR.MOGUL method to tune MF.  
#
# @title GFS.FR.MOGUL: The step of MF tuning 
# @param method.type a type of method
# @param params a list of parameters based on type of method
# @return a frbs object.
# @export
tune.MF <- function(method.type = "GFS.FR.MOGUL", params = list()){
	if (method.type == "GFS.FR.MOGUL") {
		mod <- params$mod
		data.train <- params$data.train
		max.iter <- params$max.iter
		
		ii <- 1
		old.rule <- mod$rule
		best.rules <- old.rule
		RMSE <- 1000
		while (ii < max.iter){
			new.mod <- list(rule = best.rules)
			new.RMSE <- calc.MSE(new.mod, data.train, method.type = "GFS.FR.MOGUL")
			if (new.RMSE < RMSE){
				mod$rule <- best.rules
				mod$RMSE <- new.RMSE
				RMSE <- new.RMSE
			}
			
			for (i in 1 : nrow(best.rules)){
				for (j in 1 : ncol(best.rules)){
					if (j %% 3 == 1){
						mean.delta <- ((best.rules[i, j + 1] - best.rules[i, j]))/2
						left.bound <- (best.rules[i, j] - mean.delta)
						right.bound <- (best.rules[i, j] + mean.delta)
					}
					else if (j %% 3 == 0){
						mean.delta <- ((best.rules[i, j] - best.rules[i, j - 1]))/2
						left.bound <- (best.rules[i, j] - mean.delta)
						right.bound <- (best.rules[i, j] + mean.delta)
					}
					else {
						left.bound <- (best.rules[i, j] - ((best.rules[i, j] - best.rules[i, j - 1]))/2)
						right.bound <- (best.rules[i, j] + ((best.rules[i, j + 1] - best.rules[i, j]))/2)
					}
					
					r.val <- runif(1, min = left.bound, max = right.bound)
					best.rules[i, j] <- r.val
				}
			}
					
			ii <- ii + 1
		}		
	}
	
	else if (method.type == "GFS.LT.RS"){
		rule.selection <- params$rule.selection
		popu <- params$popu
		var.mf <- params$var.mf
		mode.tuning <- params$mode.tuning
		num.rule <- params$num.rule
		## mode.tuning "global"
		if (rule.selection == TRUE){
			if (mode.tuning == "GLOBAL"){
				popu.lateral <- popu[1, 1 : (ncol(popu) - num.rule), drop = FALSE]
				
				## change parameter of triangular
				var.mf.lt <- var.mf
				for (j in 2 : 4){
					var.mf.lt[j, ] <- var.mf[j, , drop = FALSE] + popu.lateral[1, drop = FALSE]					
				}
				mod <- list(var.mf = var.mf.lt, var.mf.tune = NULL)
			} else {
				popu.lateral <- popu[1, 1 : (ncol(popu) - num.rule), drop = FALSE]
				mod <- list(var.mf = var.mf, var.mf.tune = popu.lateral)
				
			}
		} else {
			if (mode.tuning == "GLOBAL"){
				popu.lateral <- popu
				
				## change parameter of triangular
				var.mf.lt <- var.mf
				for (j in 2 : 4){
					var.mf.lt[j, ] <- var.mf[j, , drop = FALSE] + popu.lateral[1, drop = FALSE]					
				}
				mod <- list(var.mf = var.mf.lt, var.mf.tune = NULL)
				
			} else {
				popu.lateral <- popu
				mod <- list(var.mf = var.mf, var.mf.tune = popu.lateral)
				
			}
		}
	}
	return(mod)
}

# This function is the internal function in order to generate a population
#
# @title Population generating function
# @param method.type a type of method
# @param data.train a matrix(m x n) of data for the testing process, where m is the number of instances and 
# n is the number of variables.
# @param num.var a number of variables
# @param num.labels a number of linguistic terms
# @param popu.size a size of population
# @param range.data a matrix representing range of data
# @param type.mf a type of membership function
# @param name.class a value on consequent part as a class
# @return a matrix or list of population
# @export
generate.popu <- function(method.type = "GFS.FR.MOGUL", data.train = NULL, 
                         num.var = NULL, num.labels = NULL, popu.size = NULL, 
						 range.data = NULL, type.mf = NULL, name.class = NULL, params = list()){						 
	if (method.type == "GFS.FR.MOGUL"){
		data.i <- data.train
		## now, this method just implement triangular membership function
		rule.gen <- matrix()
		for (i in 1 : ncol(data.i)){					
			dt.i <- data.i[1, i]
			delta.data <- min(dt.i, 1 - dt.i)
			rand.delta <- runif(1, min = 0, max = delta.data)			
			temp.rule.i <- matrix(c((dt.i - rand.delta), dt.i, (dt.i + rand.delta)), nrow = 1)			
			rule.gen <- cbind(rule.gen, temp.rule.i)
		}
		popu <- rule.gen[1, 2 : ncol(rule.gen)]		
	}
	else if (any(method.type == c("GFS.GCCL", "FH.GBML"))){			
		rule.data.num <- matrix(nrow = popu.size, ncol = num.var)
		for (i in 1 : popu.size){
			k <- 0
			for (j in 1 : num.var){
				rule.data.num[i, j] = as.integer(runif(1, min = 0, max = (num.labels[1,j] + 0.99)))
				
				if (rule.data.num[i, j] != 0){
					rule.data.num[i, j] = rule.data.num[i, j] + (num.labels[1, j] * k)
					k = k + 1
				} else {
					k = k + 1
				}
			}
		}		
		popu <- rule.data.num
	}
	else if (method.type == "SLAVE"){
		data.sample <- data.train[which(data.train[, ncol(data.train)] == name.class), ]
		data.input <- data.sample[, -ncol(data.sample), drop = FALSE]
		## define the number of sample
		num.popu.sample <- 5
		if (num.popu.sample < nrow(data.input)){
			num.popu.sample <- nrow(data.input)
		}
		
		## generate initial model using WM
		mod <- WM(data.input[1 : num.popu.sample, ,drop = FALSE], num.labels, type.mf = "TRIANGLE", type.tnorm = "MIN", type.implication.func = "ZADEH")
		
		## antecedent part on rules
		ant.rule <- mod$rule.data.num
		
		## parameter values of MF of input and output variables
		temp1.mf <- mod$varinp.mf
		temp2.mf <- mod$varout.mf
		varinp.mf <- cbind(temp1.mf, temp2.mf)	
		
		## complete rules: antecedent and consequent part
		comp.rule <- cbind(ant.rule, name.class)
		
		## population of VAR
		## It defines whether the variables have linguistic values or don't care
		popu.var <- matrix(1, nrow = nrow(comp.rule), ncol = num.var)
		popu <- list(popu.var = popu.var, popu.val = comp.rule, varinp.mf = varinp.mf)
	}
	else if (method.type == "GFS.LT.RS"){
		mode.tuning <- params$mode.tuning
		rule.selection <- params$rule.selection
		type.tnorm <- params$type.tnorm
		type.implication.func <- params$type.implication.func
		## generate initial rules and membership function from  Wang & Mendel's method
		type.mf = "TRIANGLE"
		mod.init <- NULL
		mod.init <- WM(data.train, num.labels, type.mf, type.tnorm = type.tnorm, type.implication.func = type.implication.func, 
		               classification = FALSE, range.data = range.data)
		rule.data.num <- mod.init$rule.data.num
		rule <- mod.init$rule
		num.rule <- nrow(rule)
		varinp.mf <- mod.init$varinp.mf
		varout.mf <- mod.init$varout.mf
		var.mf <- cbind(varinp.mf, varout.mf)
		
		## manipulate points at left and right side of triangular
		seqq.left <- seq(from = 1, to = ncol(var.mf), by = num.labels[1, 1])
		seqq.right <- seq(from = num.labels[1, 1], to = ncol(var.mf), by = num.labels[1, 1])
		var.mf[2, seqq.left] <- var.mf[2, seqq.left] - 0.5
		var.mf[4, seqq.right] <- var.mf[4, seqq.right] + 0.5
		
		## coding scheme and initial gene pool as population based on model.tuning and rule.selection
		if (mode.tuning == "GLOBAL") {
			num.col.lateral <- sum(num.labels)
			if (rule.selection == TRUE){
				popu.lateral <- matrix(runif(num.col.lateral * popu.size, min = -0.5, max = 0.5), nrow = popu.size, ncol = num.col.lateral)
				popu.lateral <- rbind(matrix(rep(0, num.col.lateral), nrow = 1), popu.lateral)
				popu.rule <- matrix(round(runif(popu.size * num.rule, min = 0, max = 1)), nrow = popu.size, ncol = num.rule)
				popu.rule <- rbind(matrix(rep(1, num.rule), nrow = 1), popu.rule)
				popu <- cbind(popu.lateral, popu.rule)
			} else {
				popu.lateral <- matrix(runif(num.col.lateral * popu.size, min = -0.5, max = 0.5), nrow = popu.size, ncol = num.col.lateral)
				popu.lateral <- rbind(matrix(rep(0, num.col.lateral), nrow = 1), popu.lateral)
				popu <- popu.lateral
			}
		}
		else if (mode.tuning == "LOCAL") {
			num.col.lateral <- num.rule * ncol(num.labels)
			if (rule.selection == TRUE){
				popu.lateral <- matrix(runif(num.col.lateral * popu.size, min = -0.5, max = 0.5), nrow = popu.size, ncol = num.col.lateral)
				popu.lateral <- rbind(matrix(rep(0, num.col.lateral), nrow = 1), popu.lateral)
				popu.rule <- matrix(round(runif(popu.size * num.rule, min = 0, max = 1)), nrow = popu.size, ncol = num.rule)
				popu.rule <- rbind(matrix(rep(1, num.rule), nrow = 1), popu.rule)
				popu <- cbind(popu.lateral, popu.rule)
			} else {
				popu.lateral <- matrix(runif(num.col.lateral * popu.size, min = -0.5, max = 0.5), nrow = popu.size, ncol = num.col.lateral)
				popu.lateral <- rbind(matrix(rep(0, num.col.lateral), nrow = 1), popu.lateral)
				popu <- popu.lateral
			}
		}
		
		popu <- list(var.mf = var.mf, popu = popu, num.rule = num.rule, rule.data.num = rule.data.num)
	}
	return (popu)
}

# This function is the internal function in order to calculate small membership function width
#
# @title Membership function width
# @param rule a set of rules
# @param num.var a number of variables
MF.width <- function(rule, num.var){
	DW <- 1
	WVR <- matrix(nrow = nrow(rule), ncol = num.var)
	MWR <- matrix(nrow = nrow(rule), ncol = 1)
	a <- 1
	for (i in 1 : nrow(rule)){
		for (j in 1 : num.var){
			k <- 3 * (j - 1) + 1
			WVR[i, j] <- rule[i, k] - rule[i, k + 2]
		}	
	}	
	RW <- (rowSums(WVR)/DW)/num.var 
	MWR <- as.matrix(exp(-abs(1 - a * RW)), ncol = 1)
	return(MWR)
}





