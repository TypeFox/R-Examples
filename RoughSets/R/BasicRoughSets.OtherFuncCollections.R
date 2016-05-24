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

# This function is used to calculate the implication function
#
# @title A implication function
# @param antecedent a degree of antecedent part
# @param consequent a degree of consequent part
# @param type.implication.func a type of implication function
calc.implFunc <- function(antecedent, consequent, type.implication.func = "lukasiewicz"){

	if (class(type.implication.func) == "character"){
		if (type.implication.func == "kleene_dienes"){
			temp.rule.degree <- pmax(1 - antecedent, consequent)
		}
		else if (type.implication.func == "lukasiewicz"){
			temp.rule.degree <- pmin(1 - antecedent + consequent, 1)
		}
		else if (type.implication.func == "zadeh"){
			temp.rule.degree <- pmax((1 - antecedent), min(antecedent, consequent))
		}
		else if (type.implication.func == "gaines"){
			temp.rule.degree <- matrix(nrow = nrow(antecedent), ncol = ncol(antecedent));

			## create function for vectorizing
			func.imp <- function(i, j, antecedent, consequent){
				if (antecedent[i, j] <= consequent[j])
					temp.rule.degree[j, i] <<- 1
				else
					temp.rule.degree[j, i] <<- consequent[j] / antecedent[i,j]
			}
			
			## vectorize func.def.discern
			VecFun <- Vectorize(func.imp, vectorize.args=list("i","j"))
			outer(1 : nrow(antecedent), 1 : ncol(antecedent), VecFun, antecedent, consequent)	
		}
		else if (type.implication.func == "godel"){
			temp.rule.degree <- matrix(nrow = nrow(antecedent), ncol = ncol(antecedent));

			## create function for vectorizing
			func.imp <- function(i, j, antecedent, consequent){
				if (antecedent[i, j] <= consequent[j])
					temp.rule.degree[j, i] <<- 1
				else
					temp.rule.degree[j, i] <<- consequent[j]
			}
			
			## vectorize func.def.discern
			VecFun <- Vectorize(func.imp, vectorize.args=list("i","j"))
			outer(1 : nrow(antecedent), 1 : ncol(antecedent), VecFun, antecedent, consequent)		
		}
		else if (type.implication.func == "kleene_dienes_lukasiewicz"){
			temp.rule.degree <- (1 - antecedent + antecedent * consequent)
		}
		else if (type.implication.func == "mizumoto"){
			temp.rule.degree <- 1 - antecedent + antecedent * consequent
		}
		else if (type.implication.func == "dubois_prade"){
			temp.rule.degree <- matrix(nrow = nrow(antecedent), ncol = ncol(antecedent))
			
			## create function for vectorizing
			func.imp <- function(i, j, antecedent, consequent){
				if (consequent[j] == 0)
					temp.rule.degree[j, i] <- 1 - antecedent[i, j]
				else {
					if (antecedent[i, j] == 1)
						temp.rule.degree[j, i] <<- consequent[j]
					else 
						temp.rule.degree[j, i] <<- 1
				}
			}
			
			## vectorize func.def.discern
			VecFun <- Vectorize(func.imp, vectorize.args=list("i","j"))
			outer(1 : nrow(antecedent), 1 : ncol(antecedent), VecFun, antecedent, consequent)			
		}
	}
	else if (class(type.implication.func) == "function"){
		temp.rule.degree <- matrix(nrow = nrow(antecedent), ncol = ncol(antecedent))
		
		## create function for vectorizing
		func.imp <- function(i, j, antecedent, consequent){
			temp.rule.degree[j, i] <<- type.implication.func(antecedent[i, j], consequent[j])
		}
		
		## vectorize func.def.discern
		VecFun <- Vectorize(func.imp, vectorize.args=list("i","j"))
		outer(1 : nrow(antecedent), 1 : ncol(antecedent), VecFun, antecedent, consequent)		
	}

	return (temp.rule.degree)
}

# This function is used to calculate fuzzy similarity relation using fuzzy tolerance: implicator and t-norm
#
# @title A fuzzy similarity relation
# @param decision.table decision table
# @param attributes variables
# @param t.similarity a type of implication function
# @param t.tnorm a type of t-norm
# @param FUN a user-customized function
# @param disc.mat a boolean value showing whether it produces discernibility matrix completely or not
# @param delta a numeric value for gaussian
# @param alpha a numeric value between [0, 1]
# @param t.aggregation a type of aggregator
# @param type.MVCompletion a boolean value shows whether the decision table contains missing value or not 
calc.fsimilarity <- function(decision.table, attributes, t.similarity = "eq.1", t.tnorm = "lukasiewicz",
                            FUN = NULL, disc.mat = FALSE, delta = NULL, alpha = 0.5, t.aggregation = "t.tnorm", type.MVCompletion = FALSE){
				   
	## get data
	objects <- decision.table
	nominal.att <- attr(decision.table, "nominal.attrs")
	temp <- matrix()
	num.objects <- nrow(objects)
	
	## calculate variance and range (only for real-values attributes)
	res.varRange <- cal.var.range(objects, attributes, nominal.att)
	variance.data <- res.varRange$variance.data
	range.data <- res.varRange$range.data
	
	## initialization which depends on value of disc.mat
	if (disc.mat == TRUE){	
		miu.Ra <- data.frame()
	}
	else {
		if (type.MVCompletion == TRUE){
			miu.Ra <- list()
		}
		else {
			miu.Ra <- matrix(nrow = 1, ncol = num.objects)
		}
	}
	
	## get attributes
	obj <- objects[, c(attributes), drop = FALSE]
	nom <- nominal.att[c(attributes)]	
	all.nom <- unique(nom) %in% TRUE
	t.range.data <- matrix(range.data, ncol = 1)
	t.variance.data <- matrix(variance.data, ncol = 1)

	## disc.mat == FALSE means we do not create discernibility matrix
	if (disc.mat == FALSE){		
		if (all.nom == FALSE && t.aggregation == "kernel.frst") {
			options(warn=-1)			
			distance <- as.matrix(stats::dist(obj, method = "euclidean", diag = TRUE, upper = TRUE))

			miu.Ra.temp <- do.call(t.tnorm, list(distance, delta))
			
			## if the decision table contains missing values
			if (type.MVCompletion == TRUE){
				miu.Ra <- list()
				miu.Ra$lower <- miu.Ra.temp^2
				miu.Ra$upper <- miu.Ra.temp^0.5
				for (j in 1 : length(miu.Ra)){
					miu.Ra[[j]][which(is.na(miu.Ra[[j]]), arr.ind=TRUE)] <- 0
					miu.Ra[[j]][which(is.na(miu.Ra[[j]]), arr.ind=TRUE)] <- 1
				}
			}
			else {
				miu.Ra <- miu.Ra.temp
			}
		}
		else {
			init.IND <- list()
			for (i in 1 : length(attributes)){
				if (nom[i] == TRUE){
					## overwrite type of similarity when crisp
					t.similarity.new <- "boolean"
				}
				else {
					t.similarity.new <- t.similarity
				}
				## calculate indiscernibility relation
				init.IND[[i]] <- outer(c(obj[, i]), c(obj[, i]), 
							 function(x, y) similarity.equation(x, y, t.similarity.new, delta, t.range.data[i], t.variance.data[i], FUN, decision.table, alpha))	
			}
			
			if (type.MVCompletion == TRUE){				
				init.IND.lower <- lapply(init.IND, function(x) x^2)
				init.IND.upper <- lapply(init.IND, function(x) x^0.5)
				for (j in 1 : length(init.IND.lower)){
					init.IND.lower[[j]][which(is.na(init.IND.lower[[j]]), arr.ind=TRUE)] <- 0
					init.IND.upper[[j]][which(is.na(init.IND.upper[[j]]), arr.ind=TRUE)] <- 1
				}	
				## Reduce into single matrix by considering t.tnorm operator			
				miu.Ra$lower <- Reduce.IND(init.IND.lower, t.tnorm, t.aggregation = t.aggregation, delta = delta)
				miu.Ra$upper <- Reduce.IND(init.IND.upper, t.tnorm, t.aggregation = t.aggregation, delta = delta)
			}
			else {
				## Reduce into single matrix by considering t.tnorm operator			
				miu.Ra <- Reduce.IND(init.IND, t.tnorm, t.aggregation = t.aggregation, delta = delta)
			}
		}
	}
	else {
		temp.init <- NULL
		if (type.MVCompletion == TRUE){
			stop("The decision table contains missing values, please estimate them by using missing value preprocessing.")
		}
		for (i in 1 : length(attributes)){
			if (nom[i] == TRUE){
				t.similarity.new <- "boolean"
			}
			else {
				t.similarity.new <- t.similarity
			}
			## construct indiscernibility relation
			temp <- 1 - outer(c(obj[, i]), c(obj[, i]), 
						 function(x, y) similarity.equation(x, y, t.similarity.new, delta, t.range.data[i], t.variance.data[i], FUN, decision.table, alpha))			 
			
			if (is.null(temp.init)){
				temp.init <- as.list(temp)
			}
			else {
				temp.init <- rbind(temp.init, as.list(temp))
			}
		}
		
		## construct discernibility matrix
		miu.Ra <- data.frame()
		k <- 1
		for (i in 1: nrow(objects)){
			for (j in 1 : nrow(objects)){
				miu.Ra[i, j] <- data.frame(I(list(temp.init[,k])))
				k <- k + 1
			}
		}			
	}
	
	if (type.MVCompletion == FALSE){
		rownames(miu.Ra) <- seq(1 : nrow(miu.Ra))
		colnames(miu.Ra) <- seq(1 : nrow(miu.Ra))
	}

	return(miu.Ra)	
}

# It is used to express similarity equation
# @param x objects
# @param temp.obj objects
# @param t.similarity a type of fuzzy similarity equation
# @param delta a numeric value for gaussian method
# @param range.data range of data
# @param variance.data variance of data
# @param FUN a custom function
# @param decision.table a decision table
# @param alpha a numeric value between [0, 1]
similarity.equation <- function(x, temp.obj, t.similarity = "eq.1", delta = 0.2, range.data = NULL, 
                        variance.data = NULL, FUN = NULL, decision.table = NULL, alpha = 1){
	## calculate fuzzy similarity based on type of similarity
	if (!is.null(FUN)){
		temp.miu.Ra <- FUN(decision.table = decision.table, x = x, y = temp.obj)
	}
	
	else if (t.similarity == "eq.1"){							
		## calculate for all column/attributes
		temp.miu.Ra <- c(1 - abs(x - temp.obj) / abs(range.data))
	}
	else if (t.similarity == "eq.2"){	
		## calculate for all column/attributes
		temp.miu.Ra <- exp(-((x - temp.obj)^2) / (2 * variance.data^2))	
	}
	else if (t.similarity == "eq.3"){
		left.part <- (temp.obj - x + variance.data) / variance.data									  
		right.part <- (x - temp.obj +  variance.data) / variance.data
		min.part <- apply(rbind(left.part, right.part), 2, min)
		temp.miu.Ra <- apply(rbind(min.part, 0), 2, max)		
	}
	else if (t.similarity == "eq.instance.selection"){
		## calculate for all column/attributes
		temp <- 1 - alpha * abs(x - temp.obj) /	range.data
		temp.miu.Ra <- apply(rbind(temp, 0), 2, max)
		#temp.miu.Ra <- pmax(0, 1 - alpha * abs(x - temp.obj) /	range.data)
	}
	## fuzzy similarity relations proposed by Eric Tsang, et.al.
	else if (t.similarity == "min"){		
		## calculate for all column/attributes
		if (x == temp.obj){
			temp.miu.Ra <- 1
		}
		else {
			temp.miu.Ra <- min(x, temp.obj)
		}								
	}
	## Kernel fuzzy rough set: gaussian kernel
	else if (t.similarity == "gaussian"){	
		if (is.null(delta)){
			## calculate for all column/attributes
			temp.miu.Ra <- exp(-((x - temp.obj)^2) / (2 * variance.data^2))	
		}
		else {
			## calculate for all column/attributes
			temp.miu.Ra <- exp(- (abs(x - temp.obj)^2)/delta)								
		}
	}
	## Kernel fuzzy rough set: exponential kernel
	else if (t.similarity == "exponential"){											
		## calculate for all column/attributes
		temp.miu.Ra <- exp(- (abs(x - temp.obj))/delta)									
	}

	## Kernel fuzzy rough set: rational quadratic kernel
	else if (t.similarity == "rational"){											
		## calculate for all column/attributes
		temp.miu.Ra <- exp(1 - (abs(x - temp.obj)^2)/ ((abs(x - temp.obj)^2) + delta))									
	}

	## Kernel fuzzy rough set: circular kernel
	else if (t.similarity == "circular"){			
		if (abs(x - temp.obj) >= delta){
			warnings("the condition is not satisfied, please increase the value of delta")
		}
		## calculate for all column/attributes
		temp.miu.Ra <- 2/pi * acos(abs(x - temp.obj)/delta) - 2/pi * 
						(abs(x - temp.obj)/delta) * sqrt(1 - (abs(x - temp.obj)/delta)^2)		
	}			

	## Kernel fuzzy rough set: spherical kernel
	else if (t.similarity == "spherical"){			
		if (abs(x - temp.obj) >= delta){
			warnings("the condition is not satisfied, please increase the value of delta")
		}
		## calculate for all column/attributes
		temp.miu.Ra <- 1 - 3/2 *  abs(x - temp.obj) / delta + 1/2 * 
						(abs(x - temp.obj) / delta)^3		
	}
	## Kernel fuzzy rough set
	else if (any(t.similarity == c("gaussian.kernel", "exponential.kernel", "rational.kernel", "circular.kernel", "spherical.kernel"))){			
			## calculate for all column/attributes
			temp.miu.Ra <- abs(x - temp.obj)/range.data
	}
	else if (t.similarity == "boolean"){
		temp.miu.Ra <- apply(cbind(x, temp.obj), 1, function(x) all(x[1] == x) - 0)
	}
	
	return(temp.miu.Ra)
}

# it is used to calculate variance and range data
# @param objects data frame of objects
# @param attributes attributes are considered
# @param nominal.att a list of attribute types
cal.var.range <- function(objects, attributes, nominal.att){
	variance.data <- c()
	range.data <- c()
	desc.attrs <- attr(objects, "desc.attrs")
	
	for (i in c(attributes)){
		if (nominal.att[i] == FALSE){
			temp.variance.data <- sqrt(stats::var(objects[, i, drop = FALSE], na.rm = TRUE))
			temp.range.data <- c(desc.attrs[[i]][2] - desc.attrs[[i]][1])
			
			## condition to avoid "divide by zero"
			if (temp.range.data < 0.000000001){
				temp.range.data <- 0.000000001
			}
			if (temp.variance.data < 0.000000001){
				temp.variance.data <- 0.000000001
			}
		}
		else {
			temp.variance.data <- NA
			temp.range.data <- NA
		}
		variance.data <- cbind(variance.data, temp.variance.data)
		range.data <- cbind(range.data, temp.range.data)
	}
	
	return(list(variance.data = variance.data, range.data = range.data))
}

# This function is used to calculate t-norm of two variables
#
# @title t-norm calculation on two variables
# @param right.val a value of one attribute on right side
# @param init.val a value of attribute
# @param t.tnorm a type of t-norm
func.tnorm <- function(right.val, init.val, t.tnorm = "min"){

	if (class(t.tnorm) == "character"){
		if (t.tnorm == "lukasiewicz"){
			return (pmax(right.val + init.val - 1, 0))
		}
		else if (t.tnorm == "product"){
			return (right.val * init.val)
		}
		else if (t.tnorm == "yager"){
			return (1 - pmin(1, ((1 - init.val) + (1 - right.val))))
		}
		else if (t.tnorm == "t.cos"){
			return (pmax(init.val * right.val - sqrt(1 - init.val^2)*sqrt(1 - right.val^2), 0))
		}
		else if (t.tnorm == "hamacher"){
			return ((init.val * right.val)/(init.val + right.val - init.val * right.val))
		}
		else {
			return (pmin(right.val, init.val))
		}
	}
	else if (class(t.tnorm) == "function"){
		temp <- matrix(nrow = nrow(init.val), ncol = ncol(init.val))
		
		## create function for vectorizing
		f.tnorm <- function(i, j, init.val, right.val){
			temp[j, i] <<- t.tnorm(init.val[i, j], right.val[j])
		}
		
		## vectorize func.def.discern
		VecFun <- Vectorize(f.tnorm, vectorize.args=list("i","j"))
		outer(1 : nrow(init.val), 1 : ncol(init.val), VecFun, init.val, right.val)	
		return (temp)		
	}
}

# This is a function that implements a part of the fundamental concept of fuzzy rough set theory: positive, negative, boundary and degree of dependency.
# 
# @title Fuzzy region based on fuzzy rough set
#
# @param decision.table a data frame representing decision table
# @param fuzzyroughset a fuzzy rough set. See \code{\link{BC.LU.approximation.FRST}}.
BC.def.region.FRST <- function(decision.table, fuzzyroughset){
	## get the data
	objects <- decision.table
	num.object <- nrow(objects)
	nominal.att <- attr(decision.table, "nominal.attrs")
	if (is.null(attr(decision.table, "decision.attr"))){
		decision.attr <- ncol(objects)
	}
	else {
		decision.attr <- attr(decision.table, "decision.attr")
		## check if index of decision attribute is not in last index, then change their position
		if (decision.attr != ncol(objects)){
			objects <- cbind(objects[, -decision.attr, drop = FALSE], objects[, decision.attr, drop = FALSE])
			nominal.att <- c(nominal.att[-decision.attr], nominal.att[decision.attr])
		}
	}
	
	## get lower and upper apprx.
	fuzzy.lower <- fuzzyroughset$fuzzy.lower
	fuzzy.upper <- fuzzyroughset$fuzzy.upper
	
	if (nominal.att[decision.attr] == FALSE){
		fuzzy.lower <- apply(fuzzy.lower, 1, max)
		fuzzy.upper <- apply(fuzzy.upper, 1, max)
		
		## get positive and negative region values
		positive.reg <- fuzzy.lower
		negative.reg <- 1 - fuzzy.upper
		
		boundary.reg <- matrix()		
		## get boundary region for each decision concept
		for (i in 1 : length(fuzzy.upper)){
			temp <- fuzzy.upper[i] - fuzzy.lower[i]
			#temp.neg <- 1 - fuzzy.upper[i]
			boundary.reg[i] <- temp
		}	
		
		names(boundary.reg) <- c(as.character(unique(objects[, decision.attr])))
		
		## calculate degree of dependency
		degree.depend <- sum(positive.reg)/num.object
		res <- list(positive.freg = positive.reg, negative.freg = negative.reg, boundary.freg = boundary.reg, 
					 degree.dependency = degree.depend)		
	}
	else {
		## get positive and negative region values
		positive.reg <- apply(matrix(unlist(fuzzy.lower), nrow = length(fuzzy.lower), byrow = TRUE), 2, max)	
		negative.reg <- 1 - apply(matrix(unlist(fuzzy.upper), nrow = length(fuzzy.upper), byrow = TRUE), 2, max)
		
		boundary.reg <- c()
		
		## get boundary region for each decision concept
		for (i in 1 : length(fuzzy.upper)){
			temp <- fuzzy.upper[[i]] - fuzzy.lower[[i]]
			#temp.neg <- 1 - fuzzy.upper[[i]]
			boundary.reg <- append(boundary.reg, list(temp))
		}	
		names(boundary.reg) <- c(as.character(unique(objects[, decision.attr])))
		
		## calculate degree of dependency
		degree.depend <- sum(positive.reg)/num.object
		res <- list(positive.freg = positive.reg, negative.freg = negative.reg, boundary.freg = boundary.reg, 
					 degree.dependency = degree.depend)
	}
	return(res)
}


# calculate VQRS
# @param data a value
# @param alpha.beta a parameter representing alpha and beta for some and most.
calc.vqrs.relation <- function(data, alpha.beta = c(0.2, 1)){
	if (data <= alpha.beta[1])
		return (0)
	else if (data >= alpha.beta[1] && data <= ((alpha.beta[1] + alpha.beta[2])/2))
		return (2*(data - alpha.beta[1])^2/(alpha.beta[2] - alpha.beta[1])^2)
	else if (data >= ((alpha.beta[1] + alpha.beta[2])/2) && data <= alpha.beta[2])
		return (1 - (2*(data - alpha.beta[2])^2/(alpha.beta[2] - alpha.beta[1])^2))
	else if (data >= alpha.beta[2])
		return (1)
}

## it is used to calculate OWA
# @param data a dataset
# @param m.owa a value for OWA
# @param type.method a type of OWA: owa, weight, rfrs
# @param type.apprx a type of approximations: lower, upper
# @param type.rfrs a type of robust fuzzy rough set models: k.trimmed.min, k.mean.min, k.median.min, etc.
# @param k.rfrs a value for rfrs
# @param w.owa weighted vector for OWA
calc.OWA <- function(data, m.owa = 2, type.method = "owa", type.apprx = "lower", type.rfrs = "k.trimmed.min", 
                     k.rfrs = 1, w.owa = NULL){
	data <- matrix(data, nrow = 1)
	
	## sort the data
	if (any(type.method == c("owa", "weight"))){
		new.data.sorted <- data[, order(data, decreasing = TRUE)]
	}
	else if (type.method == "rfrs"){
		new.data.sorted <- data[, order(data, decreasing = FALSE)]
	}
	
	## lower/upper approximations based on ordered weighted average (OWA)
	if (any(type.method == c("owa"))){
		if (is.null(w.owa)){
			w <- matrix(0, nrow = length(data))
			for (i in 1 : m.owa){
				w[i] <- (2^(m.owa - i))/(2^m.owa - 1)
			}
		}
		else {
			w <- matrix(sort(w.owa, decreasing = TRUE), ncol = 1)
		}
		if (type.apprx == "lower"){
			w.lower <- rev(w)
			return(sum(new.data.sorted * w.lower))
		}
		else {
			w.upper <- w
			return(sum(new.data.sorted * w.upper))
		}
	}
	## OWA for robust fuzzy rough sets
	else if (type.method == "rfrs"){
		w <- matrix(0, nrow = length(data))
		if (type.rfrs == "k.trimmed.min"){
			w[k.rfrs + 1] <- 1
		}
		else if (type.rfrs == "k.mean.min"){
			w[1 : k.rfrs] <- 1 / k.rfrs
		}
		else if (type.rfrs == "k.median.min"){
			if (k.rfrs %% 2 != 0){
				w[(k.rfrs + 1)/2] <- 1
			}
			else {
				w[k.rfrs/2] <- 1/2
			}
		}
		else if (type.rfrs == "k.trimmed.max"){
			w[length(data) - k.rfrs - 1] <- 1
		}
		else if (type.rfrs == "k.mean.max"){
			w[length(data) - k.rfrs: length(data)] <- 1/k.rfrs
		}
		else if (type.rfrs == "k.median.max"){
			if (k.rfrs %% 2 != 0){
				w[length(data) - (k.rfrs + 1)/2] <- 1
			}
			else {
				w[length(data) - k.rfrs/2] <- 1/2
			}
		}

		if (type.apprx == "lower"){
			w.lower <- w
			return(sum(new.data.sorted * w.lower))
		}
		else {
			w.upper <- rev(w)
			return(sum(new.data.sorted * w.upper))
		}
	}
	## It is used to IS.FRIS.FRST
	else if(type.method == "weight"){
		length.dt <- length(data)
		w <- matrix(0, nrow = length.dt)
		for (i in 1 : length.dt){
			w[i] <- 2 * (length.dt - i + 1)/(length.dt * (length.dt + 1))
		}
		
		return(sum(new.data.sorted * w))
	}
}

# It is used to calculate soft distance, but in this case we express it in relation (not distance)
#
# @title soft distance
# @param IND a vector representing indiscernibility relation
# @param penalty.fact a real number representing penalty factor
# @param type a type of approximations: "lower" or "upper"
calc.sd <- function(IND, penalty.fact = 0.06, type = "lower"){
	
	## convert into distance: d = 1 - R
	dist.xy <- (1 - IND)
	
	## get hard distance
	indx.hd <- which.min(dist.xy)
	sd <- matrix()
	
	## get soft distance
	for (i in 1 : length(dist.xy)){		
		if (i == indx.hd) {
			## hard distance
			sd[i] <- min(dist.xy, na.rm = TRUE)
		}
		else if (is.na(dist.xy[i])){
			## it will be ignored
			sd[i] <- NA
		}
		else{
			sd[i] <- dist.xy[i] - penalty.fact * length(which(dist.xy < dist.xy[i]))
		}		
	}
	
	## get results
	indx.sd <- which.max(sd)
	
	return(indx.sd)
}

# It is used to calculate lower and upper beta.pfrs
#
# @title calculate lower and upper appr. based on beta.pfrs
# @param imp.val.i a data frame resulted by indiscernibility relation
# @param tnorm.val.i a data frame resulted by indiscernibility relation
# @param beta.quasi a parameter beta precision quasi
calc.LU.betaPFRS <- function(imp.val, tnorm.val, beta.quasi){
	fuzzy.lower.Ra <- as.data.frame(matrix(nrow = 1, ncol = nrow(imp.val)))
	fuzzy.upper.Ra <- as.data.frame(matrix(nrow = 1, ncol = nrow(tnorm.val)))
	
	for (hh in 1 : nrow(imp.val)){
		n.max.lower = 0
		n.max.upper = 0
		## compute lower approximation	
		imp.val.desc <-  matrix(sort(imp.val[, hh,drop = FALSE], decreasing = TRUE), nrow = 1)

		k.t = 0
		right.side = 0
		for (ii in 1 : ncol(imp.val.desc)){
			right.side = right.side + (imp.val.desc[ii] * (1 - beta.quasi))
			if (k.t <= right.side){
				k.t = k.t + 1
			}
			else {
				n.max.lower = ii
				break
			}
		}
		if (n.max.lower == 0){
			n.max.lower = ncol(imp.val)
		}
		imp.val.real <- imp.val.desc[1:n.max.lower, drop = FALSE]	
		fuzzy.lower.Ra[1, hh] <- min(imp.val.real)	
		
		## compute upper approximation				
		tnorm.val.asce <- matrix(sort(tnorm.val[, hh,drop = FALSE], decreasing = FALSE), nrow = 1)
		k.s = 0
		right.side = 0
		for (ii in 1 : ncol(tnorm.val.asce)){
			right.side = right.side + ((1 - tnorm.val.asce[ii]) * (1 - beta.quasi))
			if (k.s <= right.side){
				k.s = k.s + 1
			}
			else {
				n.max.upper = ii
				break
			}
		}
		if (n.max.upper == 0)
			n.max.upper = ncol(tnorm.val)
			
		tnorm.val.real <- tnorm.val.asce[1, 1:n.max.upper, drop = FALSE]					
		fuzzy.upper.Ra[1, hh] <- max(tnorm.val.real) 
	}
	
	return(list(fuzzy.lower.Ra = fuzzy.lower.Ra, fuzzy.upper.Ra = fuzzy.upper.Ra))
}

## It is used to check type.aggregation parameter
# @param type.aggregation a type.aggregation parameter
ch.typeAggregation <- function(type.aggregation){
	if (type.aggregation[[1]] == "t.tnorm"){
		if (!is.na(list(type.aggregation[[2]]))){
			t.tnorm <- type.aggregation[[2]]
		}
		else t.tnorm <- "lukasiewicz"
		t.aggregation <- type.aggregation[[1]]
	}
	else if (type.aggregation[[1]] == "crisp"){
		t.tnorm <- "min"
		t.aggregation <- c("crisp")
	}
	else if (type.aggregation[[1]] == "custom"){
		t.tnorm <- type.aggregation[[2]]
		t.aggregation <- c("custom")
	}	
	return(c(t.aggregation = t.aggregation, t.tnorm = t.tnorm))
}

# It is used to Reduce indiscernibility relation 
# @param IND.relation a list indiscernibility relation of considered attributes
# @param t.tnorm a triangular norm operator
Reduce.IND <- function(IND.relation, t.tnorm, t.aggregation = "t.tnorm", delta = 0.5, variance.data = 0.5){	
	
	if (length(IND.relation) == 1){
		new.IND = IND.relation
	}
	else {
		if (t.aggregation == "custom"){
			new.IND <- t.tnorm(IND.relation)
		}
		else if (t.aggregation == "kernel.frst"){
			distance <- Reduce("+", IND.relation)/length(IND.relation)
			new.IND <- do.call(t.tnorm, list(distance, delta))
		}
		else {	
			temp.init <- IND.relation[[1]]
			if (class(t.tnorm) == "character"){ 				
				for (i in 2 : length(IND.relation)){
					## perform t.tnorm	
					temp.init <- do.call(func.tnorm, list(temp.init, IND.relation[[i]], t.tnorm))		
				}
			}
			else if (class(t.tnorm) == "function"){							
				for (i in 1 : nrow(IND.relation[[1]])){
					for (j in 1 : ncol(IND.relation[[1]])){
						temp.left <- IND.relation[[1]][i,j]
						for (k in 2 : length(IND.relation)){
							temp.left <- t.tnorm(temp.left, IND.relation[[k]][i, j])	
						}
						temp.init [i,j] <- temp.left
					}					
				}
			}
			else {
				stop("please enter the correct tnorm operator")
			}
			new.IND <- temp.init
		}
	}
	
	if (class(new.IND) != "matrix"){
		new.IND = new.IND[[1]]
	}
	
	return (new.IND)
}


# This is a function that builds the decision-relative discernibility matrix based on FRST.
# 
# @title The decision-relative discernibility matrix based on fuzzy rough sets
#
# @param decision.table a data frame representing a decision table. See \code{\link{BC.IND.relation.FRST}}. 
# @param type.discernibility a type of methoda to build discernibility matrix.
# @param num.red a type showing whether we want to produce all reducts or near-optimal reduction.
#         Currently, the near-optimal red is only for FVPRS.
# @param alpha.precision a numeric value representing variable precision of FVPRS. 
# @param type.relation a type of fuzzy relation. See \code{\link{BC.IND.relation.FRST}.
# @param t.tnorm a type of t-norm. See \code{\link{BC.IND.relation.FRST}.
# @param t.implicator a type of implicator operators. See \code{\link{BC.LU.approximation.FRST}}.
# @param type.LU a type of lower & upper approximations. See \code{\link{BC.LU.approximation.FRST}}.
# @param show.discernibilityMatrix a boolean values whether we want to show the matrix or not
# @param epsilon a numeric value for constructing matrix in "standard" method
# @param delta a value for gaussian equation.
build.discMatrix.FRST <- function(decision.table, type.discernibility = "fuzzy.disc", num.red = "all", 
                                        alpha.precision = 0.05, type.relation = c("tolerance", "eq.1"), 
		                               t.implicator = "lukasiewicz", 
									   type.LU = "implicator.tnorm", show.discernibilityMatrix = FALSE, 
									   epsilon = 0.001, delta = 2, type.aggregation = c("t.norm", "lukasiewicz")){
	## get the data
	objects <- decision.table
	num.att <- ncol(objects)
	desc.attrs <- attr(decision.table, "desc.attrs")
	nominal.att <- attr(decision.table, "nominal.attrs")
	if (is.null(attr(decision.table, "decision.attr"))){
		decision.attr <- ncol(objects)
	}
	else {
		decision.attr <- attr(decision.table, "decision.attr")
		## check if index of decision attribute is not in last index, then change their position
		if (decision.attr != ncol(objects)){
			objects <- cbind(objects[, -decision.attr, drop = FALSE], objects[, decision.attr, drop = FALSE])
			nominal.att <- c(nominal.att[-decision.attr], nominal.att[decision.attr])
		}
	}
	indx.univ <- c(seq(1, nrow(objects)))
	num.object <- nrow(objects)
	names.attr <- t(colnames(objects))
	
	###################################
	#### Perform fuzzy indiscernibility relation for all conditional attributes ####
	###################################
	attributes <- c(seq(1, (num.att - 1)))
	if (any(type.discernibility == c("fuzzy.disc", "standard.red", "consistence.degree", "alpha.red"))){
		## NOTE: disc.mat = TRUE means we produce discernibility matrix which is every element in matrix could have multiple values (attributes)
		## It should be noted that in this case, IND.cond = (1 - IND)
		#control.ind <- list(type.aggregation = c(t.aggregation, t.tnorm), type.relation = type.relation, disc.mat = TRUE)
		control.ind <- list(type.aggregation = type.aggregation, type.relation = type.relation, disc.mat = TRUE)
		IND.cond <- BC.IND.relation.FRST(decision.table, attributes = attributes, control = control.ind)$IND.relation 
	}
	else if (type.discernibility == "gaussian.red"){
		## NOTE: disc.mat = TRUE means we produce discernibility matrix which is every element in matrix could have multiple values (attributes)
		## In this case, we are using common gaussian "eq.2" for gaussian.kernel
		## It should be noted that in this case, IND.cond = (1 - IND)
		control.ind <- list(type.aggregation = c("t.tnorm", "t.cos"), type.relation = c("transitive.kernel", "gaussian"), disc.mat = TRUE)
		IND.cond <- BC.IND.relation.FRST(decision.table, attributes = attributes, control = control.ind)$IND.relation 
	}
	
	if (type.discernibility == "fuzzy.disc"){				
		##  Perform fuzzy indiscernibility relation for decision attributes ####
		#control.ind <- list(type.aggregation = c(t.aggregation, t.tnorm), type.relation = type.relation, disc.mat = FALSE)
		control.ind <- list(type.aggregation = type.aggregation, type.relation = type.relation, disc.mat = FALSE)
		IND.dec <- 1 - BC.IND.relation.FRST(decision.table, attributes = c(num.att), control = control.ind)$IND.relation 
		return (list(IND.conditionAttr = IND.cond, IND.decisionAttr = IND.dec))
	}	
	else if (any(type.discernibility == c("standard.red", "gaussian.red", "consistence.degree", "alpha.red"))){
#		req.suc <- require("sets", quietly=TRUE)
#		if(!req.suc) stop("In order to use this function, you need to install the package sets.")

		## Perform fuzzy indiscernibility relation for all conditional attributes ####
		attributes <- c(seq(1, (num.att - 1)))
		
		## standard.red is based on Eric C. C. Tsang et.al.: "Attributes Reduction Using Fuzzy Rough Sets"
		## consistence.degree is based on Eric C. C. Tsang and Zhao Suyun: "Decision table reduction in KDD: FR based Approach"
		if (any(type.discernibility == c("standard.red", "consistence.degree", "alpha.red"))){			
			control.ind <- list(type.aggregation = type.aggregation, type.relation = type.relation, disc.mat = FALSE)
			#control.ind <- list(type.aggregation = c(t.aggregation, t.tnorm), type.relation = type.relation, disc.mat = FALSE)
			IND.cond.all <- BC.IND.relation.FRST(decision.table, attributes = attributes, control = control.ind)
			IND.decAttr <- BC.IND.relation.FRST(decision.table, attributes = c(num.att), control = control.ind)

			## calculate lower approximation
			decision.attr = num.att

			#control.LU <- list(t.implicator = t.implicator, t.tnorm = t.tnorm, alpha = alpha.precision)
			control.LU <- list(t.implicator = t.implicator, t.tnorm = type.aggregation[2], alpha = alpha.precision)
			if (type.discernibility == "alpha.red"){
				FRST.all <- BC.LU.approximation.FRST(decision.table, IND.cond.all, IND.decAttr, 
                             type.LU = "fvprs", control = control.LU)	
			}
			else {
				FRST.all <- BC.LU.approximation.FRST(decision.table, IND.cond.all, IND.decAttr, 
                             type.LU = type.LU, control = control.LU)
			}
			fuzzy.lower <- FRST.all$fuzzy.lower	
			
			## get consistence degree (positive region) for "consistence.degree"
			if (type.discernibility == "consistence.degree"){
				cons.degree <- BC.positive.reg.FRST(decision.table, FRST.all)$positive.freg
			}
		}
			
		## gaussian.red is based on Chen et. al's technique "Parameterized attribute reduction with gaussian kernel based FRST"
		else if (type.discernibility == "gaussian.red"){
			## fuzzy similarity relation = gaussian; tnorm = T.cos
			control.ind <- list(type.aggregation = c("t.tnorm", "t.cos"), type.relation = c("transitive.kernel", "gaussian"), disc.mat = FALSE)
			
			## perform indiscernibility relation of all conditional attributes
			IND.cond.all <- BC.IND.relation.FRST(decision.table, attributes = attributes, control = control.ind)
			
			## perform indiscernibility relation of decision attributes
			IND.decAttr <- BC.IND.relation.FRST(decision.table, attributes = c(num.att), control = control.ind)
			
			## calculate lower approximation as "lambda" (based on the paper)
			## Note: we are using type.LU : gaussian.kernel with implicator cos (I.cos)
			control.LU <- list(type.relation = c("transitive.kernel", "gaussian")) 
			FRST.all <- BC.LU.approximation.FRST(decision.table, IND.cond.all, IND.decAttr, 
               type.LU = "gaussian.kernel", control = control.LU)
			fuzzy.lower <- FRST.all$fuzzy.lower	
		}
	
		## initialize the discernibility matrix
		disc.mat <- data.frame(matrix(NA, nrow = num.object, ncol = num.object))
		num.allAttr <- ncol(objects)
		names.fuzzy.lower <- names(fuzzy.lower)
		
		## create a function for checking
		check.attrDiscernibilityMat <- function(i, j){
			## construct the discernibility matrix by select discernible attributes
			disc.list <- c()
											
			## select different values on decision attribute only
			if (objects[i, num.allAttr] !=  objects[j, num.allAttr]) {				
				disc.attr <- formula.discernibilityMatrix(IND.cond, type.discernibility, fuzzy.lower, objects, names.attr, num.att, cons.degree, alpha.precision,
                                   				          num.allAttr, names.fuzzy.lower, epsilon = epsilon, i, j)				
				if (!is.null(disc.attr)){	
					
					if (show.discernibilityMatrix == TRUE){
						## construct one element has multi value/list which object is discerned.
						disc.mat[j, i] <<- data.frame(I(list(disc.attr)))
							
						## build a list for decision reduct
						disc.list <- append(disc.list, disc.attr)						
					}
					else {
						## build a list for decision reduct
						disc.list <- append(disc.list, disc.attr)

						disc.mat <<- NULL
					}
				}
			}				
			return(disc.list)
		}

		disc.list <- mapply(check.attrDiscernibilityMat, row(IND.cond), col(IND.cond))
		disc.list <- disc.list[!sapply(disc.list, is.null)]	
					
	}
	
	return(list(disc.mat = disc.mat, disc.list = disc.list))
}

## It is used to put formula constructing discernibility matrix for vectorization 
# @param IND.cond indiscernibility relation of all conditional attributes
# @param type.discernibility a type of discernibility matrix
# @param fuzzy.lower a list representing fuzzy lower approximation
# @param objects decision table
# @param names.attr a names of all attributes
# @param num.att the number of conditional attributes
# @param cons.degree consisten degree
# @param alpha.precision a numeric value between [0,1]
# @param num.allAttr the number of all attributes
# @param names.fuzzy.lower the names of decision concept
# @param epsilon a numeric values for constructing discernibility matrix
# @param i index of object
# @param j index of object
formula.discernibilityMatrix <- function(IND.cond, type.discernibility, fuzzy.lower, objects, names.attr, num.att, cons.degree = NULL, 
                                         alpha.precision = NULL, num.allAttr, names.fuzzy.lower, epsilon, i, j){
	disc.attr <- c()
	for (k in 1 : (num.att - 1)){
		## select attributes which discerned based on standard.red
		## E. C. C. Tsang, D. G. Chen, D. S. Yeung, X. Z. Wang, and J. W. T. Lee, 
		## "Attributes Reduction Using Fuzzy Rough Sets"
		if (type.discernibility == "standard.red"){
			## check every values in each element of matrix
			if (fuzzy.lower[[which(names.fuzzy.lower == objects[i, num.allAttr])]][j] < fuzzy.lower[[which(names.fuzzy.lower == objects[i, num.allAttr])]][i]){
				## It should be noted that in this case, IND.cond = (1 - IND)
				if ((unlist(IND.cond[i, j])[k]) >= fuzzy.lower[[which(names.fuzzy.lower == objects[i, num.allAttr])]][i]){
					disc.attr <- c(disc.attr, names.attr[k])
				}
			}
		}
		## select attributes which discerned based on gaussian.red
		## D. G. Chen, Q. H. Hu, and Y. P. Yang, "Parameterized Attribute Reduction with
		## Gaussian Kernel Based Fuzzy Rough Sets"
		else if (type.discernibility == "gaussian.red"){
			## R(x,y) means 1 - IND.cond
			## lambda means fuzzy.lower - epsilon
			if ((1 - (unlist(IND.cond[i, j])[k])) <= c(sqrt(1 - (fuzzy.lower[[which(names.fuzzy.lower == objects[i, num.allAttr])]][i] - epsilon)^2))){
				disc.attr <- c(disc.attr, names.attr[k])
			}
		}
		## select attributes which discerned based on consistence degree: 
		## E. C. C. Tsang and S. Y. Zhao, "Decision Table Reduction in KDD: Fuzzy Rough Based Approach"
		else if (type.discernibility == "consistence.degree"){
			## the formula: T(a(x, y), lambda)
			## where a(x, y) is indiscernibility relation, that means (1 - IND.cond) 
			## lambda is consistence degree, that means positive region
			## based on the paper, we fix the t-norm which is "lukasiewicz" 
			tnorm.val <- func.tnorm((1 - (unlist(IND.cond[i, j])[k])), cons.degree[i], t.tnorm = "lukasiewicz")
			if (tnorm.val == 0){
				disc.attr <- c(disc.attr, names.attr[k])
			}
		}
		## select attributes which discerned based on alpha reduction of FVPRS: 
		## S. Zhao, E. C. C. Tsang, and D. Chen, "The Model of Fuzzy Variable Precision Rough Sets"
		else if (type.discernibility == "alpha.red"){
			tnorm.val <- func.tnorm((1 - (unlist(IND.cond[i, j])[k])), (fuzzy.lower[[which(names.fuzzy.lower == objects[i, num.allAttr])]][i]), t.tnorm = "lukasiewicz")
			if (tnorm.val <= alpha.precision){
				disc.attr <- c(disc.attr, names.attr[k])
			}
		}
	}
	return(disc.attr)
}

# It is used to calculate minimal element 
# @param decision.table a list of decision table
# @param t.tnorm a triangular norm
# @param type.relation a type of relation
# @param t.implicator a type of implicator operator
# @param type.LU a type of lower/upper approximation
min.disc.mat.FRST <- function(decision.table, t.tnorm = "lukasiewicz", type.relation = c("tolerance", "eq.1"), t.implicator = "lukasiewicz", type.LU = "implicator.tnorm"){
#	req.suc <- require("sets", quietly=TRUE)
#	if(!req.suc) stop("In order to use this function, you need to install the package sets.")	  
		
	## get data
	objects <- decision.table
	desc.attrs <- attr(decision.table, "desc.attrs")
	nominal.att <- attr(decision.table, "nominal.attrs")
	if (is.null(attr(decision.table, "decision.attr"))){
		decision.attr <- ncol(objects)
	}
	else {
		decision.attr <- attr(decision.table, "decision.attr")
		## check if index of decision attribute is not in last index, then change their position
		if (decision.attr != ncol(objects)){
			objects <- cbind(objects[, -decision.attr, drop = FALSE], objects[, decision.attr, drop = FALSE])
			nominal.att <- c(nominal.att[-decision.attr], nominal.att[decision.attr])
		}
	}

	num.object <- nrow(objects)
	temp <- matrix()
	num.att <- ncol(objects)
	t.similarity <- type.relation[2]
	delta <- 0.2
	
	## initialization
	miu.Ra <- matrix(nrow = num.object, ncol = num.object)
	miu.Ra.df <- data.frame()
	
	attributes <- c(seq(1, (num.att - 1)))
    control.ind <- list(type.aggregation = c("t.tnorm", t.tnorm), type.relation = type.relation)
	
    #### Perform fuzzy indiscernibility relation for all attributes ####
	IND <- BC.IND.relation.FRST(decision.table, attributes = attributes, control = control.ind)
	IND.decAttr <- BC.IND.relation.FRST(decision.table, attributes = c(num.att), control = control.ind)
	
	#### Perform fuzzy lower and upper approximation using type.LU : "implicator.tnorm" ####	 
	decision.attr = ncol(objects)
	control <- list(t.implicator = t.implicator ,t.tnorm = t.tnorm)
	
	res.lower <- BC.LU.approximation.FRST(decision.table, IND, IND.decAttr, 
               type.LU = type.LU, control = control)
	fuzzy.lower <- res.lower$fuzzy.lower	

	positive.region <- matrix(BC.positive.reg.FRST(decision.table, res.lower)$positive.freg)
	
	att.dis.r <- list()
	dis.r <- matrix(nrow = 1, ncol = 2)	
	obj.decisionAttr <- objects[, ncol(objects), drop = FALSE]
	
	## Create function to be vectorized
	formula.min.discernibilityMatrix <- function(i, j, obj, obj.decisionAttr, 
                                           nominal.att, t.similarity, delta, range.data, variance.data, fuzzy.lower, 
											h){
		## loop for each objects
		temp <- obj[i]
			if (obj.decisionAttr[i, 1] !=  obj.decisionAttr[j, 1]) {
			## calculate the similarity between two objects based on type of similarity				
			if (nominal.att[h] == TRUE){
				## overwrite type of similarity when crisp
				t.similarity.new <- "boolean"
			}	else {
			  t.similarity.new <- t.similarity
			}		
			temp.miu.Ra <- unlist(similarity.equation(temp, obj[j], t.similarity.new, delta, range.data, variance.data))

			## check discernibility among objects					
			if ((1 - temp.miu.Ra) >= (fuzzy.lower[[which(names(fuzzy.lower) == obj.decisionAttr[i, 1])]][i])){	
				## collect them
				if (is.na(dis.r[1, 1])){
					 ## get a pair of objects
					 ## update the outer variable !!!!
					 dis.r[1, ] <<- matrix(c(i, j), nrow = 1)
					 
					 ## get its attribute
					 ## update the outer variable !!!!
					 att.dis.r <<- list(h)
				}
				else {
					## perform this for avoiding duplication
					match.pt <- which(apply(dis.r, 1, function(x) all(x == c(i,j))))
					
					## save index representing a name of attribute R (if they are the same)
					if (length(match.pt) == 1){			
						## update the outer variable !!!!
						att.dis.r[[match.pt]] <<- append(att.dis.r[[match.pt]], h)
					}
					else {
						## save objects (x, y) and its index representing a name of attribute R (on new rows)
						## dis.r is used to save pairs of objects
						## att.dis.r is used to save the attributes of each pair of objects
						## update the outer variable !!!!
						dis.r <<- rbind(dis.r, c(i, j))
						att.dis.r <<- append(att.dis.r, h)
					}
				}					
			}					
		}
		else {
			temp.miu.Ra <- NA
		}		
		return(temp.miu.Ra)
	}											 

	## vectorize the function
	VecFun <- Vectorize(formula.min.discernibilityMatrix, vectorize.args=list("i","j"))
											   
	for (h in 1 : (num.att - 1)){
		## calculate variance and range of data
		res.varRange <- cal.var.range(objects, c(h), nominal.att)
		variance.data <- res.varRange$variance.data
		range.data <- res.varRange$range.data
	
		outer(1 : num.object, 1 : num.object, VecFun, objects[, h], obj.decisionAttr, 
                                           nominal.att, t.similarity, delta, range.data, variance.data, fuzzy.lower, 
                                           h)

	}
		
	## create matrix to save the number of attributes (N)
	dis.N.r <- matrix(nrow = nrow(dis.r), ncol = 4)
	
	## get N as the number of attribute
	for (ii in 1 : nrow(dis.r)){
		num <- length(att.dis.r[[ii]])
		## we construct a matrix where 
		## 1st column is sequence index
		## 2nd and 3rd column are a pair of objects
		## 4th column is the number of the attribute
		dis.N.r[ii, ] <- cbind(ii, dis.r[ii, ,drop = FALSE], num)
	}
	
	## sort the matrix based on the number of attribute (4th column)
	dis.N.r <- dis.N.r[order(dis.N.r[, 4], dis.N.r[, 2], dis.N.r[, 3]), ]
	## initialization
	sseq <- seq(1, nrow(dis.N.r))
	red.att <- list()
	exit <- FALSE
	## iterate until stopping criteria satified
	while (nrow(dis.N.r > 0) && exit == FALSE){
		## get reduct
		red <- att.dis.r[[dis.N.r[1, 1]]]	
		## detect the same attribute to be deleted
		del.pt <- lapply(att.dis.r, function(x) x %in% red)
		
		new.dis.N.r <- matrix(nrow = 1, ncol = ncol(dis.N.r))
		new.att.dis.r <- list()		
		for (kk in 1 : length(del.pt)){
			## if the attribute is not the same
			if (all(del.pt[[kk]] == FALSE)){
				temp <- dis.N.r[which(dis.N.r[, 1] == kk), ,drop = FALSE]				
				temp.l <- att.dis.r[[kk]]
				## update/collect them
				if (is.na(new.dis.N.r[1,1])){
					new.dis.N.r <- temp
					new.att.dis.r <- list(temp.l)
				}
				else {
					new.dis.N.r <- rbind(new.dis.N.r, temp)
					new.att.dis.r <- append(new.att.dis.r, list(temp.l))
				}
			}
		}
		## update matrix of objects and list of attributes
		dis.N.r <- new.dis.N.r		
		att.dis.r <- new.att.dis.r	
		
		## update reducts
		if (nrow(dis.N.r) == 1 && !is.na(dis.N.r[1,1])){
			red.att <- c(red, list(att.dis.r[[1]]))
			exit <- TRUE
		}
		else if (is.na(dis.N.r[1,1])){
			red.att <- c(red.att, list(red))
			exit <- TRUE
		}
		else {
			up.seq <- seq(1, nrow(dis.N.r))
			dis.N.r[, 1] <- up.seq
			red.att <- c(red.att, list(red))
		}
	}
	
	## change indexes into names of attributes
	disc.list <- names(objects[red.att[[1]]])
	for (i in 2 : length(red.att)){
		disc.list <- list(disc.list, names(objects[red.att[[i]]]))
	}
  
	discernibilityMatrix = list(disc.mat = NULL, disc.list = disc.list, type.discernibility = "min.element", type.model = "FRST")
	discernibilityMatrix = ObjectFactory(discernibilityMatrix, classname = "DiscernibilityMatrix")
	return(discernibilityMatrix)

#	return(calc.discernibility.func(disc.mat = NULL, disc.list = disc.list, type.discernibility = "min.element", type.model = "FRST"))
}

# an auxiliary function for dividing a discernibility class into subclasses corresponding to values of a given nominal vector
# @param INDclass equivalence class
# @param splitVec a vector of splitted INDclass
splitINDclass <- function(INDclass, splitVec)  {
  split(INDclass, splitVec[INDclass], drop = TRUE)
}

# It is used to calculate transitive closure
#
# @param IND.relation a indiscernibility relation
calc.transitive.closure <- function(IND.relation, type.MVCompletion){

	new.IND.relation <- IND.relation
	exit <- FALSE
	
	## vectorize the function
	func <- function(i, j, IND.relation){
		x.z <- IND.relation[i]
		z.y <- IND.relation[j]
		new.IND <- max(IND.relation[i, j], max(pmin(x.z, z.y)))
	}	
	vec.func <- Vectorize(func, vectorize.args=list("i","j"))

	while (exit == FALSE){
		## do calculation
		if (type.MVCompletion == TRUE){
			stop("The decision table contains missing values. For current version, they have not supported yet. 
			     Please perform missing value preprocessing.")
		}
		else {
			new.IND.relation <- outer(1 : nrow(IND.relation), 1 : ncol(IND.relation), vec.func, IND.relation)
			
			## stopping criteria
			if (isTRUE(all.equal(new.IND.relation, IND.relation))){
				exit <- TRUE
				return(new.IND.relation)
			}
			else {
				IND.relation <- new.IND.relation
			}
		}
	}
}
## It is used to calculate kernilized FRST: gaussian
func.gaussian.kernel <- function(distance, delta){
	## calculate for all column/attributes
	Rg <- exp(- (distance^2)/delta)								

	return (Rg)
}

## It is used to calculate kernilized FRST: exponential
func.exponential.kernel <- function(distance, delta){
	## calculate for all column/attributes
	Re <- exp(- distance/delta)	
	
	return (Re)
}

## It is used to calculate kernilized FRST: rational quadratic kernel
func.rational.kernel <- function(distance, delta){
	## calculate for all column/attributes
	Rr <- 1 - (distance^2/ (distance^2 + delta))	
	
	return (Rr)
}

## It is used to calculate kernilized FRST: circular kernel
func.circular.kernel <- function(distance, delta){
	if (distance >= delta){
		warnings("the condition is not satisfied, please increase the value of delta")
	}
	## calculate for all column/attributes
	Rc <- 2/pi * acos(distance/delta) - 2/pi * 
						(distance/delta) * sqrt(1 - (distance/delta)^2)	
	
	return (Rc)
}

## It is used to calculate kernilized FRST: spherical kernel
func.spherical.kernel <- function(distance, delta){
	if (distance >= delta){
		warnings("the condition is not satisfied, please increase the value of delta")
	}
	## calculate for all column/attributes
	Rs <- 1 - 3/2 *  distance / delta + 0.5 * 
						(distance / delta)^3	
	
	return (Rs)
}		
