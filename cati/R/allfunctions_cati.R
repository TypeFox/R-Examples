# traits is the Matrix of traits and factors are the nested factors to take into account in the partition of variance

partvar <- function(traits, factors, printprogress = TRUE){

	message("The partvar function decomposes the variance accross nested scales. Thus choose the order of the factors very carefully!")

	traits <- as.matrix(traits)
	factors <- as.matrix(factors)
	nfactors <- ncol(factors)
	ntraits <- ncol(traits)
	namestraits <- colnames(traits)
	res <- matrix(0, nrow = nfactors+1, ncol = ntraits)
	colnames(res) <- colnames(traits)

	if (!is.null(colnames(factors)))
		{rownames(res) <- c(colnames(factors), "within")
	}

	else {
		rownames(res) <- c(paste("factor",1:(nfactors),sep = ""), "within")
		colnames(factors) <- c(paste("factor", 1:(nfactors),sep = ""))
	}

	factors <- as.data.frame(factors)

	for (t in 1 : ntraits) {
		if(sum(is.na(traits[,t])) > 0){
			trait <- as.vector(na.omit(traits[,t]))
			warning(paste("All individuals with one NA (", sum(is.na(traits[,t])),"individual value(s)) are excluded for the trait", t, ":", namestraits[t], sep = " "))
			fact <-  as.data.frame(factors[!is.na(traits[,t]),])
		}
		else{
			trait <- traits[,t]
			fact <- as.data.frame(factors)
		}

		functionlme = paste('varcomp(lme(trait~1, random = ~1|', paste(colnames(factors), collapse = '/'), ",na.action = na.omit),1)", sep = "")
		res[,t] <- as.vector(eval(parse(text = functionlme), envir = fact))

		if (printprogress == TRUE)
			{print(paste(round(t/ntraits*100, 2), "%", sep = " ")) }
		else{}
	}

	class(res) <- "partvar"

	return(res)
}

#Pie of variance partitioning
piePartvar <- function(partvar, col.pie = NA, ...){
	nfactors <- nrow(partvar)
	ntraits <- ncol(partvar)

	if(any(is.na(col.pie))){
		col.pie <- funky.col(nfactors)
	}

	for (t in 1 : ntraits) {
		pie(partvar[,t], main = colnames(partvar)[t], col = col.pie , labels = rownames(partvar), ...)
	}
}

#barplot of variance partitioning
barPartvar <- function(partvar, col.bar = NA, ...){

	nfactors <- nrow(partvar)

	oldpar <- par(no.readonly = TRUE)
	par(mar = c(5,6.5,4,2), cex = 0.7)

	if(any(is.na(col.bar))){
		col.bar <- funky.col(nfactors)
	}

	barplot(partvar, col = col.bar, las = 1, horiz = T, xlab = "% of variance", ...)

  par(oldpar)
}


#______________#______________#______________#______________#______________#______________#______________#______________
#______________#______________#______________#______________#______________#______________#______________#______________
#__Tstats

### Function to calculation Tstats
Tstats <- function(traits, ind.plot, sp, SE = 0, reg.pool = NULL, SE.reg.pool = NULL, nperm = 99, printprogress = TRUE){
	#6 variances: I: individual, P: population, C: community, R: region
	#IP; IC; IR; PC; PR; CR

	#traits is the matrix of individual traits, ind.plot is the name of the plot in which the individual is (factor type), and sp is the species name of each individual

	if (sum(is.na(traits))>0) {
		warning(paste("This function exclude", sum(is.na(traits)),"Na values", sep=" "))
	}

	if (!is.matrix(traits)) {
		traits <- as.matrix(traits)
	}

	names_sp_ind.plot <- as.factor(paste(sp, ind.plot, sep = "@"))
	Tplosp <- unlist(strsplit(levels(names_sp_ind.plot), split = "@"))[2*(1:nlevels(names_sp_ind.plot))]
	names(Tplosp) <- levels(names_sp_ind.plot)
	#Tplosp is the plot in wich the population is

	if(!is.null(nperm)){
		if (nperm == 0) {nperm = NULL}
	}

	if(length(SE) == 1){
		SE <- rep(SE, times = ncol(traits))
	}

	########################################
	####	calculation of observed values	####
	########################################

	#________________________________________
	#Objects creation

	mean_IP <- matrix(nrow = nlevels(names_sp_ind.plot), ncol = ncol(traits))
	rownames(mean_IP) = levels(names_sp_ind.plot)
	mean_PC <- matrix(nrow = nlevels(ind.plot), ncol = ncol(traits))

	var_IP <- matrix(nrow = nlevels(names_sp_ind.plot), ncol = ncol(traits))
	var_PC <- matrix(nrow = nlevels(ind.plot), ncol = ncol(traits))

	var_CR <- vector()
	var_IC <- matrix(nrow = nlevels(ind.plot), ncol = ncol(traits))
	var_PR <- vector()
	var_IR <- vector()

	T_IP.IC <- matrix(nrow = nlevels(ind.plot), ncol = ncol(traits))
	T_IC.IR <- matrix(nrow = nlevels(ind.plot), ncol = ncol(traits))
	T_PC.PR <- matrix(nrow = nlevels(ind.plot), ncol = ncol(traits))


	for (t in 1: ncol(traits)){
		mean_IP[,t] <- tapply(traits[,t], names_sp_ind.plot ,mean, na.rm = T)
		mean_PC[,t] <- tapply(mean_IP[,t], Tplosp , mean, na.rm = T)

		var_IP[,t] <- tapply(traits[,t], names_sp_ind.plot, var, na.rm = T)
		var_PC[,t] <- tapply(mean_IP[,t], Tplosp ,var, na.rm = T)
		var_CR[t] <- var(mean_PC[,t], na.rm = T)
		var_IC[,t] <- tapply(traits[,t], ind.plot ,var, na.rm = T)
		var_PR[t] <- var(as.vector(mean_IP[,t]), na.rm = T)
		var_IR[t] <- var(traits[,t], na.rm = T)

		for(s in 1 : nlevels(ind.plot)){
			T_IP.IC[s,t] <- mean(var_IP[grepl(levels(ind.plot)[s],Tplosp),t], na.rm = T)/var_IC[s,t]
			T_IC.IR[s,t] <- var_IC[s,t]/var_IR[t]
			T_PC.PR[s,t] <- var_PC[s,t]/var_PR[t]
		}
	}

	#________________________________________

	#########################################
	#### 	 Creating null models 	 ####
	#########################################

	#null model 1 = local
	#null model 2 = regional.ind
	#null model 2.sp = regional.pop

	if (is.numeric(nperm)){

		var_IP_nm1 <- array(dim = c(nperm,ncol(traits),nrow = length(Tplosp)))
		var_PC_nm2sp <- array(dim = c(nperm,ncol(traits),nlevels(ind.plot)))
		var_IC_nm1 <- array(dim = c(nperm,ncol(traits),nlevels(ind.plot)))
		var_IC_nm2 <- array(dim = c(nperm,ncol(traits),nlevels(ind.plot)))
		var_PR_nm2sp <- array(dim = c(nperm,ncol(traits)))
		var_IR_nm2 <- array(dim = c(nperm,ncol(traits)))

		mean_IP_nm2sp <- array(dim = c(nperm,ncol(traits),length(Tplosp)))
		mean_PC_nm2sp <- array(dim = c(nperm,ncol(traits),nlevels(ind.plot)))

		traits.nm1 <- list()
		traits.nm2 <- list()
		traits.nm2sp <- list()

		T_IP.IC_nm1 <- array(dim = c(nperm,ncol(traits),nlevels(ind.plot)))
		T_IC.IR_nm2 <- array(dim = c(nperm,ncol(traits),nlevels(ind.plot)))
		T_PC.PR_nm2sp <- array(dim = c(nperm,ncol(traits),nlevels(ind.plot)))


		#########################################
		#### 	 Creating regional pools 	 ####
		#########################################

		if(!is.null(reg.pool)){isnullregpool <- FALSE} else{isnullregpool <- TRUE}

		# If the regional pool is the same for all communitiers:
		# creating a list of regional pool (one regional pool by community)
		if (is.null(reg.pool)) {
			reg.pool <- rep(list(traits), nlevels(ind.plot))
		}

		if (is.data.frame(reg.pool) | is.matrix(reg.pool) ) {
			reg.pool <- rep(list(reg.pool), nlevels(ind.plot))
		}

		# Standard error for regional pool
		if(is.null(SE.reg.pool)){
			SE.reg.pool <- SE
		}

		if(length(SE.reg.pool) == 1){
			SE.reg.pool <- rep(SE.reg.pool, times = ncol(traits))
		}

		if (is.vector(SE.reg.pool)) {
			SE.reg.pool <- rep(list(SE.reg.pool), nlevels(ind.plot))
		}

		if(length(reg.pool) != nlevels(ind.plot)){
			stop("reg.pool need to be either a matrix or a list of length equal to the number of communities")
		}

		#________________________________________
		# Warnings and stop about standard errors parameter
		if (length(SE) != ncol(traits) & length(SE) != 1) {
			stop("The vector SE need to have a length of one or equal to the number of traits")
		}

		if(!isnullregpool){
			if (length(SE.reg.pool) != length(reg.pool)) {
				stop("The vector SE.reg.pool need to have the same dimension as reg.pool")
			}
		}

		#########################################
		####  End creating regional pools  ####
		#########################################

		# Creation of three null models
		if (printprogress == T){print("creating null models")}


		#________________________________________
		# Null model 1: Sample individual traits values within communities
		for (t in 1: ncol(traits)){
			traits.nm1[[t]] <- list()
			for(s in 1: nlevels(ind.plot)) {
				traits.nm1[[t]][[s]] <- list()
				for(i in 1:nperm){

					# Take measurement standard error into account
					trait.intern <- traits[,t]
					if (SE[t] != 0) {
						trait.intern <- rnorm( length(trait.intern), mean = trait.intern, sd = SE[t])
					}

					# Sample individual traits values within communities
					if (length(traits[ind.plot == levels(ind.plot)[s], t]) != 1) {
						perm_ind.plot1 <- sample(trait.intern[ind.plot == levels(ind.plot)[s]], table(ind.plot)[s])
						traits.nm1[[t]][[s]][[i]] <- perm_ind.plot1
					}
					else {traits.nm1[[t]][[s]][[i]] <- "NA"}
				}
			}
			if (printprogress == T){print(paste(round(t/ncol(traits)/3*100,2),"%")) } else {}
		}

		#________________________________________
		# Null model 2: Sample individual traits values in the region
		for (t in 1: ncol(traits)){
			traits.nm2[[t]] <- list()
			for(s in 1: nlevels(ind.plot)) {
				traits.nm2[[t]][[s]] <- list()
				for(i in 1:nperm){

					# Take measurement standard error into account
					trait.intern <- reg.pool[[s]][, t]
					SE.reg.pool.intern <- SE.reg.pool[[s]][t]
					if (SE.reg.pool.intern != 0) {
						trait.intern <- rnorm(length(trait.intern), mean = trait.intern, sd = SE.reg.pool.intern)
					}

					# Sample individual traits values in the region
					perm_ind.plot2 <- sample(trait.intern, table(ind.plot)[s])
					traits.nm2[[t]][[s]][[i]] <- perm_ind.plot2
				}
			}
			if (printprogress == T){print(paste(round(33.3+t/ncol(traits)/3*100, 2),"%"))} else {}
		}

		#________________________________________
		# Null model 2sp: Sample populationnal traits values in the region

		traits_by_sp <- apply(traits, 2, function(x) tapply(x, names_sp_ind.plot, mean))
		traits_by_pop <- traits_by_sp[match(names_sp_ind.plot, rownames(traits_by_sp)), ]
		#traits_by_sp <- aggregate(traits, by = list(names_sp_ind.plot), mean, na.rm = T)[,-1]

		if (!is.matrix(traits_by_pop)) {
			traits_by_pop <- as.matrix(traits_by_pop)
		}

		for (t in 1: ncol(traits)){
			traits.nm2sp[[t]] <- list()
			for(s in 1: nlevels(ind.plot)){
				traits.nm2sp[[t]][[s]] <- list()
				for(i in 1:nperm){

					########################## Not trivial to use measurement error at the individual level for calculation of populationnal mean

					# Take measurement standard error into account ???????????
					#trait.intern <- traits_by_pop[,t]
					#if (SE[t] != 0) {
					#	trait.intern <- rnorm(length(trait.intern), mean = trait.intern, sd = SE[t])
					#}

					# Sample populationnal traits values in the region
					perm_ind.plot2sp <- sample(traits_by_pop[,t], table(ind.plot)[s])
					traits.nm2sp[[t]][[s]][[i]] <- perm_ind.plot2sp
				}
			}
			if (printprogress == T){print(paste(round(66.6+t/ncol(traits)/3*100, 2),"%"))} else {}
		}

		#________________________________________

		#########################################
		#### calculation of Tstats on null models ####
		#########################################

		if (printprogress == T){print("calculation of Tstats using null models")}

		yy <- length(names_sp_ind.plot)
		for (t in 1: ncol(traits)){
			for(i in 1:nperm){
				mean_IP_nm2sp[i,t, ] <- tapply(unlist(traits.nm2sp[[t]])[(1+(i-1)*yy) : (i*yy)], names_sp_ind.plot ,function(x) mean(x, na.rm = T))
				mean_PC_nm2sp[i,t, ] <- tapply(mean_IP_nm2sp[i,t, ], Tplosp, mean, na.rm = T)
			}
			if (printprogress == T){print(paste(round(t/ncol(traits)/3*100, 2),"%"))} else {}
		}


		for (t in 1: ncol(traits)){
			for(i in 1:nperm){
				var_IP_nm1[i,t, ] <- tapply(unlist(traits.nm1[[t]])[(1+(i-1)*yy) : (i*yy)], names_sp_ind.plot ,function(x) var(x, na.rm = T))
				var_PC_nm2sp[i,t, ] <- tapply(mean_IP_nm2sp[i,t, ], Tplosp ,var, na.rm = T)
				var_IC_nm1[i,t, ] <- tapply(unlist(traits.nm1[[t]])[(1+(i-1)*yy) : (i*yy)], ind.plot ,function(x) var(x, na.rm = T))
				var_IC_nm2[i,t, ] <- tapply(unlist(traits.nm2[[t]])[(1+(i-1)*yy) : (i*yy)], ind.plot ,function(x) var(x, na.rm = T))
				var_PR_nm2sp[i,t] <- var(as.vector(mean_IP_nm2sp[i,t, ]), na.rm = T)
				var_IR_nm2[i,t] <- var(unlist(traits.nm2[[t]])[(1+(i-1)*yy) : (i*yy)], na.rm = T)
			}
			if (printprogress == T){print(paste(round(33.3+t/ncol(traits)/3*100, 2),"%"))} else {}
		}


		for (t in 1: ncol(traits)){
			for(i in 1:nperm){
				for(s in 1 : nlevels(ind.plot)){
					T_IP.IC_nm1[i,t,s] <- mean(var_IP_nm1[i,t,grepl(levels(ind.plot)[s],Tplosp)], na.rm = T)/var_IC_nm1[i,t,s]
					T_IC.IR_nm2[i,t,s] <- var_IC_nm2[i,t,s]/var_IR_nm2[i,t]
					T_PC.PR_nm2sp[i,t,s] <- var_PC_nm2sp[i,t,s]/var_PR_nm2sp[i,t]
				}
			}
			if (printprogress == T){print(paste(round(66.6+t/ncol(traits)/3*100, 2),"%"))} else {}
		}

	}#end of calculation of Tstats using null models

	colnames(T_IP.IC) <- colnames(traits)
  colnames(T_IC.IR) <- colnames(traits)
  colnames(T_PC.PR) <- colnames(traits)

	if (is.numeric(nperm)){
		colnames(T_IP.IC_nm1) <- colnames(traits)
		colnames(T_IC.IR_nm2) <- colnames(traits)
		colnames(T_PC.PR_nm2sp) <- colnames(traits)
	}

	rownames(T_IP.IC) <- levels(as.factor(Tplosp))
 	rownames(T_IC.IR) <- levels(as.factor(Tplosp))
 	rownames(T_PC.PR) <- levels(as.factor(Tplosp))


	#________________________________________
  res <- list()

  res$Tstats <- list()

  res$Tstats$T_IP.IC <- T_IP.IC
  res$Tstats$T_IC.IR <- T_IC.IR
  res$Tstats$T_PC.PR <- T_PC.PR

  res$variances <- list()

  res$variances$var_IP <- var_IP
  res$variances$var_PC <- var_PC
  res$variances$var_CR <- var_CR
	res$variances$var_IC <- var_IC
  res$variances$var_PR <- var_PR
  res$variances$var_IR <- var_IR

	if (is.numeric(nperm)){
		res$variances$var_IP_nm1 <- var_IP_nm1
		res$variances$var_PC_nm2sp <- var_PC_nm2sp
		res$variances$var_IC_nm1 <- var_IC_nm1
		res$variances$var_IC_nm2 <- var_IC_nm2
		res$variances$var_PR_nm2sp <- var_PR_nm2sp
		res$variances$var_IR_nm2 <- var_IR_nm2

		res$Tstats$T_IP.IC_nm <- T_IP.IC_nm1
    	res$Tstats$T_IC.IR_nm <- T_IC.IR_nm2
    res$Tstats$T_PC.PR_nm <- T_PC.PR_nm2sp
  }
  else{}

 	res$traits <- traits
 	res$ind.plot <- ind.plot
 	res$sp <- sp
 	res$sites_richness <- table(ind.plot)
 	res$namestraits <- colnames(traits)
 	res$call <- match.call()

	class(res) <- "Tstats"
	invisible(res)
}

print.Tstats <- function(x, ...){

	if (!inherits(x, "Tstats"))
		{stop("x must be a list of objects of class Tstats")
	}

	cat("\t##################\n")
	cat("\t# T-statistics #\n")
	cat("\t##################\n")
	cat("class: ")
	cat(class(x))

	cat("\n$call: ")
	print(x$call)

	cat("\n###############")
	cat("\n$Tstats: list of observed and null T-statistics")
	cat("\n\nObserved values")

	cat("\n\t$T_IP.IC: ratio of within-population variance to total within-community variance")
	cat("\n\t$T_IC.IR: community-wide variance relative to the total variance in the regional pool")
	cat("\n\t$T_PC.PR: inter-community variance relative to the total variance in the regional pool\n")

	if(length(x$Tstats) == 6){
	  cat("\nNull values, number of permutations:", dim(x$Tstats$T_IP.IC_nm)[1])
	  cat("\n\t$T_IP.IC_nm: distribution of T_IP.IC value under the null model local")
	  cat("\n\t$T_IC.IR_nm: distribution of T_IC.IR value under the null model regional.ind ")
	  cat("\n\t$T_PC.PR_nm: distribution of T_PC.PR value under the null model regional.pop\n")
	}

	cat("\n###############")
	interm <- ""
	if(length(x$Tstats) == 6){interm <- "and null"}
	cat("\n$variances: list of observed", interm, "variances")

	TABDIM <- 3

	sumry <- array("", c(TABDIM, 4), list(1:TABDIM, c("data", "class", "dim", "content")))
	sumry[1, ] <- c("$traits", class(x$traits), paste(dim(x$traits)[1], dim(x$traits)[2], sep = ",") , "traits data")
	sumry[2, ] <- c("$ind.plot", class(x$ind.plot), length(x$ind.plot), "name of the plot in which the individual is")
	sumry[3, ] <- c("$sp", class(x$sp), length(x$sp), "individuals' groups (e.g. species)")

	class(sumry) <- "table"
	cat("\n\n###############")
	cat("\ndata used\n")
	print(sumry)

	cat("\n###############\n")
	cat("others")
	if(length(x$namestraits)>0) {
		cat("\n\t$namestraits:", length(x$namestraits), "traits\n")
		print(head(x$namestraits))
		if(length(x$namestraits)>5){print("...")}
	}
	else { cat("\n\t$namestraits:", dim(x$traits)[1], "traits\n") }

	rich <- x$sites_richness
	class(rich) <- "table"
	cat("\n\t$sites_richness:\n\t")
	print(rich)
	cat("\n")
}

print.ComIndex <- function(x, ...){

	if (!inherits(x, "ComIndex"))
		{stop("x must be a list of objects of class ComIndex")
	}

	cat("\t#################################\n")
	cat("\t# Community metrics calculation #\n")
	cat("\t#################################\n")
	cat("class: ")
	cat(class(x))

	cat("\n$call: ")
	print(x$call)

	cat("\n###############")
	cat("\n$obs: list of observed values\n")
	seq.interm <- seq(1,length(names(x$list.index)), by = 2)
	cat(paste("\t$", names(x$list.index)[seq.interm], sep = ""), sep = "\n")

	if(!is.null(x$null)){
	  cat("\n###############")
	  cat("\n$null: list of null values, number of permutations:", dim(x$null[[1]])[3], "\n")
	cat(paste(paste("\t$", names(x$list.index)[-seq.interm], sep = ""), paste("null model", "=", x$nullmodels, sep=" "), sep=" ... "), sep = "\n")
	}

	TABDIM <- 3

	sumry <- array("", c(TABDIM, 4), list(1:TABDIM, c("data", "class", "dim", "content")))
	sumry[1, ] <- c("$traits", class(x$traits), paste(dim(x$traits)[1], dim(x$traits)[2], sep = ",") , "traits data")
	sumry[2, ] <- c("$ind.plot", class(x$ind.plot), length(x$ind.plot), "name of the plot in which the individual is")
	sumry[3, ] <- c("$sp", class(x$sp), length(x$sp), "individuals' groups (e.g. species)")

	class(sumry) <- "table"
	cat("\n###############")
	cat("\ndata used\n")
	print(sumry)

	cat("\n###############\n")
	cat("others")
	cat("\n\t$namestraits:", length(x$namestraits),"traits\n")
	print(head(x$namestraits))
	if(length(x$namestraits)>5){print("...")}

	rich <- x$sites_richness
	class(rich) <- "table"
	cat("\n\t$sites_richness:\n\t")
	print(rich)
	cat("\n")
}

print.ComIndexMulti <- function(x, ...){

	if (!inherits(x, "ComIndexMulti"))
		{stop("x must be a list of objects of class ComIndexMulti")
	}

	cat("\t#################################\n")
	cat("\t# Community metrics calculation #\n")
	cat("\t#################################\n")
	cat("class: ")
	cat(class(x))

	cat("\n$call: ")
	print(x$call)

	cat("\n###############")
	cat("\n$obs: list of observed values\n")
	seq.interm <- seq(1,length(names(x$list.index)), by = 2)
	cat(paste("\t$", names(x$list.index)[seq.interm], sep = ""), sep = "\n")

	if(!is.null(x$null)){
		cat("\n###############")
		cat("\n$null: list of null values, number of permutations:", dim(x$null[[1]])[3], "\n")
		cat(paste(paste("\t$", names(x$list.index)[-seq.interm], sep = ""), paste("null model", "=", x$nullmodels, sep=" "), sep=" ... "), sep = "\n")
	}

	TABDIM <- 3

	sumry <- array("", c(TABDIM, 4), list(1:TABDIM, c("data", "class", "dim", "content")))
	sumry[1, ] <- c("$traits", class(x$traits), paste(dim(x$traits)[1], dim(x$traits)[2], sep = ",") , "traits data")
	sumry[2, ] <- c("$ind.plot", class(x$ind.plot), length(x$ind.plot), "name of the plot in which the individual is")
	sumry[3, ] <- c("$sp", class(x$sp), length(x$sp), "individuals' groups (e.g. species)")

	class(sumry) <- "table"
	cat("\n###############")
	cat("\ndata used\n")
	print(sumry)

	cat("\n###############\n")
	cat("others")
	cat("\n\t$namestraits:", length(x$namestraits),"traits\n")
	print(head(x$namestraits))
	if(length(x$namestraits)>5){print("...")}

	rich <- x$sites_richness
	class(rich) <- "table"
	cat("\n\t$sites_richness:\n\t")
	print(rich)
	cat("\n")
}

summary.Tstats <- function(object, ...){

	if (!inherits(object, "Tstats"))
		{stop("object must be a list of objects of class Tstats")
	}

	print("Observed values", ...)
	print(lapply(object$Tstats[c(1:3)], function(x) summary(x, ...)), ...)

	print("null values", ...)
	print(lapply(object$Tstats[c(4:6)], function(x) summary(x, ...)), ...)

}

summary.ComIndex <- function(object, ...){

	if (!inherits(object, "ComIndex"))
		{stop("object must be a list of objects of class ComIndex")
	}

	print("Observed values", ...)
	print(lapply(object$obs, function(x) summary(x, ...)), ...)

	print("null values", ...)
	print(lapply(object$null, function(x) summary(x, ...)), ...)

}

summary.ComIndexMulti <- function(object, ...){

	if (!inherits(object, "ComIndexMulti"))
		{stop("object must be a list of objects of class ComIndexMulti")
	}

	print("Observed values", ...)
	print(lapply(object$obs, function(x) summary(x, ...)), ...)

	print("null values", ...)
	print(lapply(object$null, function(x) summary(x, ...)), ...)

}

### Function to represent standardised effect size of Tstats using null models
plot.Tstats <- function(x, type = "normal", col.index = c("red","purple","olivedrab3"), add.conf = TRUE, color.cond = TRUE, val.quant = c(0.025,0.975), ...){

	type <- match.arg(type, c("simple", "simple_range", "normal", "barplot", "bysites", "bytraits"))

	if (!inherits(x, "Tstats")) {
		stop("x must be of class Tstats")
	}

	if(length(x$Tstats) != 6){
		stop("The object x of class Tstats does not contain null value from null models.")
	}

	plot(as.listofindex(x, namesindex = c("T_IP.IC", "T_IC.IR", "T_PC.PR")), type = type, col.index = col.index, add.conf = add.conf, color.cond = color.cond, val.quant = val.quant, ...)
}


### Function to summarize traits and community which show a significant difference between observed and simulated value
sum_Tstats <- function(x, val.quant = c(0.025,0.975), type = "all") {

	res <- x
	Tst <- res$Tstats

	#________________________________________
	ses.T_IP.IC <- (Tst$T_IP.IC-apply(Tst$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm = T))
	ses.T_IC.IR <- (Tst$T_IC.IR-apply(Tst$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm = T))
	ses.T_PC.PR <- (Tst$T_PC.PR-apply(Tst$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm = T))

	ses.T_IP.IC.inf <- (apply(Tst$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm = T, prob = val.quant[1]))-apply(Tst$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm = T))
	ses.T_IC.IR.inf <- (apply(Tst$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm = T, prob = val.quant[1]))-apply(Tst$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm = T))
	ses.T_PC.PR.inf <- (apply(Tst$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm = T, prob = val.quant[1]))-apply(Tst$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm = T))

	ses.T_IP.IC.sup <- (apply(Tst$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm = T, prob = val.quant[2]))-apply(Tst$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm = T))
	ses.T_IC.IR.sup <- (apply(Tst$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm = T, prob = val.quant[2]))-apply(Tst$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm = T))
	ses.T_PC.PR.sup <- (apply(Tst$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm = T, prob = val.quant[2]))-apply(Tst$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm = T))

	ses.T_IP.IC.mean <- t(colMeans((Tst$T_IP.IC-apply(Tst$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm = T)), na.rm = T))
	ses.T_IC.IR.mean <- t(colMeans((Tst$T_IC.IR-apply(Tst$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm = T)), na.rm = T))
	ses.T_PC.PR.mean <- t(colMeans((Tst$T_PC.PR-apply(Tst$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm = T)), na.rm = T))

	ses.T_IP.IC.inf.mean <- apply(ses.T_IP.IC.inf,2, mean)
	ses.T_IC.IR.inf.mean <- apply(ses.T_IC.IR.inf,2, mean)
	ses.T_PC.PR.inf.mean <- apply(ses.T_PC.PR.inf,2, mean)

	ses.T_IP.IC.sup.mean <- apply(ses.T_IP.IC.sup,2, mean)
	ses.T_IC.IR.sup.mean <- apply(ses.T_IC.IR.sup,2, mean)
	ses.T_PC.PR.sup.mean <- apply(ses.T_PC.PR.sup,2, mean)

	#________________________________________
	#Condition to be significantly different from null models with respect to values of quantile choosen
	cond.T_IP.IC.inf <- ses.T_IP.IC<ses.T_IP.IC.inf
	cond.T_IC.IR.inf <- ses.T_IC.IR<ses.T_IC.IR.inf
	cond.T_PC.PR.inf <- ses.T_PC.PR<ses.T_PC.PR.inf

	cond.T_IP.IC.sup <- ses.T_IP.IC>ses.T_IP.IC.sup
	cond.T_IC.IR.sup <- ses.T_IC.IR>ses.T_IC.IR.sup
	cond.T_PC.PR.sup <- ses.T_PC.PR>ses.T_PC.PR.sup

	cond.T_IP.IC.inf.mean <- ses.T_IP.IC.mean<ses.T_IP.IC.inf.mean
	cond.T_IC.IR.inf.mean <- ses.T_IC.IR.mean<ses.T_IC.IR.inf.mean
	cond.T_PC.PR.inf.mean <- ses.T_PC.PR.mean<ses.T_PC.PR.inf.mean

	cond.T_IP.IC.sup.mean <- ses.T_IP.IC.mean>ses.T_IP.IC.sup.mean
	cond.T_IC.IR.sup.mean <- ses.T_IC.IR.mean>ses.T_IC.IR.sup.mean
	cond.T_PC.PR.sup.mean <- ses.T_PC.PR.mean>ses.T_PC.PR.sup.mean


	#________________________________________

 	#########################################
	####		 calculation of p.value		 ####
	#########################################

	if (type == "all" | type == "p.value"){

 		p.valueT_IP.IC.sup <- res$Tstats$T_IP.IC
		p.valueT_IC.IR.sup <- res$Tstats$T_IP.IC
		p.valueT_PC.PR.sup <- res$Tstats$T_IP.IC

		p.valueT_IP.IC.inf <- res$Tstats$T_IP.IC
		p.valueT_IC.IR.inf <- res$Tstats$T_IP.IC
		p.valueT_PC.PR.inf <- res$Tstats$T_IP.IC

		for (t in 1: ncol(res$Tstats$T_IP.IC)){
			for(s in 1: nrow(res$Tstats$T_IP.IC)){
 				p.valueT_IP.IC.sup[s,t] <- (sum(res$Tstats$T_IP.IC[s,t]<res$Tstats$T_IP.IC_nm[,t,s], na.rm = T)+1)/(1+length(res$Tstats$T_IP.IC_nm[,t,s]))
 				p.valueT_IC.IR.sup[s,t] <- (sum(res$Tstats$T_IC.IR[s,t]<res$Tstats$T_IC.IR_nm[,t,s], na.rm = T)+1)/(1+length(res$Tstats$T_IC.IR_nm[,t,s]))
 				p.valueT_PC.PR.sup[s,t] <- (sum(res$Tstats$T_PC.PR[s,t]<res$Tstats$T_PC.PR_nm[,t,s], na.rm = T)+1)/(1+length(res$Tstats$T_PC.PR_nm[,t,s]))

				p.valueT_IP.IC.inf[s,t] <- (sum(res$Tstats$T_IP.IC[s,t]>res$Tstats$T_IP.IC_nm[,t,s], na.rm = T)+1)/(1+length(res$Tstats$T_IP.IC_nm[,t,s]))
				p.valueT_IC.IR.inf[s,t] <- (sum(res$Tstats$T_IC.IR[s,t]>res$Tstats$T_IC.IR_nm[,t,s], na.rm = T)+1)/(1+length(res$Tstats$T_IC.IR_nm[,t,s]))
				p.valueT_PC.PR.inf[s,t] <- (sum(res$Tstats$T_PC.PR[s,t]>res$Tstats$T_PC.PR_nm[,t,s], na.rm = T)+1)/(1+length(res$Tstats$T_PC.PR_nm[,t,s]))
			}
		}

		colnames(p.valueT_IP.IC.sup) <- colnames(res$Tstats$T_IP.IC)
		colnames(p.valueT_IC.IR.sup) <- colnames(res$Tstats$T_IP.IC)
		colnames(p.valueT_PC.PR.sup) <- colnames(res$Tstats$T_IP.IC)

		rownames(p.valueT_IP.IC.sup) <- rownames(res$Tstats$T_IP.IC)
		rownames(p.valueT_IC.IR.sup) <- rownames(res$Tstats$T_IP.IC)
		rownames(p.valueT_PC.PR.sup) <- rownames(res$Tstats$T_IP.IC)

		colnames(p.valueT_IP.IC.inf) <- colnames(res$Tstats$T_IP.IC)
		colnames(p.valueT_IC.IR.inf) <- colnames(res$Tstats$T_IP.IC)
		colnames(p.valueT_PC.PR.inf) <- colnames(res$Tstats$T_IP.IC)

		rownames(p.valueT_IP.IC.inf) <- rownames(res$Tstats$T_IP.IC)
		rownames(p.valueT_IC.IR.inf) <- rownames(res$Tstats$T_IP.IC)
		rownames(p.valueT_PC.PR.inf) <- rownames(res$Tstats$T_IP.IC)

		pval <- list()

		pval$T_IP.IC.inf <- p.valueT_IP.IC.inf
		pval$T_IC.IR.inf <- p.valueT_IC.IR.inf
		pval$T_PC.PR.inf <- p.valueT_PC.PR.inf

		pval$T_IP.IC.sup <- p.valueT_IP.IC.sup
		pval$T_IC.IR.sup <- p.valueT_IC.IR.sup
		pval$T_PC.PR.sup <- p.valueT_PC.PR.sup
	}
	else{}

	#________________________________________
	if (type == "binary"){
		summ.Tstats <- matrix("H0 not rejected",nrow = 6, ncol = dim(cond.T_IP.IC.inf)[2])
		summ.Tstats <- rbind(cond.T_IP.IC.inf.mean, cond.T_IP.IC.sup.mean ,cond.T_IC.IR.inf.mean, cond.T_IC.IR.sup.mean ,cond.T_PC.PR.inf.mean, cond.T_IC.IR.sup.mean)
		rownames(summ.Tstats) <- c("T_IP.IC.inf", "T_IP.IC.sup", "T_IC.IR.inf", "T_IC.IR.sup", "T_PC.PR.inf", "T_PC.PR.sup")
		colnames(summ.Tstats) <- colnames(Tst$T_IP.IC)
	}

	#________________________________________
	else if (type == "percent"){

		summ.Tstats <- matrix("H0 not rejected",nrow = 6, ncol = dim(cond.T_IP.IC.inf)[2])

		summ.Tstats[1, ] <- paste(round(colSums(cond.T_IP.IC.inf, na.rm = T)/colSums(!is.na(cond.T_IP.IC.inf)),2)*100, "%", sep = "")
		summ.Tstats[2, ] <- paste(round(colSums(cond.T_IP.IC.sup, na.rm = T)/colSums(!is.na(cond.T_IP.IC.sup)),2)*100, "%", sep = "")
		summ.Tstats[3, ] <- paste(round(colSums(cond.T_IC.IR.inf, na.rm = T)/colSums(!is.na(cond.T_IC.IR.inf)),2)*100, "%", sep = "")
		summ.Tstats[4, ] <- paste(round(colSums(cond.T_IC.IR.sup, na.rm = T)/colSums(!is.na(cond.T_IC.IR.sup)),2)*100, "%", sep = "")
		summ.Tstats[5, ] <- paste(round(colSums(cond.T_PC.PR.inf, na.rm = T)/colSums(!is.na(cond.T_PC.PR.inf)),2)*100, "%", sep = "")
		summ.Tstats[6, ] <- paste(round(colSums(cond.T_PC.PR.sup, na.rm = T)/colSums(!is.na(cond.T_PC.PR.sup)),2)*100, "%", sep = "")

		summ.Tstats[1, ][cond.T_IP.IC.inf.mean] <- paste(round(colSums(cond.T_IP.IC.inf, na.rm = T)/colSums(!is.na(cond.T_IP.IC.inf)),2)[cond.T_IP.IC.inf.mean]*100, "%" ,"*", sep = "")
		summ.Tstats[2, ][cond.T_IP.IC.sup.mean] <- paste(round(colSums(cond.T_IP.IC.sup, na.rm = T)/colSums(!is.na(cond.T_IP.IC.sup)),2)[cond.T_IP.IC.sup.mean]*100, "%" ,"*", sep = "")
		summ.Tstats[3, ][cond.T_IC.IR.inf.mean] <- paste(round(colSums(cond.T_IC.IR.inf, na.rm = T)/colSums(!is.na(cond.T_IC.IR.inf)),2)[cond.T_IC.IR.inf.mean]*100, "%","*", sep = "")
		summ.Tstats[4, ][cond.T_IC.IR.sup.mean] <- paste(round(colSums(cond.T_IC.IR.sup, na.rm = T)/colSums(!is.na(cond.T_IC.IR.sup)),2)[cond.T_IC.IR.sup.mean]*100, "%","*", sep = "")
		summ.Tstats[5, ][cond.T_PC.PR.inf.mean] <- paste(round(colSums(cond.T_PC.PR.inf, na.rm = T)/colSums(!is.na(cond.T_PC.PR.inf)),2)[cond.T_PC.PR.inf.mean]*100, "%","*", sep = "")
		summ.Tstats[6, ][cond.T_PC.PR.sup.mean] <- paste(round(colSums(cond.T_PC.PR.sup, na.rm = T)/colSums(!is.na(cond.T_PC.PR.sup)),2)[cond.T_PC.PR.sup.mean]*100, "%","*", sep = "")

		summ.Tstats[1,] <- paste(round(colSums(cond.T_IP.IC.inf, na.rm = T)/colSums(!is.na(cond.T_IP.IC.inf)),2)*100, "%", sep = "")
		summ.Tstats[2,] <- paste(round(colSums(cond.T_IP.IC.sup, na.rm = T)/colSums(!is.na(cond.T_IP.IC.sup)),2)*100, "%", sep = "")
		summ.Tstats[3,] <- paste(round(colSums(cond.T_IC.IR.inf, na.rm = T)/colSums(!is.na(cond.T_IC.IR.inf)),2)*100, "%", sep = "")
		summ.Tstats[4,] <- paste(round(colSums(cond.T_IC.IR.sup, na.rm = T)/colSums(!is.na(cond.T_IC.IR.sup)),2)*100, "%", sep = "")
		summ.Tstats[5,] <- paste(round(colSums(cond.T_PC.PR.inf, na.rm = T)/colSums(!is.na(cond.T_PC.PR.inf)),2)*100, "%", sep = "")
		summ.Tstats[6,] <- paste(round(colSums(cond.T_PC.PR.sup, na.rm = T)/colSums(!is.na(cond.T_PC.PR.sup)),2)*100, "%", sep = "")

		summ.Tstats[1,][cond.T_IP.IC.inf.mean] <- paste(round(colSums(cond.T_IP.IC.inf, na.rm = T)/colSums(!is.na(cond.T_IP.IC.inf)),2)[cond.T_IP.IC.inf.mean]*100, "%" ,"*", sep = "")
		summ.Tstats[2,][cond.T_IP.IC.sup.mean] <- paste(round(colSums(cond.T_IP.IC.sup, na.rm = T)/colSums(!is.na(cond.T_IP.IC.sup)),2)[cond.T_IP.IC.sup.mean]*100, "%" ,"*", sep = "")
		summ.Tstats[3,][cond.T_IC.IR.inf.mean] <- paste(round(colSums(cond.T_IC.IR.inf, na.rm = T)/colSums(!is.na(cond.T_IC.IR.inf)),2)[cond.T_IC.IR.inf.mean]*100, "%","*", sep = "")
		summ.Tstats[4,][cond.T_IC.IR.sup.mean] <- paste(round(colSums(cond.T_IC.IR.sup, na.rm = T)/colSums(!is.na(cond.T_IC.IR.sup)),2)[cond.T_IC.IR.sup.mean]*100, "%","*", sep = "")
		summ.Tstats[5,][cond.T_PC.PR.inf.mean] <- paste(round(colSums(cond.T_PC.PR.inf, na.rm = T)/colSums(!is.na(cond.T_PC.PR.inf)),2)[cond.T_PC.PR.inf.mean]*100, "%","*", sep = "")
		summ.Tstats[6,][cond.T_PC.PR.sup.mean] <- paste(round(colSums(cond.T_PC.PR.sup, na.rm = T)/colSums(!is.na(cond.T_PC.PR.sup)),2)[cond.T_PC.PR.sup.mean]*100, "%","*", sep = "")


		rownames(summ.Tstats) <- c("T_IP.IC.inf", "T_IP.IC.sup", "T_IC.IR.inf", "T_IC.IR.sup", "T_PC.PR.inf", "T_PC.PR.sup")
		colnames(summ.Tstats) <- colnames(Tst$T_IP.IC)

	}

	#________________________________________
	else if (type == "site"){

		summ.Tstats <- matrix("H0 not rejected",nrow = 6, ncol = dim(cond.T_IP.IC.inf)[2])
		for(t in 1: dim(cond.T_IP.IC.inf)[2]){

			if (sum(cond.T_IP.IC.inf[,t], na.rm = T)>0)
				{summ.Tstats[1,t] <- paste( na.exclude(rownames(cond.T_IP.IC.inf)[cond.T_IP.IC.inf[,t]]), collapse = " ") }
			else{summ.Tstats[1,t] <- "H0 not rejected"}

			if (sum(cond.T_IP.IC.sup[,t], na.rm = T)>0)
				{summ.Tstats[2,t] <- paste( na.exclude(rownames(cond.T_IP.IC.sup)[cond.T_IP.IC.sup[,t]]), collapse = " ") }
			else{summ.Tstats[2,t] <- "H0 not rejected"}

			if (sum(cond.T_IC.IR.inf[,t], na.rm = T)>0)
				{summ.Tstats[3,t] <- paste( na.exclude(rownames(cond.T_IC.IR.inf)[cond.T_IP.IC.inf[,t]]), collapse = " ") }
			else{summ.Tstats[3,t] <- "H0 not rejected"}

			if (sum(cond.T_IC.IR.sup[,t], na.rm = T)>0)
				{summ.Tstats[4,t] <- paste( na.exclude(rownames(cond.T_IC.IR.sup)[cond.T_IC.IR.sup[,t]]), collapse = " ")	}
			else{summ.Tstats[4,t] <- "H0 not rejected"}

			if (sum(cond.T_PC.PR.inf[,t], na.rm = T)>0)
				{summ.Tstats[5,t] <- paste( na.exclude(rownames(cond.T_PC.PR.inf)[cond.T_PC.PR.inf[,t]]), collapse = " ") 	}
			else{summ.Tstats[5,t] <- "H0 not rejected"}

			if (sum(cond.T_PC.PR.sup[,t], na.rm = T)>0)
				{summ.Tstats[6,t] <- paste( na.exclude(rownames(cond.T_PC.PR.sup)[cond.T_PC.PR.sup[,t]]), collapse = " ") 	}
			else{summ.Tstats[6,t] <- "H0 not rejected"}
		}
		rownames(summ.Tstats) <- c("T_IP.IC.inf", "T_IP.IC.sup", "T_IC.IR.inf", "T_IC.IR.sup", "T_PC.PR.inf", "T_PC.PR.sup")
		colnames(summ.Tstats) <- colnames(Tst$T_IP.IC)

	}


	#________________________________________
	else if (type == "p.value"){
		summ.Tstats <- matrix("H0 not rejected",nrow = 6, ncol = dim(cond.T_IP.IC.inf)[2])
		summ.Tstats <- rbind(pval$T_IP.IC.inf, pval$T_IP.IC.sup , pval$T_IC.IR.inf, pval$T_IC.IR.sup , pval$T_PC.PR.inf, pval$T_PC.PR.sup)
		rownames(summ.Tstats) <- c(paste(rep("T_IP.IC.inf",dim(Tst$T_IP.IC)[1]), rownames(Tst$T_IP.IC)), paste(rep("T_IP.IC.sup",dim(Tst$T_IP.IC)[1]), rownames(Tst$T_IP.IC)), paste(rep("T_IC.IR.inf",dim(Tst$T_IP.IC)[1]), rownames(Tst$T_IP.IC)), paste(rep("T_IC.IR.sup",dim(Tst$T_IP.IC)[1]), rownames(Tst$T_IP.IC)), paste(rep("T_PC.PR.inf",dim(Tst$T_IP.IC)[1]), rownames(Tst$T_IP.IC)), paste(rep("T_PC.PR.sup",dim(Tst$T_IP.IC)[1]), rownames(Tst$T_IP.IC)))
		colnames(summ.Tstats) <- colnames(Tst$T_IP.IC)
	}


	#________________________________________
	else if (type == "all"){
		summ.Tstats <- list()

		#__________
		##p.value
		summ.Tstats$p.value <- matrix("H0 not rejected", nrow = 6, ncol = dim(cond.T_IP.IC.inf)[2])
		summ.Tstats$p.value <- rbind(pval$T_IP.IC.inf, pval$T_IP.IC.sup , pval$T_IC.IR.inf, pval$T_IC.IR.sup , pval$T_PC.PR.inf, pval$T_PC.PR.sup)

		#__________
		##percent
		summ.Tstats$percent <- matrix("H0 not rejected",nrow = 6, ncol = dim(cond.T_IP.IC.inf)[2])

		summ.Tstats$percent[1, ] <- paste(round(colSums(cond.T_IP.IC.inf, na.rm = T)/colSums(!is.na(cond.T_IP.IC.inf)),2)*100, "%", sep = "")
		summ.Tstats$percent[2, ] <- paste(round(colSums(cond.T_IP.IC.sup, na.rm = T)/colSums(!is.na(cond.T_IP.IC.sup)),2)*100, "%", sep = "")
		summ.Tstats$percent[3, ] <- paste(round(colSums(cond.T_IC.IR.inf, na.rm = T)/colSums(!is.na(cond.T_IC.IR.inf)),2)*100, "%", sep = "")
		summ.Tstats$percent[4, ] <- paste(round(colSums(cond.T_IC.IR.sup, na.rm = T)/colSums(!is.na(cond.T_IC.IR.sup)),2)*100, "%", sep = "")
		summ.Tstats$percent[5, ] <- paste(round(colSums(cond.T_PC.PR.inf, na.rm = T)/colSums(!is.na(cond.T_PC.PR.inf)),2)*100, "%", sep = "")
		summ.Tstats$percent[6, ] <- paste(round(colSums(cond.T_PC.PR.sup, na.rm = T)/colSums(!is.na(cond.T_PC.PR.sup)),2)*100, "%", sep = "")

		summ.Tstats$percent[1, ][cond.T_IP.IC.inf.mean] <- paste(round(colSums(cond.T_IP.IC.inf, na.rm = T)/colSums(!is.na(cond.T_IP.IC.inf)),2)[cond.T_IP.IC.inf.mean]*100, "%" ,"*", sep = "")
		summ.Tstats$percent[2, ][cond.T_IP.IC.sup.mean] <- paste(round(colSums(cond.T_IP.IC.sup, na.rm = T)/colSums(!is.na(cond.T_IP.IC.sup)),2)[cond.T_IP.IC.sup.mean]*100, "%" ,"*", sep = "")
		summ.Tstats$percent[3, ][cond.T_IC.IR.inf.mean] <- paste(round(colSums(cond.T_IC.IR.inf, na.rm = T)/colSums(!is.na(cond.T_IC.IR.inf)),2)[cond.T_IC.IR.inf.mean]*100, "%","*", sep = "")
		summ.Tstats$percent[4, ][cond.T_IC.IR.sup.mean] <- paste(round(colSums(cond.T_IC.IR.sup, na.rm = T)/colSums(!is.na(cond.T_IC.IR.sup)),2)[cond.T_IC.IR.sup.mean]*100, "%","*", sep = "")
		summ.Tstats$percent[5, ][cond.T_PC.PR.inf.mean] <- paste(round(colSums(cond.T_PC.PR.inf, na.rm = T)/colSums(!is.na(cond.T_PC.PR.inf)),2)[cond.T_PC.PR.inf.mean]*100, "%","*", sep = "")
		summ.Tstats$percent[6, ][cond.T_PC.PR.sup.mean] <- paste(round(colSums(cond.T_PC.PR.sup, na.rm = T)/colSums(!is.na(cond.T_PC.PR.sup)),2)[cond.T_PC.PR.sup.mean]*100, "%","*", sep = "")

		summ.Tstats$percent[1,] <- paste(round(colSums(cond.T_IP.IC.inf, na.rm = T)/colSums(!is.na(cond.T_IP.IC.inf)),2)*100, "%", sep = "")
		summ.Tstats$percent[2,] <- paste(round(colSums(cond.T_IP.IC.sup, na.rm = T)/colSums(!is.na(cond.T_IP.IC.sup)),2)*100, "%", sep = "")
		summ.Tstats$percent[3,] <- paste(round(colSums(cond.T_IC.IR.inf, na.rm = T)/colSums(!is.na(cond.T_IC.IR.inf)),2)*100, "%", sep = "")
		summ.Tstats$percent[4,] <- paste(round(colSums(cond.T_IC.IR.sup, na.rm = T)/colSums(!is.na(cond.T_IC.IR.sup)),2)*100, "%", sep = "")
		summ.Tstats$percent[5,] <- paste(round(colSums(cond.T_PC.PR.inf, na.rm = T)/colSums(!is.na(cond.T_PC.PR.inf)),2)*100, "%", sep = "")
		summ.Tstats$percent[6,] <- paste(round(colSums(cond.T_PC.PR.sup, na.rm = T)/colSums(!is.na(cond.T_PC.PR.sup)),2)*100, "%", sep = "")

		summ.Tstats$percent[1,][cond.T_IP.IC.inf.mean] <- paste(round(colSums(cond.T_IP.IC.inf, na.rm = T)/colSums(!is.na(cond.T_IP.IC.inf)),2)[cond.T_IP.IC.inf.mean]*100, "%" ,"*", sep = "")
		summ.Tstats$percent[2,][cond.T_IP.IC.sup.mean] <- paste(round(colSums(cond.T_IP.IC.sup, na.rm = T)/colSums(!is.na(cond.T_IP.IC.sup)),2)[cond.T_IP.IC.sup.mean]*100, "%" ,"*", sep = "")
		summ.Tstats$percent[3,][cond.T_IC.IR.inf.mean] <- paste(round(colSums(cond.T_IC.IR.inf, na.rm = T)/colSums(!is.na(cond.T_IC.IR.inf)),2)[cond.T_IC.IR.inf.mean]*100, "%","*", sep = "")
		summ.Tstats$percent[4,][cond.T_IC.IR.sup.mean] <- paste(round(colSums(cond.T_IC.IR.sup, na.rm = T)/colSums(!is.na(cond.T_IC.IR.sup)),2)[cond.T_IC.IR.sup.mean]*100, "%","*", sep = "")
		summ.Tstats$percent[5,][cond.T_PC.PR.inf.mean] <- paste(round(colSums(cond.T_PC.PR.inf, na.rm = T)/colSums(!is.na(cond.T_PC.PR.inf)),2)[cond.T_PC.PR.inf.mean]*100, "%","*", sep = "")
		summ.Tstats$percent[6,][cond.T_PC.PR.sup.mean] <- paste(round(colSums(cond.T_PC.PR.sup, na.rm = T)/colSums(!is.na(cond.T_PC.PR.sup)),2)[cond.T_PC.PR.sup.mean]*100, "%","*", sep = "")


		#__________
		##sites
		summ.Tstats$sites <- matrix("H0 not rejected",nrow = 6, ncol = dim(cond.T_IP.IC.inf)[2])
		for(t in 1: dim(cond.T_IP.IC.inf)[2]){

			if (sum(cond.T_IP.IC.inf[,t], na.rm = T)>0)
				{summ.Tstats$sites[1,t] <- paste( na.exclude(rownames(cond.T_IP.IC.inf)[cond.T_IP.IC.inf[,t]]), collapse = " ") }
			else{summ.Tstats$sites[1,t] <- "H0 not rejected"}

			if (sum(cond.T_IP.IC.sup[,t], na.rm = T)>0)
				{summ.Tstats$sites[2,t] <- paste( na.exclude(rownames(cond.T_IP.IC.sup)[cond.T_IP.IC.sup[,t]]), collapse = " ") }
			else{summ.Tstats$sites[2,t] <- "H0 not rejected"}

			if (sum(cond.T_IC.IR.inf[,t], na.rm = T)>0)
				{summ.Tstats$sites[3,t] <- paste( na.exclude(rownames(cond.T_IC.IR.inf)[cond.T_IP.IC.inf[,t]]), collapse = " ") }
			else{summ.Tstats$sites[3,t] <- "H0 not rejected"}

			if (sum(cond.T_IC.IR.sup[,t], na.rm = T)>0)
				{summ.Tstats$sites[4,t] <- paste( na.exclude(rownames(cond.T_IC.IR.sup)[cond.T_IC.IR.sup[,t]]), collapse = " ")	}
			else{summ.Tstats$sites[4,t] <- "H0 not rejected"}

			if (sum(cond.T_PC.PR.inf[,t], na.rm = T)>0)
				{summ.Tstats$sites[5,t] <- paste( na.exclude(rownames(cond.T_PC.PR.inf)[cond.T_PC.PR.inf[,t]]), collapse = " ") 	}
			else{summ.Tstats$sites[5,t] <- "H0 not rejected"}

			if (sum(cond.T_PC.PR.sup[,t], na.rm = T)>0)
				{summ.Tstats$sites[6,t] <- paste( na.exclude(rownames(cond.T_PC.PR.sup)[cond.T_PC.PR.sup[,t]]), collapse = " ") 	}
			else{summ.Tstats$sites[6,t] <- "H0 not rejected"}
		}

		#__________
		##binary
		summ.Tstats$binary <- matrix("H0 not rejected",nrow = 6, ncol = dim(cond.T_IP.IC.inf)[2])
		summ.Tstats$binary <- rbind(cond.T_IP.IC.inf.mean, cond.T_IP.IC.sup.mean ,cond.T_IC.IR.inf.mean, cond.T_IC.IR.sup.mean ,cond.T_PC.PR.inf.mean, cond.T_IC.IR.sup.mean)

		#__________
		rownames(summ.Tstats$p.value) <- c(rep("T_IP.IC.inf",dim(Tst$T_IP.IC)[1]), rep("T_IP.IC.sup",dim(Tst$T_IP.IC)[1]), rep("T_IC.IR.inf",dim(Tst$T_IP.IC)[1]), rep("T_IC.IR.sup",dim(Tst$T_IP.IC)[1]), rep("T_PC.PR.inf",dim(Tst$T_IP.IC)[1]), rep("T_PC.PR.sup",dim(Tst$T_IP.IC)[1]))
		rownames(summ.Tstats$binary) <- c("T_IP.IC.inf", "T_IP.IC.sup", "T_IC.IR.inf", "T_IC.IR.sup", "T_PC.PR.inf", "T_PC.PR.sup")
		rownames(summ.Tstats$percent) <- c("T_IP.IC.inf", "T_IP.IC.sup", "T_IC.IR.inf", "T_IC.IR.sup", "T_PC.PR.inf", "T_PC.PR.sup")
		rownames(summ.Tstats$sites) <- c("T_IP.IC.inf", "T_IP.IC.sup", "T_IC.IR.inf", "T_IC.IR.sup", "T_PC.PR.inf", "T_PC.PR.sup")
		colnames(summ.Tstats$p.value) <- colnames(Tst$T_IP.IC)
		colnames(summ.Tstats$binary) <- colnames(Tst$T_IP.IC)
		colnames(summ.Tstats$sites) <- colnames(Tst$T_IP.IC)
		colnames(summ.Tstats$percent) <- colnames(Tst$T_IP.IC)
	}

	else{stop("Error: type must be 'binary', 'percent', 'p.value', 'site' or 'all'.")}

	return(summ.Tstats)
}

### Function to represent summarize Tstats
barplot.Tstats <- function(height, val.quant = c(0.025,0.975), col.index = c("red","purple","olivedrab3","white"), ylim = NULL, ...){

	Tst <- height$Tstats

	T_IP.IC.inf <- apply(Tst$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm = T, prob = val.quant[1]))
	T_IC.IR.inf <- apply(Tst$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm = T, prob = val.quant[1]))
	T_PC.PR.inf <- apply(Tst$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm = T, prob = val.quant[1]))

	T_IP.IC.sup <- apply(Tst$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm = T, prob = val.quant[2]))
	T_IC.IR.sup <- apply(Tst$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm = T, prob = val.quant[2]))
	T_PC.PR.sup <- apply(Tst$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm = T, prob = val.quant[2]))

	if (is.null(ylim)){
		ylim = c(min(c(T_IP.IC.inf,T_IC.IR.inf,T_PC.PR.inf,colMeans(na.omit(Tst$T_IP.IC))-apply(na.omit(Tst$T_IP.IC), 2,sd),colMeans(na.omit(Tst$T_IC.IR))-apply(na.omit(Tst$T_IC.IR), 2,sd),colMeans(na.omit(Tst$T_PC.PR))-apply(na.omit(Tst$T_PC.PR), 2,sd)), na.rm = T) , max(c(colMeans(na.omit(Tst$T_IP.IC))+apply(na.omit(Tst$T_IP.IC), 2,sd),colMeans(na.omit(Tst$T_IC.IR))+apply(na.omit(Tst$T_IC.IR), 2,sd),colMeans(na.omit(Tst$T_PC.PR))+apply(na.omit(Tst$T_PC.PR), 2,sd)), na.rm = T))
	}

	df.bar <- barplot(rbind(colMeans(na.omit(Tst$T_IP.IC)), colMeans(na.omit(Tst$T_IC.IR)),colMeans(na.omit(Tst$T_PC.PR)),0), beside = T, plot = F)
	barplot(rbind(colMeans(na.omit(Tst$T_IP.IC)), colMeans(na.omit(Tst$T_IC.IR)),colMeans(na.omit(Tst$T_PC.PR)),0), col = col.index, beside = T, ylim = ylim, ...)
	segments( df.bar[1, ], colMeans(na.omit(Tst$T_IP.IC))+apply(na.omit(Tst$T_IP.IC), 2,sd),df.bar[1, ],colMeans(na.omit(Tst$T_IP.IC))-apply(na.omit(Tst$T_IP.IC), 2,sd))
	segments( df.bar[2, ], colMeans(na.omit(Tst$T_IC.IR))+apply(na.omit(Tst$T_IC.IR), 2,sd),df.bar[2, ],colMeans(na.omit(Tst$T_IC.IR))-apply(na.omit(Tst$T_IC.IR), 2,sd))
	segments( df.bar[3, ], colMeans(na.omit(Tst$T_PC.PR))+apply(na.omit(Tst$T_PC.PR), 2,sd),df.bar[3, ],colMeans(na.omit(Tst$T_PC.PR))-apply(na.omit(Tst$T_PC.PR), 2,sd))

	points(type = "l", df.bar[1, ], colMeans(T_IP.IC.sup, na.rm = T), col = col.index[1])
	points(type = "l", df.bar[2, ], colMeans(T_IC.IR.sup, na.rm = T), col = col.index[2])
	points(type = "l", df.bar[3, ], colMeans(T_PC.PR.sup, na.rm = T), col = col.index[3])

	points(type = "l", df.bar[1, ], colMeans(T_IP.IC.inf, na.rm = T), col = col.index[1])
	points(type = "l", df.bar[2, ], colMeans(T_IC.IR.inf, na.rm = T), col = col.index[2])
	points(type = "l", df.bar[3, ], colMeans(T_PC.PR.inf, na.rm = T), col = col.index[3])

}




#______________#______________#______________#______________#______________#______________#______________#______________
#______________#______________#______________#______________#______________#______________#______________#______________
#__ ComIndex

#calculation of statistics (e.g. mean, range, CVNND and kurtosis) to test community assembly using null models
#For each statistic this function return observed value and correspondant null distribution
#This function implement three null models which keep unchanged the number of individual per community
#Models local (1) correspond to randomization of individual values within community
#Models regional.ind (2) correspond to randomization of individual values within region
#Models regional.pop (2sp) correspond to randomization of population values within region
#Models regional.pop.prab (2sp.prab) correspond to randomization of population values within region whithout taking abundance into account (prab = presence/absence)

#In most case, model local and regional.ind correspond to index at the individual level and models regional.pop and regional.pop.prab to index at the species (or any other aggregate variable like genus or family) level

ComIndex <- function(traits = NULL, index = NULL, nullmodels = NULL, ind.plot = NULL, sp = NULL, com = NULL, SE = 0, namesindex = NULL, reg.pool = NULL, SE.reg.pool = NULL, nperm = 99, printprogress = TRUE, type.sp.val = "count"){

	type.sp.val <- match.arg(type.sp.val , c("count", "abundance"))

	#Create/verify/sort variables to put into the function
	nindex <- length(index)

	if (is.null(namesindex)) { namesindex <- index }
	if (!is.matrix(traits)) { traits <- as.matrix(traits)}
	if (is.null(index)) { stop("There is no default index to compute. Please add metrics using the 'index' argument.") }


	ntr <- dim(traits)[2]
	namestraits <- colnames(traits)

	#If data are from species or population traits, this function (AbToInd) transform this data in a suitable format for cati (individual like data)

	if (is.null(ind.plot)){
		if (is.null(com)) {stop("Need 'ind.plot' or 'com' argument")}

		rownames(traits) <- sp
		res.interm <- AbToInd(traits, com, type.sp.val = type.sp.val)

		traits <- res.interm$traits
		sp <- res.interm$sp
		ind.plot <- res.interm$ind.plot
	}

	if (!is.null(ind.plot) & !is.null(com)){
		warning("If ind.plot and com are provide, the function use only the argument ind.plot")
	}

	traits <- traits[order(ind.plot),]

	ind.plot <- ind.plot[order(ind.plot)]
	sp <- sp[order(ind.plot)]

	name_sp_sites <- paste(sp, ind.plot, sep = "_")
	comm <- NULL
	comm <- t(table(ind.plot,1:length(ind.plot)))

	S <- colSums(comm>0)
	ncom <- length(S)

	if (sum(is.na(traits))>0) {
		warning(paste("This function exclude", sum(is.na(traits)),"Na value", sep=" "))
	}

	if(length(SE) == 1){
			SE <- rep(SE, times = ncol(traits))
	}

	#Null models argument preparation
	if (is.null(nullmodels)) {
		nperm <- NULL
		warning("nullmodels is NULL so no null model will be computed")
	}
	
	if (is.null(nperm)) {
		nullmodels <- NULL
		warning("nperm is NULL so no null model will be computed")
	}

	if(!is.null(nperm)){
		if (nperm == 0) {
			nperm <- NULL
			nullmodels <- NULL
			message("nperm = 0, so no null model will be computed")
		}
	}

	if (!is.null(nullmodels) & !is.null(nperm)){

		if (length(nullmodels) == 1){
			nullmodels <- rep(nullmodels,times = nindex)
		}

		nullmodels[nullmodels == "1"] <- "local"
		nullmodels[nullmodels == "2"] <- "regional.ind"
		nullmodels[nullmodels == "2sp"] <- "regional.pop"
		nullmodels[nullmodels == "2sp.prab"] <- "regional.pop.prab"

		nullmodels <- match.arg(nullmodels, c("local", "regional.ind", "regional.pop", "regional.pop.prab"), several.ok = TRUE)
		nullmodels_names <- nullmodels
	}


	###############################################################
	if (is.numeric(nperm)){

		#regional pool
		if (!is.null(reg.pool) & sum(nullmodels == "regional.ind") == 0) {
			warning("Custom regional pool 'reg.pool' is only used in the case of null model regional.ind")
		}

		#########################################
		#### 	 Creating regional pools 	 ####
		#########################################

		if(!is.null(reg.pool)){isnullregpool <- FALSE} else{isnullregpool <- TRUE}

		# If the regional pool is the same for all communities:
		# creating a list of regional pool (one regional pool by community)
		if (is.null(reg.pool)) {
			reg.pool <- rep(list(traits), nlevels(ind.plot))
		}

		if (is.data.frame(reg.pool) | is.matrix(reg.pool) ) {
			reg.pool <- rep(list(reg.pool), nlevels(ind.plot))
		}

		# Standard error for regional pool
		if(is.null(SE.reg.pool)){
			SE.reg.pool <- SE
		}

		if(length(SE.reg.pool) == 1){
			SE.reg.pool <- rep(SE.reg.pool, times = ncol(traits))
		}

		if (is.vector(SE.reg.pool)) {
			SE.reg.pool <- rep(list(SE.reg.pool), nlevels(ind.plot))
		}

		if(length(reg.pool) != nlevels(ind.plot)){
			stop("reg.pool need to be either a matrix or a list of length equal to the number of communities")
		}

		#________________________________________
		# Warnings and stop about standard errors parameter
		if (length(SE) != ncol(traits) & length(SE) != 1) {
			stop("The vector SE need to have a length of one or equal to the number of traits")
		}

		if(!isnullregpool){
			if (length(SE.reg.pool) != length(reg.pool)) {
				stop("The vector SE.reg.pool need to have the same dimension as reg.pool")
			}
		}

		nullmodels[nullmodels == "local"] <- "1"
		nullmodels[nullmodels == "regional.ind"] <- "2"
		nullmodels[nullmodels == "regional.pop"] <- "2sp"
		nullmodels[nullmodels == "regional.pop.prab"] <- "2sp.prab"


		#########################################
		#### 	 calculation of null models 	 ####
		#########################################
		#Creation of three null models
		if (printprogress == T){ print("creating null models")}

		if (sum(nullmodels == "1")>0){
			#________________________________________
			#null model local: sample individuals traits values within community

			traits.nm1 <- list()

			for (t in 1: ntr){
				traits.nm1[[eval(namestraits[t])]] <- matrix(NA, nrow = dim(traits)[1], ncol = nperm)
				perm_ind.plot <- list()

				for(n in 1:nperm){
					for(s in 1: ncom) {
						perm_ind.plot[[s]] <- sample(traits[ind.plot == levels(ind.plot)[s], t], table(ind.plot)[s])
					}

					traits.nm1[[eval(namestraits[t])]][,n] <- unlist(perm_ind.plot)
				}
				if (printprogress == T){
					print(paste("local",round(t/ntr*100,2),"%"))
				}
			}
		}


		if (sum(nullmodels == "2")>0){
			#________________________________________
			#null model regional.ind: sample individuals traits values within the region (among community)
			traits.nm2 <- list()

			for (t in 1: ntr){
				traits.nm2[[eval(namestraits[t])]] <- matrix(NA, nrow = dim(traits)[1], ncol = nperm)
				perm_ind.plot <- list()

				# Take measurement standard error into account
				for(s in 1: ncom) {
					trait.intern <- reg.pool[[s]][, t]
					SE.reg.pool.intern <- SE.reg.pool[[s]][t]
					if (SE.reg.pool.intern != 0) {
						trait.intern <- rnorm(length(trait.intern), mean = trait.intern, sd = SE.reg.pool.intern)
					}
				}

				for(n in 1:nperm){
					for(s in 1: ncom) {
						perm_ind.plot[[s]] <- sample(trait.intern, table(ind.plot)[s])
					}

					traits.nm2[[eval(namestraits[t])]][,n] <- unlist(perm_ind.plot)
				}
				if (printprogress == T){
					print(paste("regional.ind",round(t/ntr*100,2),"%"))
				}
			}
		}


		if (sum(nullmodels == "2sp")>0){
			#________________________________________
			#null model regional.pop: sample populational traits values within the region (among community)

			traits.nm2sp <- list()
			traits_by_sp <- apply(traits,2,function(x) tapply(x,name_sp_sites,mean, na.rm = T))

			traits_by_pop <- traits_by_sp[match(name_sp_sites, rownames(traits_by_sp)), ]

			for (t in 1: ntr){
				traits.nm2sp[[eval(namestraits[t])]] <- matrix(NA, nrow = dim(traits)[1], ncol = nperm)
				perm_ind.plot <- list()

				for(n in 1:nperm){
					for(s in 1: ncom) {
						perm_ind.plot[[s]] <- sample(traits_by_pop, table(ind.plot)[s])
					}

					traits.nm2sp[[eval(namestraits[t])]][,n] <- unlist(perm_ind.plot)
				}
				if (printprogress == T){
					print(paste("regional.pop",round(t/ntr*100,2),"%"))
				}
			}
		}

		if (sum(nullmodels == "2sp.prab")>0){
			#________________________________________
			#null model regional.pop.prab: sample one populational traits values for each population within the region (among community)
			traits.nm2sp.prab <- list()
			traits_by_sp <- apply(traits, 2, function(x) tapply(x, name_sp_sites, mean, na.rm = T))

			for (t in 1: ntr){
				traits.nm2sp.prab[[eval(namestraits[t])]] <- matrix(NA, nrow = dim(traits_by_sp)[1], ncol = nperm)
				perm_ind.plot <- list()

				for(n in 1:nperm){
					perm_ind.plot <- sample(traits_by_sp[,t], dim(traits_by_sp)[1])

					traits.nm2sp.prab[[eval(namestraits[t])]][,n] <- unlist(perm_ind.plot)
				}

				if (printprogress == T){
					print(paste("regional.pop.prab",round(t/ntr*100,2),"%"))
				}
			}
		}

		########################################
		####	 calculation of random values  	####
		########################################
		null <- list()
		nm_bypop <- list()
		nm_bypop.bis <- list()

		if (printprogress == T){print("calculation of null values using null models")}

		for(i in 1:nindex){
			if (nullmodels[i] == "1"){nm.bis <- traits.nm1[[1]]}
			else if (nullmodels[i] == "2"){nm.bis <- traits.nm2[[1]]}
			else if (nullmodels[i] == "2sp"){nm.bis <- traits.nm2sp[[1]]}
			else if (nullmodels[i] == "2sp.prab"){nm.bis <- traits.nm2sp.prab[[1]]}
			else{print("nullmodels need values local, regional.ind, regional.pop, regional.pop.prab")}

			functionindex = eval(index[i])

			if (nullmodels[i] == "2sp"){
				nm_bypop.bis[[eval(namesindex[i])]] <- apply(nm.bis, 2 , function (x) tapply(x, name_sp_sites, mean , na.rm = T))

				dim2 <- dim(apply(nm_bypop.bis[[eval(namesindex[i])]], 2, function (x) eval(parse(text = functionindex))))[1]
				null[[eval(namesindex[i])]] <- array(NA, dim = c(ntr, dim2, nperm) )
				if (is.null(dim2)) {
					null[[eval(namesindex[i])]] <- array(NA, dim = c(ntr, 1, nperm) )
				}
			}

			else if (nullmodels[i] == "2sp.prab"){
				nm_bypop.bis[[eval(namesindex[i])]] <- apply(nm.bis, 2 , function (x) tapply(x, rownames(traits_by_sp), mean , na.rm = T))

				dim2 <- dim(apply(nm_bypop.bis[[eval(namesindex[i])]], 2, function (x) eval(parse(text = functionindex))))[1]
				null[[eval(namesindex[i])]] <- array(NA, dim = c(ntr, dim2, nperm) )
				if (is.null(dim2)) {
					null[[eval(namesindex[i])]] <- array(NA, dim = c(ntr, 1, nperm) )
				}
			}

			else{
				dim2 <- dim(apply(nm.bis, 2, function (x) eval(parse(text = functionindex))))[1]
				null[[eval(namesindex[i])]] <- array(NA, dim = c(ntr, dim2, nperm) )
				if (is.null(dim2)) {
					null[[eval(namesindex[i])]] <- array(NA, dim = c(ntr, 1, nperm) )
				}
			}

			for (t in 1: ntr){

				if (nullmodels[i] == "1"){nm <- traits.nm1[[t]]}
				else if (nullmodels[i] == "2"){nm <- traits.nm2[[t]]}
				else if (nullmodels[i] == "2sp"){nm <- traits.nm2sp[[t]]}
				else if (nullmodels[i] == "2sp.prab"){nm <- traits.nm2sp.prab[[t]]}
				else{print("nullmodels need values local, regional.ind, regional.pop, regional.pop.prab")}

				if (nullmodels[i] == "2sp"){
					nm_bypop[[eval(namesindex[i])]] <- apply(nm, 2 , function (x) tapply(x, name_sp_sites, mean , na.rm = T))

					null[[eval(namesindex[i])]] [t,, ] <- apply(nm_bypop[[eval(namesindex[i])]], 2, function (x) eval(parse(text = functionindex)))
				}

				else if (nullmodels[i] == "2sp.prab"){
					nm_bypop[[eval(namesindex[i])]] <- apply(nm, 2 , function (x) tapply(x, rownames(traits_by_sp), mean , na.rm = T))

					null[[eval(namesindex[i])]] [t,, ] <- apply(nm_bypop[[eval(namesindex[i])]], 2, function (x) eval(parse(text = functionindex)))
				}

				else{
					null[[eval(namesindex[i])]] [t,, ] <- apply(nm, 2, function (x) eval(parse(text = functionindex)))
				}

				if (printprogress == T){
					print(paste(eval(namesindex[i]), round(t/ntr*100,2),"%"))
				}
			}
		}
	}

	########################################
	####	calculation of observed values	####
	########################################
	obs <- list()
	traits_by_pop <- c()

	if (printprogress == T){print("calculation of observed values")}

	for(i in 1:nindex){
		functionindex = eval(index[i])

		if (is.null(nullmodels)) {
			obs[[eval(namesindex[i])]] <- array(dim = c(ntr, dim(apply(traits, 2, function (x) eval(parse(text = functionindex))))[1]))
			obs[[eval(namesindex[i])]] <- apply(traits, 2, function (x) eval(parse(text = functionindex)))
			#obs[[eval(namesindex[i])]] [ !is.finite(obs[[eval(namesindex[i])]] )] <- NA
		}

		else if(!is.null(nullmodels)){
			if (nullmodels[i] == "2sp") {
				traits_by_pop <- apply(traits, 2 , function (x) tapply(x, name_sp_sites, mean , na.rm = T))
				obs[[eval(namesindex[i])]] <- array(dim = c(ntr, dim(apply(traits_by_pop, 2, function (x) eval(parse(text = functionindex))))[1]))
				obs[[eval(namesindex[i])]] <- apply(traits_by_pop, 2, function (x) eval(parse(text = functionindex)))
			}

			if (nullmodels[i] == "2sp.prab") {
				traits_by_sp <- apply(traits, 2, function(x) tapply(x, name_sp_sites, mean, na.rm = T))
				obs[[eval(namesindex[i])]] <- array(dim = c(ntr, dim(apply(traits_by_sp, 2, function (x) eval(parse(text = functionindex))))[1]))
				obs[[eval(namesindex[i])]] <- apply(traits_by_sp, 2, function (x) eval(parse(text = functionindex)))
			}

			else if (nullmodels[i] == "1" | nullmodels[i] == "2") {
				obs[[eval(namesindex[i])]] <- array(dim = c(ntr, dim(apply(traits, 2, function (x) eval(parse(text = functionindex))))[1]))
				obs[[eval(namesindex[i])]] <- apply(traits, 2, function (x) eval(parse(text = functionindex)))
				#obs[[eval(namesindex[i])]] [ !is.finite(obs[[eval(namesindex[i])]] )] <- NA
			}
		}

		if (printprogress == T){
			print(paste(round(i/nindex*100,2),"%"))
		}
	}


	########################################
	####		Create results list		####
	########################################

	ComIndex <- list()
	ComIndex$obs <- obs

	if (is.numeric(nperm)){
		ComIndex$null <- null
	}

	ComIndex$list.index <- list()
	ComIndex$list.index.t <- list()
	name.ComIndex_list.index <- vector()

	for(i in 1:nindex){
		ComIndex$list.index.t[[seq(1,nindex*2,by = 2)[i]]] <- t(obs[[i]])
		ComIndex$list.index[[seq(1,nindex*2,by = 2)[i]]] <- obs[[i]]
		name.ComIndex_list.index[seq(1,nindex*2,by = 2)[i]] <- names(obs)[i]

		if (is.numeric(nperm)){
			ComIndex$list.index[[seq(1,nindex*2,by = 2)[i]+1]] <- null[[i]]
			ComIndex$list.index.t[[seq(1,nindex*2,by = 2)[i]+1]] <- null[[i]]
			name.ComIndex_list.index[seq(1,nindex*2,by = 2)[i]+1] <- paste(names(null)[i], "nm", sep = "_")
		}
	}

	names(ComIndex$list.index.t) <- name.ComIndex_list.index
	names(ComIndex$list.index) <- name.ComIndex_list.index

	ComIndex$sites_richness <- S
	ComIndex$namestraits <- namestraits

	ComIndex$traits <- traits
 	ComIndex$ind.plot <- ind.plot
 	ComIndex$sp <- sp

 	if (is.numeric(nperm)){
 		ComIndex$nullmodels <- nullmodels_names
		names(ComIndex$nullmodels) <- namesindex
 	}

 	ComIndex$call <- match.call()

 	class(ComIndex) <- "ComIndex"

  invisible(ComIndex)
}

ComIndexMulti <- function(traits = NULL, index = NULL, by.factor = NULL, nullmodels = NULL, ind.plot = NULL, sp = NULL, com = NULL, SE = 0, namesindex = NULL, reg.pool = NULL, SE.reg.pool = NULL, nperm = 99, printprogress = TRUE, type.sp.val = "count"){
	
	if (is.null(index)) { stop("There is no default index to compute. Please add metrics using the 'index' argument.") }

	#___________________________
	#Create/verify/sort variables to put into the function
	type.sp.val <- match.arg(type.sp.val, c("count", "abundance"))
	if (is.null(namesindex)) { namesindex <- index }
	if (!is.matrix(traits)) { traits <- as.matrix(traits) }
	ntr <- dim(traits)[2]
	namestraits <- colnames(traits)
	nindex <- length(index)
	#___________________________

	#___________________________
	#Null models argument preparation
	if (is.null(nullmodels)) {
		nperm <- NULL
		warning("nullmodels is NULL so no null model will be computed")
	}

	if (is.null(nperm)) {
		nullmodels <- NULL
		warning("nperm is NULL so no null model will be computed")
	}

	if(!is.null(nperm)){
		if (nperm == 0) {
			nperm <- NULL
			nullmodels <- NULL
			warning("nperm = 0 or nperm is NULL, so no null model will be computed")
		}
	}

	if (!is.null(nullmodels) & !is.null(nperm)){

		if (length(nullmodels) == 1){
		nullmodels <- rep(nullmodels,times = nindex)
		}

		nullmodels[nullmodels == "1"] <- "local"
		nullmodels[nullmodels == "2"] <- "regional.ind"
		nullmodels[nullmodels == "2sp"] <- "regional.pop"
		nullmodels[nullmodels == "2sp.prab"] <- "regional.pop.prab"

		nullmodels <- match.arg(nullmodels, c("local", "regional.ind", "regional.pop", "regional.pop.prab"), several.ok = TRUE)
		nullmodels_names <- nullmodels
	}
	#___________________________


	#___________________________
	if (!is.null(ind.plot) & !is.null(com)){
		warning("If ind.plot and com are provide, the function use only the argument ind.plot")
	}

	#regional pool
	if (!is.null(reg.pool) & sum(nullmodels == "regional.ind") == 0 & is.numeric(nperm)) {
		warning("Custom regional pool 'reg.pool' is only used in the case of null model regional.ind")
	}

	#If data are from species or population traits, this function (AbToInd) transform this data in a suitable format for cati
	if (is.null(ind.plot)){
		if (is.null(com)) {stop("Need 'ind.plot' or 'com' argument")}

		rownames(traits) <- sp
		res.interm <- AbToInd(traits, com, type.sp.val = type.sp.val)

		traits <- res.interm$traits
		sp <- res.interm$sp
		ind.plot <- res.interm$ind.plot
	}
	#___________________________

	#___________________________
	traits <- traits[order(ind.plot),]

	ind.plot <- ind.plot[order(ind.plot)]
	sp <- sp[order(ind.plot)]

	name_sp_sites <- paste(sp, ind.plot, sep = "_")
	comm <- NULL
	comm <- t(table(ind.plot,1:length(ind.plot)))

	S <- colSums(comm>0)
	ncom <- length(S)

	names_sp_ind.plot <- as.factor(paste(sp, ind.plot, sep = "@"))
	#___________________________

	#___________________________
	if (sum(is.na(traits))>0) {
		warning(paste("This function exclude", sum(is.na(traits)),"Na value", sep=" "))
	}

	if(length(SE) == 1){
			SE <- rep(SE, times = ncol(traits))
	}
	#___________________________

	if (is.null(by.factor)) {by.factor = rep("Region", length(ind.plot))}


	###############################################################
	if (is.numeric(nperm)){

		#########################################
		#### 	 Creating regional pools 	 ####
		#########################################

		if(!is.null(reg.pool)){isnullregpool <- FALSE} else{isnullregpool <- TRUE}

		# If the regional pool is the same for all communitiers:
		# creating a list of regional pool (one regional pool by community)
		if (is.null(reg.pool)) {
			reg.pool <- rep(list(traits), nlevels(ind.plot))
		}

		if (is.data.frame(reg.pool) | is.matrix(reg.pool) ) {
			reg.pool <- rep(list(reg.pool), nlevels(ind.plot))
		}

		# Standard error for regional pool
		if(is.null(SE.reg.pool)){
			SE.reg.pool <- SE
		}

		if(length(SE.reg.pool) == 1){
			SE.reg.pool <- rep(SE.reg.pool, times = ncol(traits))
		}

		if (is.vector(SE.reg.pool)) {
			SE.reg.pool <- rep(list(SE.reg.pool), nlevels(ind.plot))
		}

		if(length(reg.pool) != nlevels(ind.plot)){
			stop("reg.pool need to be either a matrix or a list of length equal to the number of communities")
		}

		#________________________________________
		# Warnings and stop about standard errors parameter
		if (length(SE) != ncol(traits) & length(SE) != 1) {
			stop("The vector SE need to have a length of one or equal to the number of traits")
		}

		if(!isnullregpool){
			if (length(SE.reg.pool) != length(reg.pool)) {
				stop("The vector SE.reg.pool need to have the same dimension as reg.pool")
			}
		}

		nullmodels[nullmodels == "local"] <- "1"
		nullmodels[nullmodels == "regional.ind"] <- "2"
		nullmodels[nullmodels == "regional.pop"] <- "2sp"
		nullmodels[nullmodels == "regional.pop.prab"] <- "2sp.prab"

		#########################################
		####  calculation of null models 	 ####
		#########################################
		#Creation of three null models
		if (printprogress == T){ print("creating null models")}

		if (sum(nullmodels == "1")>0){
			#________________________________________
			#null model local: sample individuals traits values within community

			traits.nm1 <- list()

			for (t in 1: ntr){
				traits.nm1[[eval(namestraits[t])]] <- matrix(NA, nrow = dim(traits)[1], ncol = nperm)
				perm_ind.plot <- list()

				for(n in 1:nperm){
					for(s in 1: ncom) {
						perm_ind.plot[[s]] <- sample(traits[ind.plot == levels(ind.plot)[s], t], table(ind.plot)[s])
					}

					traits.nm1[[eval(namestraits[t])]][,n] <- unlist(perm_ind.plot)
				}
				if (printprogress == T){
					print(paste("local",round(t/ntr*100,2),"%"))
				}
			}
		}


		if (sum(nullmodels == "2")>0){
			#________________________________________
			#null model regional.ind: sample individuals traits values within the region (among community)
			traits.nm2 <- list()

			for (t in 1: ntr){
				traits.nm2[[eval(namestraits[t])]] <- matrix(NA, nrow = dim(traits)[1], ncol = nperm)
				perm_ind.plot <- list()

				for(s in 1: ncom) {
				# Take measurement standard error into account
					trait.intern <- reg.pool[[s]][, t]
					SE.reg.pool.intern <- SE.reg.pool[[s]][t]
					if (SE.reg.pool.intern != 0) {
						trait.intern <- rnorm(length(trait.intern), mean = trait.intern, sd = SE.reg.pool.intern)
					}
				}

				for(n in 1:nperm){
					for(s in 1: ncom) {
						perm_ind.plot[[s]] <- sample(trait.intern, table(ind.plot)[s])
					}

					traits.nm2[[eval(namestraits[t])]][,n] <- unlist(perm_ind.plot)
				}
				if (printprogress == T){
					print(paste("regional.ind",round(t/ntr*100,2),"%"))
				}
			}
		}


		if (sum(nullmodels == "2sp")>0){
			#________________________________________
			#null model regional.pop: sample populational traits values within the region (among community)
			traits.nm2sp <- list()
			traits_by_sp <- apply(traits,2,function(x) tapply(x,name_sp_sites,mean, na.rm = T))

			traits_by_pop <- traits_by_sp[match(name_sp_sites,rownames(traits_by_sp)), ]

			for (t in 1: ntr){
				traits.nm2sp[[eval(namestraits[t])]] <- matrix(NA, nrow = dim(traits)[1], ncol = nperm)
				perm_ind.plot <- list()

				for(n in 1:nperm){
					for(s in 1: ncom) {
						perm_ind.plot[[s]] <- sample(traits_by_pop, table(ind.plot)[s])
					}

					traits.nm2sp[[eval(namestraits[t])]][,n] <- unlist(perm_ind.plot)
				}
				if (printprogress == T){
					print(paste("regional.pop",round(t/ntr*100,2),"%"))
				}
			}
		}


		if (sum(nullmodels == "2sp.prab")>0){
			#________________________________________
			#null model regional.pop.prab: sample one populational traits values for each population within the region (among community)
			traits.nm2sp.prab <- list()
			traits_by_sp <- apply(traits, 2, function(x) tapply(x, name_sp_sites, mean, na.rm = T))

			for (t in 1: ntr){
				traits.nm2sp.prab[[eval(namestraits[t])]] <- matrix(NA, nrow = dim(traits_by_sp)[1], ncol = nperm)
				perm_ind.plot <- list()

				for(n in 1:nperm){
					perm_ind.plot <- sample(traits_by_sp[,t], dim(traits_by_sp)[1])

					traits.nm2sp.prab[[eval(namestraits[t])]][,n] <- unlist(perm_ind.plot)
				}

				if (printprogress == T){
					print(paste("regional.pop.prab",round(t/ntr*100,2),"%"))
				}
			}
		}


		########################################
		#### calculation of random values  ####
		########################################
		null <- list()

		if (printprogress == T) {print("calculation of null values using null models")}

		for(i in 1:nindex){

			if (nullmodels[i] == "1"){nm <- array(unlist(traits.nm1),dim = c(ncol(traits), dim(traits)[1], nperm) )}
			else if (nullmodels[i] == "2"){nm <- array(unlist(traits.nm2),dim = c(ncol(traits), dim(traits)[1], nperm) )}
			else if (nullmodels[i] == "2sp"){nm <- array(unlist(traits.nm2sp),dim = c(ncol(traits), dim(traits)[1], nperm) )}
			else if (nullmodels[i] == "2sp.prab"){nm <- array(unlist(traits.nm2sp.prab),dim = c(ncol(traits), dim(traits_by_sp)[1], nperm) )}
			else{print("nullmodels need 1, 2, 2sp or 2sp.prab")}

			nm_n <- nm[,,n]
			colnames(nm_n) <- rownames(comm)
			rownames(nm_n) <- colnames(traits)

			functionindex = eval(index[i])

			dim2 <- dim(by(t(nm_n), by.factor, function (x) eval(parse(text = functionindex))))[1]
			null[[eval(namesindex[i])]] <- array(NA, dim = c(dim2, nperm) )

			if (is.null(dim2)) {
				null[[eval(namesindex[i])]] <- array(NA, dim = c(1, nperm) )
			}

			for(n in 1:nperm){
				null[[eval(namesindex[i])]][,n] <- as.vector(by(t(nm[,,n]), by.factor, function (x) eval(parse(text = functionindex))))
			}

			if (printprogress == T){
				print(paste(eval(namesindex[i]), round(i/nindex*100,2),"%"))
			}
		}

	}

	########################################
	####	calculation of observed values	####
	########################################
	obs <- list()

	if (printprogress == T){print("calculation of observed values")}

	for(i in 1:nindex){

		functionindex = eval(index[i])

		dim2 <- dim(by(traits, by.factor, function (x) eval(parse(text = functionindex))))[1]
		obs[[eval(namesindex[i])]] <- rep(NA, times = dim2)

		if (is.null(nullmodels)) {
			obs[[eval(namesindex[i])]] <- as.vector(by(traits, by.factor, function (x) eval(parse(text = functionindex))))
		}
		else if (!is.null(nullmodels)) {

			if (nullmodels[i] == "2sp") {
				traits_by_sp <- apply(traits, 2, function(x) tapply(x, name_sp_sites, mean, na.rm = T))
				traits_by_pop <- traits_by_sp[match(name_sp_sites, rownames(traits_by_sp)), ]


				obs[[eval(namesindex[i])]] <- as.vector(by(traits_by_pop, by.factor, function (x) eval(parse(text = functionindex))))
			}

			if (nullmodels[i] == "2sp.prab") {
				traits_by_sp <- apply(traits, 2, function(x) tapply(x, name_sp_sites, mean, na.rm = T))
				obs[[eval(namesindex[i])]] <- as.vector(by(traits_by_sp, by.factor, function (x) eval(parse(text = functionindex))))
			}

			else if (nullmodels[i] == "1" | nullmodels[i] == "2") {
				obs[[eval(namesindex[i])]] <- as.vector(by(traits, by.factor, function (x) eval(parse(text = functionindex))))
			}
		}

		obs[[eval(namesindex[i])]] <- as.vector(obs[[eval(namesindex[i])]])

		if (printprogress == T){
			print(paste(round(i/nindex*100,2),"%"))
		}
	}


	########################################
	####		Create results list		####
	########################################

	ComIndex <- list()
	ComIndex$obs <- obs

	if (is.numeric(nperm)){
		ComIndex$null <- null
	}

	ComIndex$list.index <- list()
	ComIndex$list.index.t <- list()
	name.ComIndex_list.index <- vector()

	for(i in 1:nindex){
		ComIndex$list.index.t[[seq(1,nindex*2,by = 2)[i]]] <- t(obs[[i]])
		ComIndex$list.index[[seq(1,nindex*2,by = 2)[i]]] <- obs[[i]]
		name.ComIndex_list.index[seq(1,nindex*2,by = 2)[i]] <- names(obs)[i]

		if (is.numeric(nperm)){
			ComIndex$list.index[[seq(1,nindex*2,by = 2)[i]+1]] <- null[[i]]
			ComIndex$list.index.t[[seq(1,nindex*2,by = 2)[i]+1]] <- null[[i]]
			name.ComIndex_list.index[seq(1,nindex*2,by = 2)[i]+1] <- paste(names(null)[i], "nm", sep = "_")
		}
	}

	names(ComIndex$list.index.t) <- name.ComIndex_list.index
	names(ComIndex$list.index) <- name.ComIndex_list.index

	ComIndex$sites_richness <- S
	ComIndex$namestraits <- namestraits

	ComIndex$traits <- traits
 	ComIndex$ind.plot <- ind.plot
 	ComIndex$sp <- sp

 	if (is.numeric(nperm)){
 		ComIndex$nullmodels <- nullmodels_names
		names(ComIndex$nullmodels) <- namesindex
 	}

 	ComIndex$call <- match.call()

	class(ComIndex) <- "ComIndexMulti"

	return(ComIndex)
}

as.listofindex <- function(x, namesindex = NULL) {

	if (class(x) != "list"){x <- list(x)}

	nlist <- length(x)
	nindex <- vector()
	for(i in 1: nlist){
		if (inherits(x[[i]], "Tstats")) {
			nindex[i] <- 3
		}
		else if (inherits(x[[i]], "ComIndex") | inherits(x[[i]], "ComIndexMulti")) {
			nindex[i] <- length(x[[i]]$obs)
		}
		else{stop("x must be a list of objects of class Tstats, ComIndex or ComIndexMulti")}
	}

	res <- list()

	for(l in 1: nlist){
		if (inherits(x[[l]], "Tstats")) {
			for(i in c(1,4,2,5,3,6) ){
				res <- c(res, list(x[[l]][[1]][[i]]))
			}
		}

		else if (inherits(x[[l]], c("ComIndex", "ComIndexMulti"))) {
			for(i in 1: nindex[l]){
				res <- c(res, list(x[[l]]$obs[[i]]), list(x[[l]]$null[[i]]) )
			}
		}
		else{stop("x must be a list of objects of class Tstats, ComIndex or ComIndexMulti")}
	}

	if (is.null(namesindex)) {
		for(l in 1: nlist){
			namesindex <- c(namesindex, paste( rep("index", nindex[l]), l, 1:nindex[l], sep = "_") )
		}
	}

	if (length(namesindex) == sum(nindex)) {
		interm <- c()

		for(i in seq(1, sum(nindex))){
			interm <- c(interm, i, i+sum(nindex))
		}

		namesindex <- c(namesindex, paste(namesindex, "nm"))[interm]
	}

	names(res) <- namesindex

	class(res) <- "listofindex"
	return(res)
}

### Function to represent standardised effect size of all indices using null models
# index.list is a list of index associate with a list of corresponding null models in this order: [1] index 1 obs - [2] index 1 null model - [3] index 2 obs - [4] index 2 null model ...
#e.g. index.list <- list(T_IP.IC = res.finch$T_IP.IC, T_IP.IC_nm = res.finch$T_IP.IC_nm, T_PC.PR = res.finch$T_PC.PR, T_PC.PR_nm = res.finch$T_PC.PR_nm)
#observed matrix of values need to be of the same dimension
#You can transpose the observed matrix to represent either the ses by traits or by plots
plot.listofindex <- function(x, type = "normal", col.index = c("red","purple","olivedrab3"), add.conf = TRUE, color.cond = TRUE, val.quant = c(0.025,0.975), grid.v = TRUE, grid.h = TRUE, xlim = NULL, ylim = NULL, cex.text = 0.8, plot.ask = FALSE, srt.text = 90, alpha = 0.4, ...){

	type <- match.arg(type, c("simple", "simple_range", "normal", "barplot", "bysites", "bytraits"))

	if (!inherits(x, "listofindex")) {
		if (inherits(x[[1]], "Tstats") | inherits(x[[2]], "ComIndex") | inherits(x[[3]], "ComIndexMulti")) {
			x <- as.listofindex(x)
		}
		else{stop("x must be a list of objects of class Tstats, ComIndex, ComIndexMulti or listofindex")}
	}

	#function from adegenet (function transp)
	transpa<-function (col, alpha = 0.5) {
		res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
		return(res)
	}

	index.list <- x

	#Graphical parameters
	par(ask = plot.ask)
	par(mar = c(5, 7, 4, 2))

	namesindex.all <- names(index.list)
	nindex <- length(names(index.list))/2
	namesindex <- names(index.list)[seq(1,nindex*2, by = 2)]
	namestraits <- colnames(index.list[[1]])
	namescommunity <- rownames(index.list[[1]])

	ncom <- c()
	ntr <- c()
	for(i in seq(1, 2*nindex, by = 2)){
		ncom <- c(ncom,dim(as.matrix(index.list[[i]]))[1])
		ntr <- c(ntr,dim(as.matrix(index.list[[i]]))[2])
	}

	if (is.null(namescommunity)) {warning("rownames of index.list[[1]] is empty so names of plots cannot be plot")}
	if (is.null(namestraits)) {warning("colnames of index.list[[1]] is empty so names of traits cannot be plot")}

	if (is.null(ncom)) {ncom <- dim(as.matrix(index.list[[1]]))[1]}
	if (is.null(ntr)) {ntr <- dim(as.matrix(index.list[[1]]))[2]}

	if (is.null(ncom)) {ncom = 1}
	if (is.null(ntr)) {ntr = 1}

	if (length(col.index)<nindex){
		col.index <- funky.col(nindex)
	}

	#________________________________________
	#calculation of standardised effect size

	res <- list()
	for (i in seq(1,nindex*2, by = 2)){
		res[[eval(namesindex.all[i])]] <- ses(obs = index.list[[i]], nullmodel = index.list[[i+1]], val.quant = val.quant)
	}

	res <- lapply(res, function(x) lapply(x, as.matrix)  )

	nfactor <- c()
	for(i in 1:nindex){
		if (ntr[i]>1) nfactor <- c(nfactor, dim(as.matrix(res[[eval(namesindex[i])]]$ses))[1])
		if (ntr[i] == 1) nfactor <- c(nfactor, length(res[[eval(namesindex[i])]]$ses))
	}

	bysite <- FALSE
	if(type == "bysites"){
		type <- "bytraits"
		bysite <- TRUE
	}

	#________________________________________
	#plot : possible type = "simple", "simple_range", "normal", "barplot" and "bytraits"

	#__________
	if (type == "bytraits"){
		par(mar = c(5, 4, 4, 2) + 0.1)

		if (is.null(ylim)) { ylim = c(0,nindex+1) }
		res.total <- unlist(res)
		res.total <- res.total[is.finite(res.total)]
		if (is.null(xlim)) {xlim = c(min(c(res.total,-2), na.rm = T), max(c(res.total,2), na.rm = T))}

		if (is.logical(color.cond)) {color.cond = c("blue","orange")}

		if (bysite){

			if (length(unique(ncom)) != 1){stop("if type = 'bysites', all index need the same number of community")}

			for (s in 1: ncom[1]){
				plot(mean(res[[eval(namesindex.all[1])]]$ses[s, ], na.rm = T),(1:nindex)[1] ,bty = "n", cex.lab = 0.8, yaxt = "n", xlab = paste("SES", namescommunity[s]), ylim = ylim, xlim = xlim, pch = 16, type = "n", ...)
				abline(v = 0)

				for(i in 1:nindex){
					abline(h = (1:nindex)[i], lty = 2, col = "lightgray")

					segments(mean(res[[eval(namesindex[i])]]$ses.sup[s, ], na.rm = T), (1:nindex)[i], mean(res[[eval(namesindex[i])]]$ses.inf[s, ], na.rm = T), (1:nindex)[i])

					points(mean(res[[eval(namesindex[i])]]$ses.sup[s, ], na.rm = T), (1:nindex)[i], pch = "|")
					points(mean(res[[eval(namesindex[i])]]$ses.inf[s, ], na.rm = T), (1:nindex)[i], pch = "|")

					points(res[[eval(namesindex[i])]]$ses[s, ], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[s, ]) ), pch = "*")

					cond.sup <- res[[eval(namesindex[i])]]$ses[s, ]>res[[eval(namesindex[i])]]$ses.sup[s, ]
					points(res[[eval(namesindex[i])]]$ses[s, ][cond.sup], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[s, ][cond.sup]) ), pch = "*", cex = 3, col = color.cond[2])

					cond.inf <- res[[eval(namesindex[i])]]$ses[s, ]<res[[eval(namesindex[i])]]$ses.inf[s, ]
					points(res[[eval(namesindex[i])]]$ses[s, ][cond.inf], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[s, ][cond.inf]) ), pch = "*", cex = 3, col = color.cond[1])

					cond.sup <- res[[eval(namesindex[i])]]$ses[s,]>res[[eval(namesindex[i])]]$ses.sup[s,]
					points(res[[eval(namesindex[i])]]$ses[s,][cond.sup], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[s,][cond.sup]) ), pch = "*", cex = 3, col = color.cond[2])

					cond.inf <- res[[eval(namesindex[i])]]$ses[s,]<res[[eval(namesindex[i])]]$ses.inf[s,]
					points(res[[eval(namesindex[i])]]$ses[s,][cond.inf], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[s,][cond.inf]) ), pch = "*", cex = 3, col = color.cond[1])

					points(mean(res[[eval(namesindex[i])]]$ses[s, ], na.rm = T), (1:nindex)[i], col = "red", pch = 16)

					text(1, (1:nindex)[i]+0.3, namesindex[i], cex = cex.text, pos = 4, font = 2)

					chh <- par()$cxy[ 2 ] ## character height
					text(res[[eval(namesindex[i])]]$ses[s, ], chh + rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[s, ]) ), namestraits, cex = cex.text, srt = srt.text,)

				}
			}
		}


		else {
			if (length(unique(ntr)) != 1){stop("if type = 'bysites', all index need the same number of community")}

			for (t in 1: ntr[1]){

				plot(mean(res[[eval(namesindex[i])]]$ses[,t], na.rm = T), (1:nindex)[i] ,bty = "n", cex.lab = 0.8, yaxt = "n", xlab = paste("SES", namestraits[t]), ylim = ylim, xlim = xlim, pch = 16, type = "n", ...)
				abline(v = 0)

				for(i in 1:nindex){

					abline(h = (1:nindex)[i], lty = 2, col = "lightgray")

					segments(mean(res[[eval(namesindex[i])]]$ses.sup[,t], na.rm = T), (1:nindex)[i], mean(res[[eval(namesindex[i])]]$ses.inf[,t], na.rm = T), (1:nindex)[i])

					points(mean(res[[eval(namesindex[i])]]$ses.sup[,t], na.rm = T), (1:nindex)[i], pch = "|")
					points(mean(res[[eval(namesindex[i])]]$ses.inf[,t], na.rm = T), (1:nindex)[i], pch = "|")

					points(res[[eval(namesindex[i])]]$ses[,t], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[,t]) ), pch = "*")

					cond.sup <- res[[eval(namesindex[i])]]$ses[,t]>res[[eval(namesindex[i])]]$ses.sup[,t]
					points(res[[eval(namesindex[i])]]$ses[,t][cond.sup], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[,t][cond.sup]) ), pch = "*", cex = 3, col = color.cond[2])

					cond.inf <- res[[eval(namesindex[i])]]$ses[,t]<res[[eval(namesindex[i])]]$ses.inf[,t]
					points(res[[eval(namesindex[i])]]$ses[,t][cond.inf], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[,t][cond.inf]) ), pch = "*", cex = 3, col = color.cond[1])


					points(mean(res[[eval(namesindex[i])]]$ses[,t], na.rm = T), (1:nindex)[i], col = "red", pch = 16)

					text(1, (1:nindex)[i]+0.3, namesindex[i], cex = 0.8, pos = 4, font = 2)

					chh <- par()$cxy[ 2 ] ## character height
					text(res[[eval(namesindex[i])]]$ses[,t], chh + rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[,t]) ), namescommunity, cex = cex.text, srt = srt.text,)
				}
			}
		}
	}

	else if (type == "simple" | type == "simple_range" | type == "normal" | type == "barplot"){

		if (is.null(ylim)) { ylim = c(0,5.5+(nindex+1)*max(ntr)) }
		res.total <- unlist(res)
		res.total <- res.total[is.finite(res.total)]
		if (is.null(xlim)) {xlim = c(min(c(res.total,-2), na.rm = T), max(c(res.total,2), na.rm = T))}


		plot(0, ylab = "Traits", yaxt = "n", xlab = "Standardized Effect Size", ylim = ylim, xlim = xlim, col = "black", type = "l", ...)
		try(axis(side = 2, seq(from = 5.5, to = 4.5+(nindex+1)*max(ntr), by = nindex+1)+(nindex+1)/2, labels = namestraits, las = 1, cex.axis = 0.7 ) )
		abline(v = 0)

		if (grid.v == T) {
			range. <- max(c(res.total,2), na.rm = T)-min(c(res.total,-2), na.rm = T)

			vect.grid <- seq(min(res.total, na.rm = T),max(res.total, na.rm = T), by = round(range.,2)/9)
			for(j in vect.grid){
				abline(v = j, lty = 2, col = "lightgray")
			}
		}

		if (grid.h == T) {
			for(j in seq(5.5,5.5+(nindex+1)*max(ntr)) ){
				abline(h = j, lty = 2, col = "lightgray")
			}
		}

		#__________
		if (type == "simple"){

			for(i in 1:nindex){
				for(t in 1:ntr[i]){

					if (add.conf == T){
						ypolyg<-c( 5.1 + (nindex+1)*t-i, 5.9 + (nindex+1)*t-i, 5.9 + (nindex+1)*t-i, 5.1 + (nindex+1)*t-i, 5.1 + (nindex+1)*t-i )
						xpolyg<-c(mean(res[[eval(namesindex[i])]]$ses.inf[,t], na.rm = T), mean(res[[eval(namesindex[i])]]$ses.inf[,t], na.rm = T), mean(res[[eval(namesindex[i])]]$ses.sup[,t], na.rm = T), mean(res[[eval(namesindex[i])]]$ses.sup[,t], na.rm = T), mean(res[[eval(namesindex[i])]]$ses.inf[,t], na.rm = T) )

						polygon(xpolyg, ypolyg, col = transpa(col.index[i], alpha), border = transpa(col.index[i], alpha + (1-alpha)/2))
					}

					if (color.cond == F){
						points(res[[eval(namesindex[i])]]$ses [,t], rep(5.5+(nindex+1)*t-i, times = nfactor[i]), pch = 20, col = col.index[i])
					}

					if (color.cond == T){
						if (length(col.index) != 2*nindex) {col.index[(nindex+1):(nindex*2)] <- "grey50"}
						points(res[[eval(namesindex[i])]]$ses [,t], rep(5.5+(nindex+1)*t-i, times = nfactor[i]), pch = 20, col = col.index[nindex+i])
						condition <- res[[eval(namesindex[i])]]$ses [,t] > res[[eval(namesindex[i])]]$ses.sup [,t] |  res[[eval(namesindex[i])]]$ses [,t] < res[[eval(namesindex[i])]]$ses.inf [,t]
						condition[is.na(condition)] <- FALSE

						if (sum(condition)>0){
							points(res[[eval(namesindex[i])]]$ses [condition,t], rep(5.5+(nindex+1)*t-i, times = sum(condition)), pch = 20, col = col.index[i])
						}
					}

					points(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm = T), 5.5+(nindex+1)*t-i, pch = 17, col = col.index[i])
					abline(seq(from = 5.5, to = 4.5+(nindex+1)*ntr[i], by = nindex+1)[t],b = 0, lty = 4, lwd = 0.2)

				}
			}
		}

		#__________
		else if (type == "simple_range"){

			for(i in 1:nindex){
				for(t in 1:ntr[i]){

					if (add.conf == T){
						ypolyg<-c( 5.1 + (nindex+1)*t-i, 5.9 + (nindex+1)*t-i, 5.9 + (nindex+1)*t-i, 5.1 + (nindex+1)*t-i, 5.1 + (nindex+1)*t-i )
						xpolyg<-c(mean(res[[eval(namesindex[i])]]$ses.inf[,t], na.rm = T), mean(res[[eval(namesindex[i])]]$ses.inf[,t], na.rm = T), mean(res[[eval(namesindex[i])]]$ses.sup[,t], na.rm = T), mean(res[[eval(namesindex[i])]]$ses.sup[,t], na.rm = T), mean(res[[eval(namesindex[i])]]$ses.inf[,t], na.rm = T) )

						polygon(xpolyg, ypolyg, col = transpa(col.index[i], alpha), border = transpa(col.index[i], alpha + (1-alpha)/2))
					}

					if (color.cond == F){
						points(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm = T), 5.5+(nindex+1)*t-i, pch = 17, col = col.index[i])
					}

					if (color.cond == T){
						if (length(col.index) != 2*nindex) {col.index[(nindex+1):(nindex*2)] <- "grey50"}
						points(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm = T), 5.5+(nindex+1)*t-i, pch = 17, col = col.index[nindex+i])

						condition <- mean(res[[eval(namesindex[i])]]$ses [,t], na.rm = T) > mean(res[[eval(namesindex[i])]]$ses.sup [,t], na.rm = T) | mean( res[[eval(namesindex[i])]]$ses [,t], na.rm = T) < mean(res[[eval(namesindex[i])]]$ses.inf [,t], na.rm = T)
						if (sum(condition)>0){
							points(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm = T), 5.5+(nindex+1)*t-i, pch = 17, col = col.index[i])
						}
					}

					segments(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm = T) + sd(res[[eval(namesindex[i])]]$ses [,t],na.rm = T), 5.5+(nindex+1)*t-i , mean(res[[eval(namesindex[i])]]$ses [,t],na.rm = T) - sd(res[[eval(namesindex[i])]]$ses [,t],na.rm = T), 5.5+(nindex+1)*t-i, col = col.index[i])

					points(min(res[[eval(namesindex[i])]]$ses [,t],na.rm = T), 5.5+(nindex+1)*t-i, pch = "*", col = col.index[i])
					points(max(res[[eval(namesindex[i])]]$ses [,t],na.rm = T), 5.5+(nindex+1)*t-i, pch = "*", col = col.index[i])

					abline(seq(from = 5.5, to = 4.5+(nindex+1)*ntr[i], by = nindex+1)[t],b = 0, lty = 4, lwd = 0.2)
				}
			}
		}

		#__________
		else if (type == "normal"){

			for(i in 1:nindex){
				for(t in 1:ntr[i]){

					if (add.conf == T){
						ypolyg<-c( 5.1 + (nindex+1)*t-i, 5.9 + (nindex+1)*t-i, 5.9 + (nindex+1)*t-i, 5.1 + (nindex+1)*t-i, 5.1 + (nindex+1)*t-i )
						xpolyg<-c(mean(res[[eval(namesindex[i])]]$ses.inf[,t], na.rm = T), mean(res[[eval(namesindex[i])]]$ses.inf[,t], na.rm = T), mean(res[[eval(namesindex[i])]]$ses.sup[,t], na.rm = T), mean(res[[eval(namesindex[i])]]$ses.sup[,t], na.rm = T), mean(res[[eval(namesindex[i])]]$ses.inf[,t], na.rm = T) )

						polygon(xpolyg, ypolyg, col = transpa(col.index[i], alpha), border = transpa(col.index[i], alpha + (1-alpha)/2))
					}

					if (color.cond == F){
						points(res[[eval(namesindex[i])]]$ses [,t], rep(5.5+(nindex+1)*t-i, times = nfactor[i]), pch = 20, col = col.index[i])
					}

					if (color.cond == T){
						if (length(col.index) != 2*nindex) {col.index[(nindex+1):(nindex*2)] <- "grey50"}
						points(res[[eval(namesindex[i])]]$ses [,t], rep(5.5+(nindex+1)*t-i, times = nfactor[i]), pch = 20, col = col.index[nindex+i])
						condition <- res[[eval(namesindex[i])]]$ses [,t] > res[[eval(namesindex[i])]]$ses.sup [,t] |  res[[eval(namesindex[i])]]$ses [,t] < res[[eval(namesindex[i])]]$ses.inf [,t]
						condition[is.na(condition)] <- FALSE

						if (sum(condition)>0){
							points(res[[eval(namesindex[i])]]$ses [condition,t], rep(5.5+(nindex+1)*t-i, times = sum(condition)), pch = 20, col = col.index[i])
						}
					}
					segments(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm = T) + sd(res[[eval(namesindex[i])]]$ses [,t],na.rm = T), 5.5+(nindex+1)*t-i , mean(res[[eval(namesindex[i])]]$ses [,t],na.rm = T) - sd(res[[eval(namesindex[i])]]$ses [,t],na.rm = T), 5.5+(nindex+1)*t-i, col = col.index[i])
					points(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm = T), 5.5+(nindex+1)*t-i, pch = 10, lwd = 1.6, col = col.index[i])

					abline(seq(from = 5.5, to = 4.5+(nindex+1)*ntr[i], by = nindex+1)[t],b = 0, lty = 4, lwd = 0.2)
				}
			}
		}

		#__________
		else if (type == "barplot"){
			for(i in 1:nindex){
				for(t in 1:ntr[i]){

					if (add.conf == T){
						ypolyg<-c( 5.1 + (nindex+1)*t-i, 5.9 + (nindex+1)*t-i, 5.9 + (nindex+1)*t-i, 5.1 + (nindex+1)*t-i, 5.1 + (nindex+1)*t-i )
						xpolyg<-c(mean(res[[eval(namesindex[i])]]$ses.inf[,t], na.rm = T), mean(res[[eval(namesindex[i])]]$ses.inf[,t], na.rm = T), mean(res[[eval(namesindex[i])]]$ses.sup[,t], na.rm = T), mean(res[[eval(namesindex[i])]]$ses.sup[,t], na.rm = T), mean(res[[eval(namesindex[i])]]$ses.inf[,t], na.rm = T) )

						polygon(xpolyg, ypolyg, col = transpa(col.index[i], alpha), border = transpa(col.index[i], alpha + (1-alpha)/2))
					}

					segments(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm = T), 5.5+(nindex+1)*t-i, 0, 5.5+(nindex+1)*t-i, pch = 17, col = col.index[i], lwd = 8)
					segments(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm = T) + sd(res[[eval(namesindex[i])]]$ses [,t],na.rm = T), 5.5+(nindex+1)*t-i , mean(res[[eval(namesindex[i])]]$ses [,t],na.rm = T) - sd(res[[eval(namesindex[i])]]$ses [,t],na.rm = T), 5.5+(nindex+1)*t-i, col = col.index[i])

					abline(seq(from = 5.5, to = 4.5+(nindex+1)*ntr[i], by = nindex+1)[t],b = 0, lty = 4, lwd = 0.2)

				}
			}

		}
		legend("bottomleft", inset = .005, namesindex, fill = col.index, ncol = round(nindex/3) ,cex = 0.6, bty = "0")
	}

	else{print(paste("Error:",type,"is not a valid type of plot"))}

	#return to default parameter
	par(ask = FALSE)
	par(mar = c(5, 4, 4, 2) + 0.1)
}

plot.ComIndex <- function(x, type = "normal", col.index = c("red","purple","olivedrab3"), add.conf = TRUE, color.cond = TRUE, val.quant = c(0.025,0.975), ...){

	type <- match.arg(type, c("simple", "simple_range", "normal", "barplot", "bysites", "bytraits"))

	if (!inherits(x, "ComIndex")) {
		stop("x must be of class ComIndex")
	}

	plot.listofindex (as.listofindex(x), type = type, col.index = col.index, add.conf = add.conf, color.cond = color.cond, val.quant = val.quant, ...)
}

plot.ComIndexMulti <- function(x, type = "normal", col.index = c("red","purple","olivedrab3"), add.conf = TRUE, color.cond = TRUE, val.quant = c(0.025,0.975), ...){

	type <- match.arg(type, c("simple", "simple_range", "normal", "barplot", "bysites", "bytraits"))

	if (!inherits(x, "ComIndexMulti")) {
		stop("x must be of class ComIndexMulti")
	}

	plot.listofindex (as.listofindex(x), type = type, col.index = col.index, add.conf = add.conf, color.cond = color.cond, val.quant = val.quant, ...)
}


#______________#______________#______________#______________#______________#______________#______________#______________
#______________#______________#______________#______________#______________#______________#______________#______________
#__ Decomposition of variance

decompCTRE <- function(traits = NULL , formula = ~ 1, ind.plot = NULL, sp = NULL, printprogress = TRUE, ...) {

	form.string <- deparse(substitute(formula))

 	ntr <- dim(traits)[2]
	namestraits <- colnames(traits)

	comm <- t(table(ind.plot,1:length(ind.plot)))
	S <- colSums(comm>0)
	ncom <- length(S)

	decomp <- list()

	moy_pop <- apply(traits, 2, function(x) tapply(x, paste(ind.plot, sp, sep = "."), mean, na.rm = T))
	moy_sp <- apply(traits, 2, function(x) tapply(x, sp, mean, na.rm = T))



	for(t in 1:ntr){
		if(sum(is.na(traits[,t])) > 0){
			warning(paste("All individuals with one NA (", sum(is.na(traits[,t])),"individual value(s)) are excluded for the trait", t, ":", namestraits[t], sep = " "))
		}

		interm <- apply(comm, 2, function(x) tapply(x, paste(ind.plot, sp, sep = ".") , sum, na.rm = T ) )[!is.na(moy_pop[,t]), 1:ncom]
		interm2 <- apply(comm, 2, function(x) tapply(x, sp, sum, na.rm = T ) )[!is.na(moy_sp[,t]), 1:ncom]

		specif.avg <- t(moy_pop[,t][!is.na(moy_pop[,t])] %*% interm )/colSums(interm)
		const.avg <- t(moy_sp[,t][!is.na(moy_sp[,t])] %*% interm2 )/colSums(interm2)

		flex = paste("traitflex.anova(", form.string , ", specif.avg, const.avg, ...)", sep = "")

 		decomp[[eval(namestraits[t])]] <- as.vector(eval(parse(text = flex)))

		if (printprogress == T){print(paste(round(t/ntr*100, 2),"%"))} else {}
	}

	decomp$traits <- namestraits
	class(decomp) <- "decompCTRE"
	return(decomp)
}

barplot.decompCTRE <- function(height, resume = TRUE, ...) {

	x <- height

	oldpar <- par(no.readonly = TRUE)

	if (resume == TRUE){
		res <- matrix(ncol = dim(as.matrix(x))[1]-1, nrow = dim(x[[1]]$SumSq)[2]-1)
		for(i in 1: (dim(as.matrix(x))[1]-1)) {
			res[,i] <- plot.traitflex(x[[i]], plot.res = FALSE, cumul = TRUE, beside = T, ...)
		}
		colnames(res) <- x$traits
		barplot(res, beside = T)
		legend("topright", c("Turnover", "Intraspec.", "Covariation"), fill = c(gray.colors(3)) )
	}

	else if (resume == FALSE){
		par(cex = 2/length(x$traits))

		for(i in 1:dim(x[[1]]$SumSq)[2]){
			plot.traitflex(x[[i]], main = x$traits[i], ...)
			abline(v = 0, lty = 2)
		}
	}

	par(oldpar)

}

# Leps et al. function
traitflex.anova <- function(formula, specif.avg, const.avg, ...) {

  # Formula into string form
  form.string <- deparse(substitute(formula))

  # Check formula parameter
  form.parts <- unlist(strsplit(form.string,"~"))

  if (length(form.parts) != 2){
		stop("First parameter must be valid one-sided formula, like ~A*B")
  }

  if (nchar(form.parts[1])>0){
		warning("Left side of the formula was ignored!")
  }

  # The two average variables into string form
  spec.av <- deparse(substitute(specif.avg))
  cons.av <- deparse(substitute(const.avg))

  test.has.onelevel <- function(aov.summ)
  {
   (length(aov.summ) == 1)
  }

  test.has.resid <- function(aov.one)
  {
		if (class(aov.one)[1] != "anova"){
			warning("specified object is not aov result!")
		}

		nrows <- dim(aov.one)[1]

		if (nrows == 0){
			return( FALSE)
		}

		last.row.lbl <- dimnames(aov.one)[[1]][nrows]

		if (unlist(strsplit(last.row.lbl, " "))[1] != "Residuals"){
			return( FALSE)
		}

		if (dim(aov.one)[2] < 5){ # P and F are missing
			return( FALSE)
		}

		TRUE
  }

  # Specific averages ANOVA
  form.1 <- as.formula(paste(spec.av,form.parts[2],sep = "~"))
  res.1 <- summary( aov(form.1, ...))

  if (test.has.onelevel( res.1) == FALSE){
		stop("Cannot evaluate ANOVAs with multiple error levels!")
  }

  res.1 <- res.1[[1]]

  if (test.has.resid( res.1) == FALSE){
		stop("No residual DFs left, cannot continue")
  }

  # Constant averages ANOVA
  form.2 <- as.formula(paste(cons.av,form.parts[2],sep = "~"))
  # no need to test for multilevels by now
  res.2 <- summary( aov(form.2, ...))[[1]]

  if (test.has.resid( res.2) == FALSE){
		stop("No residual DFs left in constant averages ANOVA, cannot continue")
  }

  # Now the differences:
  spec.const.diff <- paste("I(", spec.av, "-", cons.av, ")", sep = "")

  form.3 <- as.formula(paste(spec.const.diff,form.parts[2],sep = "~"))
  # no need to test for multilevels by now
  res.3 <- summary( aov(form.3, ...))[[1]]

  if (test.has.resid( res.3) == FALSE){
		stop("No residual DFs left in (specif-const) ANOVA, cannot continue")
  }

  if ((dim(res.1) != dim(res.2)) || (dim(res.1) != dim(res.3))){
  stop("Tables from the three ANOVAs have incompatible sizes")
  }

  # Create sum of squares table: add SS(Tot) except for null models
  nrows <- dim(res.1)[1]
  ss.turn <- res.2[,2]
  ss.var <- res.3[,2]
  ss.tot <- res.1[,2]
  ss.covar <- ss.tot - ss.turn - ss.var
  ss.row.names <- dimnames(res.1)[[1]]
  if (nrows > 1)
  {
		ss.turn <- c(ss.turn, sum(ss.turn))
		ss.var <- c(ss.var, sum(ss.var))
		ss.tot <- c(ss.tot, sum(ss.tot))
		ss.covar <- c(ss.covar,sum(ss.covar))
		ss.row.names <- c(ss.row.names, "Total")
		nrows  <- nrows + 1
  } else {
		# replace row title
		ss.row.names[1] <- "Total"
  }

  SS.tab <- data.frame( Turnover = ss.turn, Intraspec. = ss.var,
             Covariation = ss.covar, Total = ss.tot,
             row.names = ss.row.names)
  # Calculate relative fractions
  TotalSS <- SS.tab[nrows, 4] # lower right corner
  SS.tab.rel <- SS.tab / TotalSS

  # Collect significances
  if (nrows > 1){ # get rid of the "Total" label again
		ss.row.names <- ss.row.names[-nrows]
  }

  P.tab <- data.frame( Turnover = res.2[,5], Intraspec. = res.3[,5],
             Total = res.1[,5], row.names = ss.row.names)

  res <- list( SumSq = SS.tab, RelSumSq = SS.tab.rel, Pvals = P.tab,
         anova.turnover = res.2, anova.total = res.1, anova.diff = res.3)
  class(res) <- "traitflex"
  res
 }

print.traitflex <- function(x, ...)  {
  op <- options()
  cat("\n Decomposing trait sum of squares into composition turnover\n")
  cat( " effect, intraspecific trait variability, and their covariation:\n")
  options(digits = 5)
  print(x$SumSq, ...)
  cat("\n Relative contributions:\n")
  options(digits = 4)
  print(x$RelSumSq, ...)
  nPvals <- dim(x$Pvals)[1]
  if (nPvals > 1)
  {
   cat("\n Significance of testable effects:\n")
   options(digits = 5)
   print(x$Pvals[-nPvals, ], ...)
  }
  options(op)
  invisible(x)
 }

plot.traitflex <- function(x, plot.total = FALSE, use.percentage = TRUE, plot.covar = FALSE, cumul = FALSE, legend.pos = if (plot.total) "topleft" else "topright", plot.res = TRUE, ...) {
  if (use.percentage)
   SumSq <- 100 * x$RelSumSq
  else
   SumSq <- x$SumSq
  if (legend.pos == "none")
   legend.txt <- NULL
  else
   legend.txt <- colnames(SumSq)[-4]

  nrows <- dim(SumSq)[1]
  plot.tab <- as.matrix(SumSq)

  if (nrows > 1)
  {
   if (plot.covar)
   {
    if (plot.total)
     plot.tab <- plot.tab[,-4]
    else
     plot.tab <- plot.tab[-nrows,-4]
   }
   else
   {
    if (plot.total)
     plot.tab <- plot.tab[,1:2]
    else
     plot.tab <- plot.tab[-nrows, 1:2]
    if (legend.pos != "none")
     legend.txt <- legend.txt[1:2]
   }
  }

  else
  {
		if (!plot.total){
			plot.tab <- t(plot.tab[,-4])
    }

		legend.pos <- "none"
		legend.txt <- NULL
  }

  if (plot.res == TRUE){
	  if (!cumul){

		  if (use.percentage)
		  {
		   xpos <- barplot( plot.tab, horiz = T,
		            ylab = "Explained variation (%)", ...)
		  }
		  else
		   xpos <- barplot( plot.tab = T, horiz = T,
		            ylab = "Sum of squares of analysed trait", ...)

		  if (plot.total == FALSE){
				abline(v = as.matrix(SumSq)[,4])
		  }

	  }


		if (cumul){

		  if (use.percentage)
		  {
		   xpos <- barplot( t(plot.tab), horiz = T,
		            ylab = "Explained variation (%)", ...)
		  }
		  else
		   xpos <- barplot( t(plot.tab), horiz = T,
		            ylab = "Sum of squares of analysed trait", ...)

		  if (plot.total == FALSE){
				abline(v = as.matrix(SumSq)[,4], lty = 2)
		  }

		  if (min(plot.tab)<0){
				abline(v = 0)
			}
	  }
	}


  if (nrows > 1)
  {
   if (!plot.covar)
   { if (length(xpos) > 1)
    line.half <- 0.4*(xpos[2] - xpos[1])
    else
     line.half <- 0.4*(xpos[1])
    segments( xpos-line.half, SumSq[,4],
         xpos+line.half, SumSq[,4], lwd = 3)
    if (legend.pos != "none")
    { NL <- length(legend.txt)
     legend( legend.pos, legend = c(legend.txt,"Total"),
         fill = c(gray.colors(NL), NA),
         lty = c(rep( 0,NL),1),
         lwd = c(rep( 0,NL),3))
    }
   }
   else
   {
    if (legend.pos != "none")
     legend( legend.pos, legend = legend.txt,
         fill = gray.colors(length(legend.txt)))
   }
  }

  if (!plot.res){
		return(plot.tab)
  }
}





#______________#______________#______________#______________#______________#______________#______________#______________
#______________#______________#______________#______________#______________#______________#______________#______________
#__ Other plotting functions

### Function to represent correlations between Tstats
plotCorTstats <- function(tstats = NULL, val.quant = c(0.025,0.975), add.text = FALSE, bysite = FALSE, col.obj = NULL, plot.ask = TRUE, multipanel = TRUE, ...) {

	oldpar <- par(no.readonly = TRUE)
	par(ask = plot.ask)
	Tst <- tstats$Tstats

	#________________________________________
	ses.T_IP.IC.moy <- t(colMeans((Tst$T_IP.IC-apply(Tst$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm = T)), na.rm = T))
	ses.T_IC.IR.moy <- t(colMeans((Tst$T_IC.IR-apply(Tst$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm = T)), na.rm = T))
	ses.T_PC.PR.moy <- t(colMeans((Tst$T_PC.PR-apply(Tst$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm = T)), na.rm = T))

	ses.T_IP.IC <- t((Tst$T_IP.IC-apply(Tst$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm = T)))
	ses.T_IC.IR <- t((Tst$T_IC.IR-apply(Tst$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm = T)))
	ses.T_PC.PR <- t((Tst$T_PC.PR-apply(Tst$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm = T)))

	ses.T_IP.IC.inf <- t((apply(Tst$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm = T, prob = val.quant[1]))-apply(Tst$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm = T)))
	ses.T_IC.IR.inf <- t((apply(Tst$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm = T, prob = val.quant[1]))-apply(Tst$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm = T)))
	ses.T_PC.PR.inf <- t((apply(Tst$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm = T, prob = val.quant[1]))-apply(Tst$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm = T)))

	ses.T_IP.IC.sup <- t((apply(Tst$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm = T, prob = val.quant[2]))-apply(Tst$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm = T)))
	ses.T_IC.IR.sup <- t((apply(Tst$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm = T, prob = val.quant[2]))-apply(Tst$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm = T)))
	ses.T_PC.PR.sup <- t((apply(Tst$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm = T, prob = val.quant[2]))-apply(Tst$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm = T)))/apply(Tst$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm = T)))

	cond.T_IP.IC.inf <- ses.T_IP.IC<ses.T_IP.IC.inf
	cond.T_IC.IR.inf <- ses.T_IC.IR<ses.T_IC.IR.inf
	cond.T_PC.PR.inf <- ses.T_PC.PR<ses.T_PC.PR.inf

	cond.T_IP.IC.sup <- ses.T_IP.IC>ses.T_IP.IC.sup
	cond.T_IC.IR.sup <- ses.T_IC.IR>ses.T_IC.IR.sup
	cond.T_PC.PR.sup <- ses.T_PC.PR>ses.T_PC.PR.sup

	#________________________________________
	if (bysite == F){

		if (is.null(col.obj)) {col.obj <- funky.col(dim(Tst$T_IP.IC)[2])}
		else{}
		if (multipanel) {par(mfrow = c(sqrt(dim(Tst$T_IP.IC)[2])+1,sqrt(dim(Tst$T_IP.IC)[2])+1)) }

		#__________
		#First panel of figures

		plot(0,0, xlim = c(-4,4), ylim = c(-4,4), cex.lab = 1.2 ,ylab = "ses.T_IP.IC", xlab = "ses.T_IC.IR", type = "n")
		abline(h = 2)
		abline(v = 2)
		abline(h = -2)
		abline(v = -2)
		text(0,0,"null \r\n model \r\n zone")

		for(t in 1:dim(Tst$T_IP.IC)[2]){
			plot(as.vector(ses.T_IP.IC)~as.vector(ses.T_IC.IR), col = "grey", pch = 20, main = rownames(ses.T_IC.IR)[t], ...)
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_IC.IR.inf)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type = "l")
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_IC.IR.sup)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type = "l")
			points(sort(as.vector(ses.T_IC.IR)),as.vector(ses.T_IP.IC.inf)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type = "l")
			points(sort(as.vector(ses.T_IC.IR)),as.vector(ses.T_IP.IC.sup)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type = "l")
			segments(rep(rowMeans(ses.T_IC.IR, na.rm = T)[t], times = dim(Tst$T_IP.IC)[1]), rep(rowMeans(ses.T_IP.IC, na.rm = T)[t], times = dim(Tst$T_IP.IC)[1]) ,ses.T_IC.IR[t, ], ses.T_IP.IC[t, ], col = col.obj[t])
		}

		plot.new()

		#__________
		#Second panel of figures

		plot(0,0, xlim = c(-4,4), ylim = c(-4,4), cex.lab = 1.2 ,ylab = "ses.T_IP.IC", xlab = "ses.T_PC.PR", type = "n")
		abline(h = 2)
		abline(v = 2)
		abline(h = -2)
		abline(v = -2)
		text(0,0,"null \r\n model \r\n zone")

		for(t in 1:dim(Tst$T_IP.IC)[2]){
			plot(as.vector(ses.T_IP.IC)~as.vector(ses.T_PC.PR), col = "grey", pch = 20, main = rownames(ses.T_PC.PR)[t], ...)
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_PC.PR.inf)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type = "l")
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_PC.PR.sup)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type = "l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IP.IC.inf)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type = "l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IP.IC.sup)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type = "l")
			segments(rep(rowMeans(ses.T_PC.PR, na.rm = T)[t], times = dim(Tst$T_IP.IC)[1]), rep(rowMeans(ses.T_IP.IC, na.rm = T)[t], times = dim(Tst$T_IP.IC)[1]) ,ses.T_PC.PR[t, ], ses.T_IP.IC[t, ], col = col.obj[t])
		}

		plot.new()

		#__________
		#Third panel of figures

		plot(0,0, xlim = c(-4,4), ylim = c(-4,4), cex.lab = 1.2 ,ylab = "ses.T_IC.IR", xlab = "ses.T_PC.PR", type = "n")
		abline(h = 2)
		abline(v = 2)
		abline(h = -2)
		abline(v = -2)
		text(0,0,"null \r\n model \r\n zone")

		for(t in 1:dim(Tst$T_IC.IR)[2]){
			plot(as.vector(ses.T_IC.IR)~as.vector(ses.T_PC.PR), col = "grey", pch = 20, main = rownames(ses.T_PC.PR)[t], ...)
			points(sort(as.vector(ses.T_IC.IR))~as.vector(ses.T_PC.PR.inf)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type = "l")
			points(sort(as.vector(ses.T_IC.IR))~as.vector(ses.T_PC.PR.sup)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type = "l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IC.IR.inf)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type = "l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IC.IR.sup)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type = "l")
			segments(rep(rowMeans(ses.T_PC.PR, na.rm = T)[t], times = dim(Tst$T_IC.IR)[1]), rep(rowMeans(ses.T_IC.IR, na.rm = T)[t], times = dim(Tst$T_IC.IR)[1]) ,ses.T_PC.PR[t, ], ses.T_IC.IR[t, ], col = col.obj[t])
		}
	}

	#________________________________________
	else if (bysite == T){

		if (is.null(col.obj)) {col.obj <- funky.col(dim(Tst$T_IP.IC)[1])}
		else{}

		if (multipanel) {par(mfrow = c(sqrt(dim(Tst$T_IP.IC)[1])+1,sqrt(dim(Tst$T_IP.IC)[1])+1))}

		#__________
		#First panel of figures

		plot(0,0, xlim = c(-4,4), ylim = c(-4,4), cex.lab = 1.2 ,ylab = "ses.T_IC.IR", xlab = "ses.T_IC.IR", type = "n")
		abline(h = 2)
		abline(v = 2)
		abline(h = -2)
		abline(v = -2)
		text(0,0,"null \r\n model \r\n zone")

		for(s in 1:dim(Tst$T_IP.IC)[1]){
			plot(as.vector(ses.T_IP.IC)~as.vector(ses.T_IC.IR), col = "grey", pch = 20, main = colnames(ses.T_PC.PR)[s], ...)
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_IC.IR.inf)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type = "l")
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_IC.IR.sup)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type = "l")
			points(sort(as.vector(ses.T_IC.IR)),as.vector(ses.T_IP.IC.inf)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type = "l")
			points(sort(as.vector(ses.T_IC.IR)),as.vector(ses.T_IP.IC.sup)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type = "l")
			segments(rep(colMeans(ses.T_IC.IR, na.rm = T)[s], times = dim(Tst$T_IP.IC)[2]), rep(colMeans(ses.T_IP.IC, na.rm = T)[s], times = dim(Tst$T_IP.IC)[2]) ,ses.T_IC.IR[,s], ses.T_IP.IC[,s], col = col.obj[s])
		}

		#__________
		#Second panel of figures

		plot(0,0, xlim = c(-4,4), ylim = c(-4,4), cex.lab = 1.2 ,ylab = "ses.T_IP.IC", xlab = "ses.T_PC.PR", type = "n")
		abline(h = 2)
		abline(v = 2)
		abline(h = -2)
		abline(v = -2)

		text(0,0,"null \r\n model \r\n zone")

		for(s in 1:dim(Tst$T_IP.IC)[1]){
			plot(as.vector(ses.T_IP.IC)~as.vector(ses.T_PC.PR), col = "grey", pch = 20, main = colnames(ses.T_PC.PR)[s], ... )
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_PC.PR.inf)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type = "l")
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_PC.PR.sup)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type = "l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IP.IC.inf)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type = "l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IP.IC.sup)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type = "l")
			segments(rep(colMeans(ses.T_PC.PR, na.rm = T)[s], times = dim(Tst$T_IP.IC)[2]), rep(colMeans(ses.T_IP.IC, na.rm = T)[s], times = dim(Tst$T_IP.IC)[2]) ,ses.T_PC.PR[,s], ses.T_IP.IC[,s], col = col.obj[s])
		}

		#__________
		#Third panel of figures

		plot(0,0, xlim = c(-4,4), ylim = c(-4,4), cex.lab = 1.2 ,ylab = "ses.T_IC.IR", xlab = "ses.T_PC.PR", type = "n")
		abline(h = 2)
		abline(v = 2)
		abline(h = -2)
		abline(v = -2)

		text(0,0,"null \r\n model \r\n zone")

		for(s in 1:dim(Tst$T_IC.IR)[1]){
			plot(as.vector(ses.T_IC.IR)~as.vector(ses.T_PC.PR), col = "grey", pch = 20, main = colnames(ses.T_PC.PR)[s], ... )
			points(sort(as.vector(ses.T_IC.IR))~as.vector(ses.T_PC.PR.inf)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type = "l")
			points(sort(as.vector(ses.T_IC.IR))~as.vector(ses.T_PC.PR.sup)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type = "l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IC.IR.inf)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type = "l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IC.IR.sup)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type = "l")
			segments(rep(colMeans(ses.T_PC.PR, na.rm = T)[s], times = dim(Tst$T_IC.IR)[2]), rep(colMeans(ses.T_IC.IR, na.rm = T)[s], times = dim(Tst$T_IC.IR)[2]) ,ses.T_PC.PR[,s], ses.T_IC.IR[,s], col = col.obj[s])
		}
	}

	else{print("Error: obj need to be either traits or sites")}

	par(oldpar)
}

# plot ses of an index against an other variable wich correspond to plot. For example species richness or a gradient variable

plotSESvar <- function(index.list, variable = NULL, ylab = "variable", color.traits = NULL, val.quant = c(0.025,0.975), resume = FALSE, multipanel = TRUE){

	y <- variable

	namesindex.all <- names(index.list)
	nindex <- length(names(index.list))/2
	namesindex <- names(index.list)[seq(1,nindex*2, by = 2)]
	namestraits <- colnames(index.list[[1]])
	namescommunity <- rownames(index.list[[1]])

	ncom <- dim(index.list[[1]])[1]
	ntr <- dim(index.list[[1]])[2]

	if (is.null(color.traits)){
		color.traits <- palette(terrain.colors(ntr))
	}

	#________________________________________
	#calculation of standardised effect size

	res <- list()
	for (i in seq(1,nindex*2, by = 2)){
		res[[eval(namesindex.all[i])]] <- ses(obs = index.list[[i]], nullmodel = index.list[[i+1]], val.quant = val.quant)
	}

	oldpar <- par(no.readonly = TRUE)
	if (multipanel) {par(mfrow = c(ceiling(sqrt(nindex))-1, ceiling(sqrt(nindex))))}

	ylim = c(min(y, na.rm = T), max(y, na.rm = T))


	for(i in seq(1,nindex*2, by = 2)){
		if (resume == FALSE){xlim = c(min(c(-4,res[[eval(namesindex.all[i])]]$ses), na.rm = T), max(c(4,res[[eval(namesindex.all[i])]]$ses), na.rm = T))}
		else{xlim = c(min(c(-4,rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm = T)-apply(res[[eval(namesindex.all[i])]]$ses, 1,sd,na.rm = T), rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm = T)+apply(res[[eval(namesindex.all[i])]]$ses, 1,sd,na.rm = T)), na.rm = T), max(c(4,rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm = T)-apply(res[[eval(namesindex.all[i])]]$ses, 1,sd,na.rm = T), rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm = T)+apply(res[[eval(namesindex.all[i])]]$ses, 1,sd,na.rm = T)), na.rm = T)) }

		plot(0, 0 ,bty = "n", cex.lab = 0.8, xlab = paste("SES", namesindex.all[i]), ylab = ylab, ylim = ylim, xlim = xlim, pch = 16, type = "n")
		abline(v = 0, lty = 1, col = "black")


		if (resume == TRUE){
			points(rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm = T), y, pch = 16)
			segments(rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm = T)-apply(res[[eval(namesindex.all[i])]]$ses, 1,sd,na.rm = T), y, rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm = T)+apply(res[[eval(namesindex.all[i])]]$ses, 1,sd,na.rm = T), y, pch = 16)
		}

		else{
			for(nco in 1:ncom){
				abline(h = y[nco], lty = 2, pch = 0.5, col = "lightgray")
			}

			for (t in 1: ntr){
				points(res[[eval(namesindex.all[i])]]$ses[,t] , y, pch = 16, col = color.traits[t])
			}
		}
	}

	if (resume != TRUE){
		plot(0, 0 ,bty = "n", cex.lab = 0.8, xlab = paste("SES", namesindex.all[i]), ylim = ylim, xlim = xlim, pch = 16, type = "n")
		legend("center", legend = namestraits, fill = color.traits, bty = "n", ncol = round(sqrt(nlevels(as.factor(namestraits)))-1 ) )
	}

	par(oldpar)

}


plotDistri <- function(traits = NULL, var.1 = NULL, var.2 = NULL, col.dens = NULL, plot.ask = TRUE, ylim.cex = 1, cex.leg = 0.8, polyg = TRUE, multipanel = TRUE, leg = TRUE, xlim = NULL, ylim = NULL, main = "default", ...) {

	var.1 <- as.factor(as.vector(var.1))
	var.2 <- as.factor(as.vector(var.2))

	namestraits <- colnames(traits)
	namescommunity <- unique(var.1)
	ncom <- length(namescommunity)
	ntr <- dim(as.matrix(traits))[2]

	#Graphical parameters
	if (is.null(col.dens)) {col.dens <- funky.col(nlevels(as.factor(var.2)))}
	if (length(leg) != ntr){
		leg <- rep_len(leg, ntr)
	}
	if(plot.ask | multipanel) {
		oldpar <- par(no.readonly = TRUE)
	}
	par(ask = plot.ask)
	if (multipanel) {
		par(mfrow = c(2,2))
	}

	# main
	main.plot <- c()

	if(length(main) != ntr*ncom){
		rep(main, times = ntr*ncom)
	}


	# start plotting
	for(co in 1:ncom){
		for(t in 1:ntr){
			# main
			if(!is.null(main)){
				if(main[1] == "default") {
					main.plot <- paste(namestraits[t], levels(as.factor(var.1))[co], " ")
				}

				else{
					main.plot <- main[t]
				}
			}

			if (length(na.omit(traits[as.factor(var.1) == levels(as.factor(var.1))[co],t]))>1){

				#______
				#Define the limit for the plot
				xli.interm <- c()
				yli.interm <- c()
				xlimin.interm <- c()
				ylimin.interm <- c()
				for(s in 1:nlevels(as.factor(var.2))) {
					if (length(na.omit(traits[as.factor(var.1) == levels(as.factor(var.1))[co] & as.factor(var.2) == levels(as.factor(var.2))[s],t]))>1) {
						xli.interm[s] <- max(density(traits[as.factor(var.1) == levels(as.factor(var.1))[co] & as.factor(var.2) == levels(as.factor(var.2))[s],t], na.rm = T)$x, na.rm = TRUE)
						yli.interm[s] <- max(density(traits[as.factor(var.1) == levels(as.factor(var.1))[co] & as.factor(var.2) == levels(as.factor(var.2))[s],t], na.rm = T)$y, na.rm = TRUE)
						xlimin.interm[s] <- min(density(traits[as.factor(var.1) == levels(as.factor(var.1))[co] & as.factor(var.2) == levels(as.factor(var.2))[s],t], na.rm = T)$x, na.rm = TRUE)
						ylimin.interm[s] <- min(density(traits[as.factor(var.1) == levels(as.factor(var.1))[co] & as.factor(var.2) == levels(as.factor(var.2))[s],t], na.rm = T)$y, na.rm = TRUE)
					}
				}
				if(is.null(xlim)){
					xli <- max(xli.interm, na.rm = TRUE)
					xlimin <- min(xlimin.interm, na.rm = TRUE)
					xlimite = c(min(c(min(density(traits[,t], na.rm = T)$x), xlimin)), max(c(max(density(traits[,t], na.rm = T)$x, na.rm = T), xli)))
				}
				else{xlimite <- xlim}

				if(is.null(ylim)){
					yli <- max(yli.interm, na.rm = TRUE)
					ylimin <- min(ylimin.interm, na.rm = TRUE)
					ylimite = c(min(c(min(density(traits[,t], na.rm = T)$y), ylimin)), max(c(max(density(traits[,t], na.rm = T)$y, na.rm = T), yli))*1.05)
				}
				else{ylimite <- ylim}

				#______
				#Plot the regional distribution
				plot(main = main.plot, density(traits[as.factor(var.1) == levels(as.factor(var.1))[co],t], na.rm = T), ylim = ylimite, xlim = xlimite, col = "black", ...)


				lines(density(traits[,t], na.rm = T), lty = 2, col = "grey")

				if (polyg == T) {
					x <- density(traits[as.factor(var.1) == levels(as.factor(var.1))[co],t], na.rm = T)$x
					y <- density(traits[as.factor(var.1) == levels(as.factor(var.1))[co],t], na.rm = T)$y
					polygon(c(x,rev(x)), c(rep(0,length(x)),rev(y)), border = NA, col = rgb(0.5,0.5,0.5,0.7))
				}

				#______
				#Add a legend
				if (leg[t]){
					if (mean(traits[as.factor(var.1) == levels(as.factor(var.1))[co],t], na.rm = T) <  mean(traits[,t], na.rm = T) ) {pos = "topright"}
					else{pos = "topleft"}
					legend(pos, inset = 0.05, levels(as.factor(var.2)), fill = col.dens, cex = cex.leg, bty = "n", ncol = round(sqrt(nlevels(as.factor((var.2))))-1))
				}

				#______
				#Plot the distribution by the factor 2
				for(s in 1:nlevels(as.factor(var.2))) {
					if (length(na.omit(traits[as.factor(var.1) == levels(as.factor(var.1))[co] & as.factor(var.2) == levels(as.factor(var.2))[s],t]))>1)
					{lines(density(traits[as.factor(var.1) == levels(as.factor(var.1))[co] & as.factor(var.2) == levels(as.factor(var.2))[s],t], na.rm = T), col = col.dens[s])}

					else if (length(na.omit(traits[as.factor(var.1) == levels(as.factor(var.1))[co] & as.factor(var.2) == levels(as.factor(var.2))[s],t])) == 1)
					{points(0,na.omit(traits[as.factor(var.1) == levels(as.factor(var.1))[co] & as.factor(var.2) == levels(as.factor(var.2))[s],t]), col = col.dens[s])}
				}
			}

		}
	}
	if(plot.ask | multipanel) {
		par(oldpar)
	}
}

# Ackerly & Cornwell 2007
# Plot populations values against species values
plotSpPop <- function(traits = NULL, ind.plot = NULL, sp = NULL, col.ind = rgb(0.5,0.5,0.5,0.5), col.pop = NULL, col.sp = NULL, col.site = NULL, resume = FALSE, p.val = 0.05, min.ind.signif = 10 , multipanel = TRUE, col.nonsignif.lm = rgb(0,0,0,0.5), col.signif.lm = rgb(1,0.1,0.1,0.8), silent = FALSE) {

	ntr <- dim(traits)[2]
	namestraits <- colnames(traits)


	traits <- traits[order(ind.plot),]

	ind.plot <- ind.plot[order(ind.plot)]
	sp <- sp[order(ind.plot)]

	name_sp_sites = paste(sp, ind.plot, sep = "_")
	comm <- t(table(ind.plot,1:length(ind.plot)))

	S <- colSums(comm>0)
	ncom <- length(S)


	plotsp = unlist(strsplit(levels(as.factor(name_sp_sites)),split = "_"))[seq(3,3*nlevels(as.factor(name_sp_sites)), by = 3)]
	#plosp is the plot in wich the population is
	plotind = unlist(strsplit(name_sp_sites,split = "_"))[seq(3,3*length(name_sp_sites), by = 3)]
	spplot = paste(unlist(strsplit(levels(as.factor(name_sp_sites)),split = "_"))[seq(1,3*nlevels(as.factor(name_sp_sites)), by = 3)], unlist(strsplit(levels(as.factor(name_sp_sites)),split = "_"))[seq(2,3*nlevels(as.factor(name_sp_sites)), by = 3)], sep = "_")

	traits_by_pop <- apply(traits,2,function(x) tapply(x, name_sp_sites,mean, na.rm = T))

	traits_by_sites <- apply(traits,2,function(x) tapply(x, ind.plot,mean, na.rm = T))


	if (multipanel){
		par(mfrow = c(sqrt(ntr), sqrt(ntr)) )
	}


	if (is.null(col.sp)){
		col.sp <- funky.col(nlevels(sp))
	}

	if (length(col.sp)<nlevels(sp)){
		col.sp <- rep(col.sp ,length.out = nlevels(sp))
	}

	if (is.null(col.pop)){
		col.pop <- col.sp[match(spplot, levels(sp))]
	}

	if (is.null(col.site)){
		col.site <- rep(rgb(0.1, 0.1, 0.1 , 0.8),ncom)
	}

	if (!resume){
		for(t in 1:ntr){
			x.ind <- traits_by_sites[match(plotind,rownames(traits_by_sites)),t]
			y.ind <- traits[,t]
			plot(x.ind, y.ind, pch = 16, col = col.ind, cex = 0.5)

			x.pop <- traits_by_sites[match(plotsp,rownames(traits_by_sites)),t]
			y.pop <- traits_by_pop[,t]
			points(x.pop, y.pop, pch = 16, col = col.pop)

			for(s in 1:nlevels(sp)){

				try( interm <- lm(y.pop[spplot == levels(sp)[s]] ~ x.pop[spplot == levels(sp)[s]]), silent = silent)

				options(warn = -1)
				if (class(try( interm <- lm(y.pop[spplot == levels(sp)[s]] ~ x.pop[spplot == levels(sp)[s]]), silent = silent)) == "lm"){
					if (!is.na(summary(interm)$coefficient[,4])){
						if (summary(interm)$coefficients[2,4]<p.val & length(interm$fitted.values)>min.ind.signif){
							lty.lm = 1
							lwd.lm = 3
							color.lm <- col.signif.lm
						}
						else{
							lty.lm = 0
							lwd.lm = 0
						}
						options(warn = 0)
					}


					else{
						lty.lm = 3
						lwd.lm = 1
						color.lm <- col.nonsignif.lm
					}

					if (!is.na(interm$coefficient[2])){
						lines(interm$model[,2] , interm$fitted.values, lty = lty.lm, lwd = lwd.lm, col = color.lm )
					}
				}
			}

			options(warn = -1)

			try( interm2 <- lm(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t] ~ traits_by_pop[,t]), silent = silent)

			color.lm2 <- col.nonsignif.lm

			if (class(try(lm(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t] ~ traits_by_pop[,t]), silent = silent)) == "lm"){
				if (!is.na(summary(interm2)$coefficient[,4])){
					if (summary(interm2)$coefficients[2,4]<p.val & length(interm2$fitted.values)>min.ind.signif){
						color.lm2 <- col.signif.lm
					}
					else{}
					options(warn = 0)
				}
			}


			points(tapply( traits_by_pop[,t], plotsp ,mean, na.rm = T) ~ tapply(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t], plotsp ,mean, na.rm = T) , col = col.site, pch = "*", cex = 3)

			points(tapply( traits_by_pop[,t], spplot ,mean, na.rm = T) ~ tapply(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t], spplot ,mean, na.rm = T) , col = col.sp, pch = 16, cex = 1.5)
			abline(try(lm( traits_by_sites[match(plotsp,rownames(traits_by_sites)),t] ~ traits_by_pop[,t]), silent = silent) , lty = 2, lwd = 2, col = color.lm2)

		}
	}

	if (resume){
		for(t in 1:ntr){
			x.pop <- traits_by_sites[match(plotsp,rownames(traits_by_sites)),t]
			y.pop <- traits_by_pop[,t]
			plot(x.pop, y.pop, pch = 16, col = col.pop)
			abline(a = 0, b = 1, lty = 3, lwd = 2)

			for(s in 1:nlevels(sp)){

				try( interm <- lm(y.pop[spplot == levels(sp)[s]] ~ x.pop[spplot == levels(sp)[s]]), silent = silent)

				options(warn = -1)
				if (class(try( interm <- lm(y.pop[spplot == levels(sp)[s]] ~ x.pop[spplot == levels(sp)[s]]), silent = silent)) == "lm"){
					if (!is.na(summary(interm)$coefficient[,4])){
						if (summary(interm)$coefficients[2,4]<p.val & length(interm$fitted.values)>min.ind.signif){
							lty.lm = 1
							lwd.lm = 3
						}
						else{
							lty.lm = 0
							lwd.lm = 0
						}
						options(warn = 0)
					}

					else{
						lty.lm = 3
						lwd.lm = 1
					}

					if (!is.na(interm$coefficient[2])){
						lines(interm$model[,2] , interm$fitted.values, lty = lty.lm, lwd = lwd.lm, col = col.sp[s] )
					}
				}
			}
			points(tapply( traits_by_pop[,t], spplot ,mean, na.rm = T) ~ tapply(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t], spplot ,mean, na.rm = T) , col = col.sp, pch = 16, cex = 1.5)
		}
	}


	par(mfrow = c(1, 1))

}



# Plot populations values against environmental variable(s)
plotSpVar <- function(traits = NULL, ind.plot = NULL, sp = NULL, variable = NULL, col.ind = rgb(0.5,0.5,0.5,0.5), col.pop = NULL, col.sp = NULL, col.site = NULL, resume = FALSE, p.val = 0.05, min.ind.signif = 10 , multipanel = TRUE, col.nonsignif.lm = rgb(0,0,0,0.5), col.signif.lm = rgb(1,0.1,0.1,0.8), silent = FALSE) {

	ntr <- dim(traits)[2]
	namestraits <- colnames(traits)

	traits <- traits[order(ind.plot),]

	ind.plot <- ind.plot[order(ind.plot)]
	sp <- sp[order(ind.plot)]

	name_sp_sites = paste(sp, ind.plot, sep = "_")
	comm <- t(table(ind.plot,1:length(ind.plot)))

	S <- colSums(comm>0)
	ncom <- length(S)


	plotsp = unlist(strsplit(levels(as.factor(name_sp_sites)),split = "_"))[seq(3,3*nlevels(as.factor(name_sp_sites)), by = 3)]
	#plosp is the plot in wich the population is
	plotind = unlist(strsplit(name_sp_sites,split = "_"))[seq(3,3*length(name_sp_sites), by = 3)]
	spplot = paste(unlist(strsplit(levels(as.factor(name_sp_sites)),split = "_"))[seq(1,3*nlevels(as.factor(name_sp_sites)), by = 3)], unlist(strsplit(levels(as.factor(name_sp_sites)),split = "_"))[seq(2,3*nlevels(as.factor(name_sp_sites)), by = 3)], sep = "_")

	traits_by_pop <- apply(traits,2,function(x) tapply(x, name_sp_sites,mean, na.rm = T))

	traits_by_sites <- variable

	if (is.vector(traits_by_sites)){
		traits_by_sites <- matrix(rep(traits_by_sites, times = ntr), ncol = ntr)
	}

	interm.for.names <- apply(traits,2,function(x) tapply(x, ind.plot,mean, na.rm = T))
	colnames(traits_by_sites) <- colnames(interm.for.names)

	rownames(traits_by_sites) <- rownames(interm.for.names)



	if (multipanel){
		par(mfrow = c(sqrt(ntr), sqrt(ntr)) )
	}


	if (is.null(col.sp)){
		col.sp <- funky.col(nlevels(sp))
	}

	if (length(col.sp)<nlevels(sp)){
		col.sp <- rep(col.sp ,length.out = nlevels(sp))
	}

	if (is.null(col.pop)){
		col.pop <- col.sp[match(spplot, levels(sp))]
	}

	if (is.null(col.site)){
		col.site <- rep(rgb(0.1, 0.1, 0.1 , 0.8),ncom)
	}

	if (!resume){
		for(t in 1:ntr){
			x.ind <- traits_by_sites[match(plotind,rownames(traits_by_sites)),t]
			y.ind <- traits[,t]
			plot(x.ind, y.ind, pch = 16, col = col.ind, cex = 0.5)

			x.pop <- traits_by_sites[match(plotsp,rownames(traits_by_sites)),t]
			y.pop <- traits_by_pop[,t]
			points(x.pop, y.pop, pch = 16, col = col.pop)

			for(s in 1:nlevels(sp)){

				try( interm <- lm(y.pop[spplot == levels(sp)[s]] ~ x.pop[spplot == levels(sp)[s]]), silent = silent)

				options(warn = -1)
				if (class(try( interm <- lm(y.pop[spplot == levels(sp)[s]] ~ x.pop[spplot == levels(sp)[s]]), silent = silent)) == "lm"){
					if (!is.na(summary(interm)$coefficient[,4])){
						if (summary(interm)$coefficients[2,4]<p.val & length(interm$fitted.values)>min.ind.signif){
							lty.lm = 1
							lwd.lm = 3
							color.lm <- col.signif.lm
						}
						else{
							lty.lm = 0
							lwd.lm = 0
						}
						options(warn = 0)
					}


					else{
						lty.lm = 3
						lwd.lm = 1
						color.lm <- col.nonsignif.lm
					}

					if (!is.na(interm$coefficient[2])){
						lines(interm$model[,2] , interm$fitted.values, lty = lty.lm, lwd = lwd.lm, col = color.lm )
					}
				}
			}

			options(warn = -1)

			try( interm2 <- lm(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t] ~ traits_by_pop[,t]), silent = silent)

			color.lm2 <- col.nonsignif.lm

			if (class(try(lm(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t] ~ traits_by_pop[,t]), silent = silent)) == "lm"){
				if (!is.na(summary(interm2)$coefficient[,4])){
					if (summary(interm2)$coefficients[2,4]<p.val & length(interm2$fitted.values)>min.ind.signif){
						color.lm2 <- col.signif.lm
					}
					else{}
					options(warn = 0)
				}
			}


			points(tapply( traits_by_pop[,t], plotsp ,mean, na.rm = T) ~ tapply(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t], plotsp ,mean, na.rm = T) , col = col.site, pch = "*", cex = 3)

			points(tapply( traits_by_pop[,t], spplot ,mean, na.rm = T) ~ tapply(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t], spplot ,mean, na.rm = T) , col = col.sp, pch = 16, cex = 1.5)
			abline(try(lm( traits_by_sites[match(plotsp,rownames(traits_by_sites)),t] ~ traits_by_pop[,t]), silent = silent) , lty = 2, lwd = 2, col = color.lm2)

		}
	}

	if (resume){
		for(t in 1:ntr){
			x.pop <- traits_by_sites[match(plotsp,rownames(traits_by_sites)),t]
			y.pop <- traits_by_pop[,t]
			plot(x.pop, y.pop, pch = 16, col = col.pop)
			abline(a = 0, b = 1, lty = 3, lwd = 2)

			for(s in 1:nlevels(sp)){

				try( interm <- lm(y.pop[spplot == levels(sp)[s]] ~ x.pop[spplot == levels(sp)[s]]), silent = silent)

				options(warn = -1)
				if (class(try( interm <- lm(y.pop[spplot == levels(sp)[s]] ~ x.pop[spplot == levels(sp)[s]]), silent = silent)) == "lm"){
					if (!is.na(summary(interm)$coefficient[,4])){
						if (summary(interm)$coefficients[2,4]<p.val & length(interm$fitted.values)>min.ind.signif){
							lty.lm = 1
							lwd.lm = 3
						}
						else{
							lty.lm = 0
							lwd.lm = 0
						}
						options(warn = 0)
					}

					else{
						lty.lm = 3
						lwd.lm = 1
					}

					if (!is.na(interm$coefficient[2])){
						lines(interm$model[,2] , interm$fitted.values, lty = lty.lm, lwd = lwd.lm, col = col.sp[s] )
					}
				}
			}
			points(tapply( traits_by_pop[,t], spplot ,mean, na.rm = T) ~ tapply(traits_by_sites[match(plotsp,rownames(traits_by_sites)),t], spplot ,mean, na.rm = T) , col = col.sp, pch = 16, cex = 1.5)
		}
	}
	par(mfrow = c(1, 1))
}



#______________#______________#______________#______________#______________#______________#______________#______________
#______________#______________#______________#______________#______________#______________#______________#______________
#__ Other functions

### Function to calculation SES on list of index
ses.listofindex <- function(index.list = NULL, val.quant = c(0.025, 0.975) ){

	namesindex.all <- names(index.list)
	nindex <- length(names(index.list))/2
	namesindex <- names(index.list)[seq(1,nindex*2, by = 2)]
	namestraits <- colnames(index.list[[1]])
	namescommunity <- rownames(index.list[[1]])

	ncom <- dim(index.list[[1]])[1]
	ntr <- dim(index.list[[1]])[2]

	#________________________________________
	#calculation of standardised effect size

	res <- list()
	for (i in seq(1,nindex*2, by = 2)){
		res[[eval(namesindex.all[i])]] <- ses(obs = index.list[[i]], nullmodel = index.list[[i+1]], val.quant = val.quant)
	}

	class(res) <- "ses.list"
	return(res)
}


ses.Tstats <- function(x, val.quant = c(0.025, 0.975)){

	if (!inherits(x, "Tstats")) {
		stop("x must be of class Tstats")
	}

	if(length(x$Tstats) != 6){
		stop("The object x of class Tstats does not contain null value from null models.")
	}

	res <- ses.listofindex(as.listofindex(x), val.quant = val.quant)
	return(res)
}

### Function to calculation SES
ses <- function(obs = NULL, nullmodel = NULL, val.quant = c(0.025, 0.975) ){

	if (is.vector(obs)){
		obs <- as.matrix(obs)
	}

	if (length(dim(obs)) != 2 ) {
		obs <- as.matrix(obs)
	}

	if (dim(obs)[1] == dim(obs)[2]) {
		warning("Observed matrix have the same number of rows and columns. The function is not able to detect automatically the correspondance between dimension of observed matrix and null model. You need to be sure that the null model is in the form of an array within the first and second dimension corresespond respectively to the first and second dimension of the observed matrix and the third dimension correspond to permutations")
		cond = c(1,2)

		if (class(nullmodel) == "list"){
			if (class(nullmodel[[1]]) == "list"){
				nullmodel <- array(unlist(nullmodel), dim = c(nrow(nullmodel[[1]][[1]]),ncol(nullmodel[[1]][[1]]), length(unlist(nullmodel))/nrow(nullmodel[[1]][[1]])/ncol(nullmodel[[1]][[1]])))
			}

			else {nullmodel <- array(unlist(nullmodel), dim = c(nrow(nullmodel[[1]]),ncol(nullmodel[[1]]), length(unlist(nullmodel))/nrow(nullmodel[[1]])/ncol(nullmodel[[1]])))}
		}
	}

	else{

		if (class(nullmodel) == "list"){
			if (class(nullmodel[[1]]) == "list"){
				nullmodel <- array(unlist(nullmodel), dim = c(nrow(nullmodel[[1]][[1]]),ncol(nullmodel[[1]][[1]]), length(unlist(nullmodel))/nrow(nullmodel[[1]][[1]])/ncol(nullmodel[[1]][[1]])))
			}

			else {nullmodel <- array(unlist(nullmodel), dim = c(nrow(nullmodel[[1]]),ncol(nullmodel[[1]]), length(unlist(nullmodel))/nrow(nullmodel[[1]])/ncol(nullmodel[[1]])))}
		}

		if (class(obs) == "list"){
			obs <- matrix(obs[[1]], nrow = nrow(obs[[1]]), ncol = ncol(obs[[1]]))
		}

		if (!is.null(dim(obs))) {
			cond <- c(NA,NA)

			if (dim(obs)[1] == dim(nullmodel)[1]){
			cond[1] <- 1
			}

			if (dim(obs)[1] == dim(nullmodel)[2]){
			cond[1] <- 2
			}

			if (length(dim(nullmodel)) == 3){
				if (dim(obs)[1] == dim(nullmodel)[3]){
				cond[1] <- 3
				}
			}

			if (dim(obs)[2] == dim(nullmodel)[1]){
			cond[2] <- 1
			}

			if (dim(obs)[2] == dim(nullmodel)[2]){
			cond[2] <- 2
			}

			if (length(dim(nullmodel)) == 3){
				if (dim(obs)[2] == dim(nullmodel)[3]){
				cond[2] <- 3
				}
			}
		}
	}

	cond <- na.omit(cond)

	res <- list()
	res$ses <- (obs-apply(nullmodel, cond, function(x) mean(x, na.rm = T)))/apply(nullmodel, cond, function(x) sd(x, na.rm = T))
	res$ses.inf <- (apply(nullmodel, cond, function(x) quantile(x, na.rm = T, prob = val.quant[1]))-apply(nullmodel, cond, function(x) mean(x, na.rm = T)))/apply(nullmodel, cond, function(x) sd(x, na.rm = T))
	res$ses.sup <- (apply(nullmodel, cond, function(x) quantile(x, na.rm = T, prob = val.quant[2]))-apply(nullmodel, cond, function(x) mean(x, na.rm = T)))/apply(nullmodel, cond, function(x) sd(x, na.rm = T))

	class(res) <- "ses"
	return(res)
}

RaoRel <- function(sample, dfunc, dphyl, weight = FALSE, Jost = FALSE, structure = NULL)  {
	####function Qdecomp by Villeger & Mouillot (J Ecol, 2008) modify by Wilfried Thuiller #####

	Qdecomp = function(functdist, abundances, w = TRUE) {

		# number and names of local communities
		c <- dim(abundances)[1]
		namescomm <- row.names(abundances)
		abundances <- as.matrix(abundances)

		# if necessary, transformation of functdist into matrix object
		if (is.matrix(functdist) == F) functdist <- as.matrix(functdist)

		# checking 'abundances' and 'functdist' dimensions
		if (dim(functdist)[1] != dim(functdist)[2]) stop("error : 'functdist' has different number of rows and columns")
		if (dim(abundances)[2] != dim(functdist)[1]) stop("error : different number of species in 'functdist' and 'abundances' ")

		# checking NA absence in 'functdist'
		if (length(which(is.na(functdist) == T)) != 0) stop("error : NA in 'functdist'")

		# replacement of NA by 0 in abundances
		if (is.na(sum(abundances)) == T) {
		for (i in 1:dim(abundances)[1])
		for (j in 1:dim(abundances)[2] )
		{ if (is.na(abundances[i,j]) == T) abundances[i,j] <- 0 } # end of i j
		} # end of if

		# species richness and total abundances in local communities
		abloc <- apply(abundances,1,sum)
		nbsploc <- apply(abundances,1,function(x) {length(which(x>0))} )

		# relative abundances inside each local community
		locabrel <- abundances/abloc

		# alpha diversity
		Qalpha = apply(locabrel, 1, function(x) t(x) %*% functdist %*% x)

		#Wc
		#Wc = abloc/sum(abloc)
		relsp <- apply(abundances,1,max)
		Wc = relsp/sum(relsp)

		# abundance-weighted mean alpha
		mQalpha <- as.numeric(Qalpha%*%abloc/sum(abloc) )
		#mQalpha <- as.numeric(Qalpha%*%relsp/sum(relsp) )

		#Villeger's correction
		if (w == T) {
			# abundance-weighted mean alpha
			mQalpha <- as.numeric(Qalpha%*%relsp/sum(relsp) )
			#mQalpha <- as.numeric(Qalpha%*%abloc/sum(abloc) )
			#totabrel <- apply(abundances,2,sum)/sum(abundances)
			totabrel <- apply(locabrel*Wc, 2, sum)
			Qalpha = Qalpha*Wc
		}

		# Rao's original definition: mean of Pi
		else {
			mQalpha <- mean(Qalpha)
			totabrel <- apply(locabrel,2,mean)
			}


		# gamma diversity
		Qgamma <- ( totabrel %*% functdist %*% totabrel ) [1]

		# beta diversity
		Qbeta <- as.numeric( Qgamma-mQalpha )

		# standardized beta diversity
		Qbetastd <- as.numeric(Qbeta/Qgamma )

		# list of results
		resQ <- list(Richness_per_plot = nbsploc, Relative_abundance = locabrel, Pi = totabrel, Wc = Wc, Species_abundance_per_plot = abloc, Alpha = Qalpha, Mean_alpha = mQalpha, Gamma = Qgamma, Beta = Qbeta, Standardize_Beta = Qbetastd )

		return(resQ)

	}


	###########function disc originally from S. Pavoine####

	disc = function (samples, dis = NULL, structures = NULL, Jost = F){
	  if (!inherits(samples, "data.frame"))
	    stop("Non convenient samples")
	  if (any(samples < 0))
	    stop("Negative value in samples")
	  if (any(apply(samples, 2, sum) < 1e-16))
	    stop("Empty samples")
	  if (!is.null(dis)) {
	    if (!inherits(dis, "dist"))
	      stop("Object of class 'dist' expected for distance")
	    # if (!is.euclid(dis))
	      #stop("Euclidean property is expected for distance")
	    dis <- as.matrix(dis)
	    if (nrow(samples) != nrow(dis))
	      stop("Non convenient samples")
	  }
		if (!is.null(structures)) {
	    if (!inherits(structures, "data.frame"))
	      stop("Non convenient structures")
	    m <- match(apply(structures, 2, function(x) length(x)),
	      ncol(samples), 0)
	    if (length(m[m == 1]) != ncol(structures))
	      stop("Non convenient structures")
	    m <- match(tapply(1:ncol(structures), as.factor(1:ncol(structures)),
	      function(x) is.factor(structures[, x])), TRUE, 0)
	    if (length(m[m == 1]) != ncol(structures))
	      stop("Non convenient structures")
	  }
	  Structutil <- function(dp2, Np, unit, Jost) {
	    if (!is.null(unit)) {
	      modunit <- model.matrix(~-1 + unit)
	      sumcol <- apply(Np, 2, sum)
	      Ng <- modunit * sumcol
	      lesnoms <- levels(unit)
	    }
	    else {
	      Ng <- as.matrix(Np)
	      lesnoms <- colnames(Np)
	    }
	    sumcol <- apply(Ng, 2, sum)
	    Lg <- t(t(Ng)/sumcol)
	    colnames(Lg) <- lesnoms
	    Pg <- as.matrix(apply(Ng, 2, sum)/nbhaplotypes)
	    rownames(Pg) <- lesnoms
	    deltag <- as.matrix(apply(Lg, 2, function(x) t(x) %*%
	      dp2 %*% x))
	    ug <- matrix(1, ncol(Lg), 1)
	    if (Jost) {
	      #dp2 <- as.matrix(as.dist(dfunct01))
	      deltag <- as.matrix(apply(Lg, 2, function(x) t(x) %*% dp2 %*% x))
	      X = t(Lg) %*% dp2 %*% Lg
	      alpha = 1/2 * (deltag %*% t(ug) + ug %*% t(deltag))
	      Gam = (X + alpha)/2
	      alpha = 1/(1-alpha) #Jost correction
	      Gam = 1/(1-Gam) #Jost correction
	      Beta_add = Gam - alpha
	      Beta_mult = 100*(Gam - alpha)/Gam
	    }
	    else {
	     deltag <- as.matrix(apply(Lg, 2, function(x) t(x) %*% dp2 %*% x))
	     X = t(Lg) %*% dp2 %*% Lg
	     alpha = 1/2 * (deltag %*% t(ug) + ug %*% t(deltag))
	     Gam = (X + alpha)/2
	     Beta_add = Gam - alpha
	     Beta_mult = 100*(Gam - alpha)/Gam
	    }
	    colnames(Beta_add) <- lesnoms
	    rownames(Beta_add) <- lesnoms
	    return(list(Beta_add = as.dist(Beta_add), Beta_mult = as.dist(Beta_mult),
	     Gamma = as.dist(Gam), Alpha = as.dist(alpha), Ng = Ng, Pg = Pg))
	  }
	  Diss <- function(dis, nbhaplotypes, samples, structures, Jost) {
	    structutil <- list(0)
	    structutil[[1]] <- Structutil(dp2 = dis, Np = samples, NULL, Jost)
	    diss <- list(structutil[[1]]$Alpha, structutil[[1]]$Gamma, structutil[[1]]$Beta_add, structutil[[1]]$Beta_mult)
	     if (!is.null(structures)) {
	      for (i in 1:length(structures)) {
	        structutil[[i + 1]] <- Structutil(as.matrix(structutil[[1]]$Beta_add),
	         structutil[[1]]$Ng, structures[, i], Jost)
	      }
	      diss <- c(diss, tapply(1:length(structures), factor(1:length(structures)),
	        function(x) as.dist(structutil[[x + 1]]$Beta_add)))
	    }
	    return(diss)
	  }
	  nbhaplotypes <- sum(samples)
	  diss <- Diss(dis, nbhaplotypes, samples, structures, Jost)
	  if (!is.null(structures)) {
	    names(diss) <- c("Alpha", "Gamma", "Beta_add", "Beta_prop", "Beta_region")
	    return(diss)
	  }
	  names(diss) <- c("Alpha", "Gamma", "Beta_add", "Beta_prop")
	  return(diss)
	}


	TD <- FD <- PD <- NULL

	#Taxonomic diversity
	dS <- matrix(1, nrow(sample), nrow(sample)) - diag(rep(1, nrow(sample)))
	temp_qdec <- Qdecomp(dS,t(sample), w = weight)  #Call the Qdecomp function for alpha, gamma and beta estimations.
	TD$Richness_per_plot = temp_qdec$Richness_per_plot
	TD$Relative_abundance = temp_qdec$Relative_abundance
	TD$Pi = temp_qdec$Pi
	TD$Wc = temp_qdec$Wc
	if (Jost){
		TD$Mean_Alpha = 1/(1-temp_qdec$Mean_alpha)
		TD$Alpha = 1/(1-temp_qdec$Alpha)
		TD$Gamma = 1/(1-temp_qdec$Gamma)
		TD$Beta_add = (TD$Gamma -TD$Mean_Alpha )
		TD$Beta_prop = 100*TD$Beta_add/TD$Gamma
		#Call the disc function for alpha, gamma and beta estimations for each pair of samples
		TD$Pairwise_samples <- disc(as.data.frame(sample), as.dist(dS), structure = structure, Jost = Jost)
		}
	else {
	  TD$Mean_Alpha = temp_qdec$Mean_alpha
		TD$Alpha = temp_qdec$Alpha
		TD$Gamma = temp_qdec$Gamma
		TD$Beta_add = (TD$Gamma -TD$Mean_Alpha )
		TD$Beta_prop = 100*TD$Beta_add/TD$Gamma
		#Call the disc function for alpha, gamma and beta estimations for each pair of samples
		TD$Pairwise_samples <- disc(as.data.frame(sample), as.dist(dS), structure = structure, Jost = Jost)
	}

	#Functional diversity estimation
	if (!is.null(dfunc)){
		FD <- list()
		if (Jost){
			if (max(dfunc)>1) dfunc <- dfunc/max(dfunc)  #Make sure the distance are between 0 and 1 for the Jost correction
			temp_qdec <- Qdecomp(dfunc,t(sample), w = weight)  #Call the Qdecomp function for alpha, gamma and beta estimations.
			# FD$Alpha = 1/(1-temp_qdec$Alpha)
			# FD$Mean_Alpha = mean(FD$Alpha)
			FD$Mean_Alpha = 1/(1-temp_qdec$Mean_alpha)
			FD$Alpha = 1/(1-temp_qdec$Alpha)
			FD$Gamma = 1/(1-temp_qdec$Gamma)
			#FD$Beta_add = (FD$Gamma -FD$Mean_Alpha )
			FD$Beta_prop = 100*FD$Beta_add/FD$Gamma
			#Call the disc function for alpha, gamma and beta estimations for each pair of samples
			FD$Pairwise_samples <- disc(as.data.frame(sample), as.dist(dfunc), structure = structure, Jost = Jost)
		}
		else {
			temp_qdec <- Qdecomp(dfunc,t(sample), w = weight) #Call the Qdecomp function for alpha, gamma and beta estimations.
			FD$Mean_Alpha = temp_qdec$Mean_alpha
			FD$Alpha = temp_qdec$Alpha
			FD$Gamma = temp_qdec$Gamma
			FD$Beta_add = (FD$Gamma -FD$Mean_Alpha )
			FD$Beta_prop = 100*FD$Beta_add/FD$Gamma
		  #FD$Beta = temp_qdec$Beta#
		  #Call the disc function for alpha, gamma and beta estimations for each pair of samples
		  FD$Pairwise_samples <- disc(as.data.frame(sample), as.dist(dfunc), structure = structure, Jost = Jost)
		}
	}

	#Phylogenetic diversity estimation
	if (!is.null(dphyl)){
	  PD <- list()
	  if (Jost){
			if (max(dphyl)>1) dphyl <- dphyl/max(dphyl)  #Make sure the distance are between 0 and 1 for the Jost correction
			temp_qdec <- Qdecomp(dphyl,t(sample), w = weight)  #Call the Qdecomp function for alpha, gamma and beta estimations.
			PD$Mean_Alpha = 1/(1-temp_qdec$Mean_alpha)
			PD$Alpha = 1/(1-temp_qdec$Alpha)
			PD$Gamma = 1/(1-temp_qdec$Gamma)
			PD$Beta_add = (PD$Gamma -PD$Mean_Alpha )
			PD$Beta_prop = 100*PD$Beta_add/PD$Gamma
			#Call the disc function for alpha, gamma and beta estimations for each pair of samples
			PD$Pairwise_samples <- disc(as.data.frame(sample), as.dist(dphyl), structure = structure, Jost = Jost)
		}
	  else {
			temp_qdec <- Qdecomp(dphyl,t(sample), w = weight) #Call the Qdecomp function for alpha, gamma and beta estimations.
			PD$Mean_Alpha = temp_qdec$Mean_alpha
			PD$Alpha = temp_qdec$Alpha
			PD$Gamma = temp_qdec$Gamma
			PD$Beta_add = (PD$Gamma -PD$Mean_Alpha )
			PD$Beta_prop = 100*PD$Beta_add/PD$Gamma
			#PD$Beta = temp_qdec$Beta
			#Call the disc function for alpha, gamma and beta estimations for each pair of samples
			PD$Pairwise_samples <- disc(as.data.frame(sample), as.dist(dphyl), structure = structure, Jost = Jost)
	  }
	}
	out <- list(TD, FD, PD)
	names(out) <- c("TD", "FD", "PD")
	return(out)
}

###################################################################################################################################
# 	The Rao function computes alpha, gamma and beta-components for taxonomic, functional and phylogenetic diversity with the Rao index
# 	The script integrates two functions: "Qdecomp", by Villeger & Mouillot (J Ecol, 2008) modify by Wilfried Thuiller, and "disc", by S. Pavoine, in the package ade4.
# 	For a regional assemblage of C local communities gamma = mean(alpha) + beta, where:
# 	gamma is the diversity of the regional pool
# 	alpha are the diversities of the local communities
# 	beta is the turn over between local communities
# 	diversity is estimated with the Rao quadratic entropy index (Rao 1982)
#
# INPUTS:
#	- "abundances": matrix of abundances (c x s) of the s species for the c local communities (or samples)
#	- "dfunct": matrix (s x s) or dist object with pairwise functional trait distances between the s species
#	- "dphyl": as dfunct but for phylogenetic distances
#	- "weight": defining if the correction by Villeger & Mouillot (J Ecol, 2008) is applied or not
#	- "Jost": defining if the Jost correction is applied (this paper and Jost 2007)
#	- "structure": a data frame containing the name of the group to which samples belong see
#   NA are not allowed in 'locabrel <- abundances/ablocist'
#   NA are automatically replaced by 0 in 'abundances'
#
# OUTPUTS:
#	- The results are organized for Taxonomic diversity ($TD), Functional diversity ($FD) and phylogenetical diversity ($PD). Beta and gamma diversities are calculated for the whole data set and for each pair of samples ("Pairwise_samples")
#	- "$Richness_per_plot"(number of species per sample)
#	- "$Relative_abundance" (species relative abundances per plot)
#	- "$Pi" (species regional relative abundance)
#	- "$Wc" (weigthing factor),
#	- "$Mean_Alpha" (mean aplpha diversity; for taxonomic diversity the Simpson index is calculated)
#	- "$Alpha" (alpha diversity for each sample; for taxonomic diversity the Simpson index is calculated)
#	- "$Gamma" (gamma diversity; for taxonomic diversity the Simpson index is calculated)
#	- "$Beta_add" (Gamma-Mean_Alpha)
#	- "$Beta_prop" (Beta_add*100/Gamma)
#	- "$Pairwise_samples$Alpha" (mean alpha for each pair of samples)
#	- "$Pairwise_samples$Gamma" (gamma for each pair of samples)
#	- "$Pairwise_samples$Beta_add" (beta for each pair of samples as Gamma-Mean_Alpha)
#	- "$Pairwise_samples$Beta_prop" (beta for each pair of samples as Beta_add*100/Gamma)

#####################################################################################################################################





#Function to plot result of observed indices values against null distribution
#Use the function plot(as.randtest(x)) from ade4
plotRandtest <- function(x, alternative = "two-sided", ...){

	alternative <- match.arg(alternative, c("greater", "less", "two-sided"), several.ok = TRUE)

	if (!inherits(x, "listofindex")) {
		if (inherits(x, "Tstats") | inherits(x, "ComIndex") | inherits(x, "ComIndexMulti")) {
			x <- as.listofindex(x)
		}
		else{stop("x must be a list of objects of class listofindex, Tstats, ComIndex or ComIndexMulti")}
	}

	index.list <- x

	oldpar <- par(no.readonly = TRUE)

	namesindex.all <- names(index.list)
	nindex <- length(names(index.list))/2
	namesindex <- names(index.list)[seq(1,nindex*2, by = 2)]
	namestraits <- colnames(index.list[[1]])
	namescommunity <- rownames(index.list[[1]])

	ncom <- c()
	ntr <- c()
	for(i in seq(1, 2*nindex, by = 2)){
		ncom <- c(ncom,dim(as.matrix(index.list[[i]]))[1])
		ntr <- c(ntr,dim(as.matrix(index.list[[i]]))[2])
	}

	if (is.null(ncom)) {ncom <- dim(as.matrix(index.list[[1]]))[1]}
	if (is.null(ntr)) {ntr <- dim(as.matrix(index.list[[1]]))[2]}

	if (is.null(ncom)) {ncom = 1}
	if (is.null(ntr)) {ntr = 1}

	for (i in seq(1, nindex*2, by = 2)){
		for (t in 1:ntr[1]){

			for(s in 1: ncom[1]){

				rt <- as.randtest(sim = na.omit(index.list[[i+1]][,t,s]), obs = na.omit(index.list[[i]][s,t]), alter = alternative)
				plot(rt, main = paste(namesindex.all[i], namestraits[t], namescommunity[s], "p.value = ", round(rt$pvalue, digits = 5)), ...)
			}

			rt <- as.randtest(sim = index.list[[i+1]][t,,], obs = index.list[[i]][t], alter = alternative)

			plot(rt, main = paste(namesindex.all[i], namestraits[t], "p.value = ", round(rt$pvalue, digits = 5)), ...)
		}
	}

}



#Replace a matrix with abundance of species and mean traits by pop in a pseudo-individual matrix
#Each individual take therefore the value of the population

AbToInd <- function(traits, com, type.sp.val = "count"){

	type.sp.val <- match.arg(type.sp.val, c("count", "abundance"))
	if (nrow(traits) != nrow(com)){
		stop("number of rows of traits and com need to be equal")
	}

	#type.sp.val is either count data or abundance
	#transform abundance data to number of individual
	#Using abundance type.sp.val is EXPERIMENTAL. This function round abundance to fit count data.

	if (type.sp.val == "abundance"){
		if (min(com)<1){
			com <- apply(com,2, function(x) round(x/min(x[x>0], na.rm = T), 1)*10 )
		}
	}

	ntr <- ncol(traits)

	#number of individual to individual traits data
	x2 <- c()
	x4 <- c()
	x6 <- c()
	for(i in 1 : nrow(com)){

		x1 <- matrix(rep(traits[i, ], rowSums(com)[i]), nrow = ntr, ncol = rowSums(com)[i])
		x2 <- cbind(x2,x1)
		x3 <- rep(rownames(com)[i], rowSums(com)[i])
		x4 <- c(x4,x3)

		for(co in 1:ncol(com)){
			x5 <- rep(colnames(com)[co], com[i,co])
			x6 <- c(x5,x6)
		}
	}

	res <- list()

	res$traits <- data.frame(t(x2))
	res$sp <- as.factor(x4)
	res$ind.plot <- as.factor(x6)

	return(res)

}


########################################################################
###### 				    Metrics								 ######
########################################################################

#calculation of the coefficient of variation of the NND (nearest neigbour distance) with our without division by the range of the trait
CVNND <- function(traits, div_range = FALSE, na.rm = FALSE, scale.tr = TRUE, method.dist = "euclidian"){

	#Uni-traits vector transforme to a matrix
	if(!is.matrix(traits)){
		if(na.rm) {
			traits <- na.omit(traits)
		}
		traits <- as.matrix(traits, ncol = 1)
	}

	#Calcul of nearest neighbourhood distance
	if(scale.tr) {
		traits <- apply(traits, 2, scale)
	}

	mat.dist <- as.matrix(dist(traits, method = method.dist))
	diag(mat.dist )<- NA
	nnd <- apply(mat.dist, 1, function(x) min (x, na.rm = na.rm))

	#Metric calculation
	CVNND <- sd(nnd[is.finite(nnd)], na.rm = na.rm)/mean(nnd[is.finite(nnd)], na.rm = na.rm)

	if (div_range) {
		CVNND <- CVNND/(max(traits, na.rm = na.rm) - min(traits, na.rm = na.rm) )
	}
	else {}

	return(CVNND)
}

#calculation of mean NND (nearest neigbour distance) for one trait or multiple traits with our without division by the range of the trait
MNND <- function(traits, div_range = FALSE, na.rm = FALSE, scale.tr = TRUE, method.dist = "euclidian"){

	#Uni-traits vector transforme to a matrix
	if(is.vector(traits)){
		traits <- as.matrix(traits, ncol = 1)
	}

	#Calcul of nearest neighbourhood distance
	if(scale.tr) {
		traits <- apply(traits, 2, scale)
	}

	mat.dist <- dist(traits, method = method.dist)
	nnd <- apply(as.matrix(mat.dist), 1, function(x) min (x[x>0], na.rm = T))

	#Metric calculation
	MNND <- mean(nnd[is.finite(nnd)], na.rm = na.rm)

	if (div_range) {
		MNND <- MNND/(max(traits, na.rm = na.rm) - min(traits, na.rm = na.rm) )
	}
	else {}

	return(MNND)
}

#calculation of sd NND for one trait or multiple traits with our without division by the range of the trait
SDNND <- function(traits, div_range = FALSE, na.rm = FALSE, scale.tr = TRUE, method.dist = "euclidian"){

	#Uni-traits vector transforme to a matrix
	if(is.vector(traits)){
		traits <- as.matrix(traits, ncol = 1)
	}

	#Calcul of nearest neighbourhood distance
	if(scale.tr) {
		traits <- apply(traits, 2, scale)
	}

	mat.dist <- dist(traits, method = method.dist)
	nnd <- apply(as.matrix(mat.dist), 1, function(x) min (x[x>0], na.rm = T))

	#Metric calculation
	SDNND <- sd(nnd[is.finite(nnd)], na.rm = na.rm)

	if (div_range) {
		SDNND <- SDNND/(max(traits, na.rm = na.rm) - min(traits, na.rm = na.rm) )
	}
	else {}

	return(SDNND)
}

#calculation of minimum NND for one trait or multiple traits with our without division by the range of the trait
MinNND <- function(traits, div_range = FALSE, na.rm = FALSE, scale.tr = TRUE, method.dist = "euclidian"){

	#Uni-traits vector transforme to a matrix
	if(is.vector(traits)){
		traits <- as.matrix(traits, ncol = 1)
	}

	#Calcul of nearest neighbourhood distance
	if(scale.tr) {
		traits <- apply(traits, 2, scale)
	}

	mat.dist <- dist(traits, method = method.dist)
	nnd <- apply(as.matrix(mat.dist), 1, function(x) min (x[x>0], na.rm = T))

	#Metric calculation
	MinNND <- min(nnd[is.finite(nnd)], na.rm = na.rm)

	if (div_range) {
		MinNND <- MinNND/(max(traits, na.rm = na.rm) - min(traits, na.rm = na.rm) )
	}
	else {}

	return(MinNND)
}


#calculation of sd ND (neigbour distance) for one trait with our without division by the range of the trait
SDND <- function(trait, div_range = FALSE, na.rm = FALSE){

	r = sort(trait)

	SDND <- sd(diff(r, na.rm = na.rm), na.rm = na.rm)

	if (div_range) {
		SDND <- SDND/(max(trait, na.rm = na.rm) - min(trait, na.rm = na.rm) )
	}
	else {}

	return(SDND)
}

#calculation of sd ND (neigbour distance) for one trait with our without division by the range of the trait
MND <- function(trait, div_range = FALSE, na.rm = FALSE){

	r = sort(trait)

	MND <- mean(diff(r, na.rm = na.rm), na.rm = na.rm)

	if (div_range) {
		MND <- MND/(max(trait, na.rm = na.rm) - min(trait, na.rm = na.rm) )
	}
	else {}

	return(MND)
}

#Sum of branch length
#All individuals with one NA are omitted

SumBL <- function(traits, gower.dist = TRUE, method.hclust = "average", scale.tr = TRUE, method.dist = "euclidian"){
	#require(vegan)

	if(sum(is.na(traits)) > 0){
		traits <- na.omit(traits)
		warning("All individuals with one NA are excluded")
	}

	if(gower.dist){
		#require(FD)
		tree.dist <- gowdis(traits)
	}

	else{
		if(scale.tr){
			traits <- apply(traits, 2, scale)
		}
		tree.dist <- dist(traits, method = method.dist)
	}

	tree.dist <- na.omit(tree.dist)
	tree <- hclust(tree.dist, method = method.hclust)

	res <- treeheight(tree)
	return(res)
}


#Min/Max MST
MinMaxMST <- function (traits, gower.dist = TRUE, scale.tr = TRUE, method.dist = "euclidian"){

	if (sum(is.na(traits))>0) {
		warning(paste("This function exclude", sum(is.na(traits)),"Na value", sep=" "))
	}

	if(gower.dist){
		#require(FD)
		tree.dist <- gowdis(traits)
	}
	else{
		if(scale.tr){
			traits <- apply(traits, 2, scale)
		}
		tree.dist <- dist(traits, method = method.dist)
	}

	traits.mst <- mst(as.matrix(tree.dist))
	res.mst <- as.matrix(tree.dist)*as.matrix(traits.mst)

	res.mst <- apply(res.mst, 1, function(x) as.numeric(sub("^0$",NA,x)))
	res <- min(res.mst, na.rm = T) / max(res.mst, na.rm = T)
	return(res)
}


IndexByGroups <- function(metrics, groups){
	res <- c()
	for (i in 1:length(metrics)){
		res <- c(res, paste("tapply(x, ", groups, ", function(x) ", metrics[i], ")", sep = ""))
	}

	return(res)
}



##########################################################################################################################
# function to compute the 3 functional diversity indices presented in Villeger et al. 2008 (Ecology 89 2290-2301):    #
#    Functional richness (FRic), Functional evenness (FEve), Functional divergence (FDiv)              #
#                                                            #
# This function requires R libraries 'geometry' and 'ape'                                #
#                                                            #
#                                                            #
# inputs:                                                        #
#    - "traits" = functional matrix (S x T) with values of the T functional traits for the S species of interest   #
# 	    NA are not allowed                                              #
# 	    for each trait, values are standardized (mean = 0 and standard deviation = 1)                #
#      for FRic computation, number of species must be higher than number of traits (S>T)              #
#                                                            #
#    - "abundances" = abundance matrix (C x S) with the abundances of the S species in the C commnuities       #
# 	    NA are automatically replaced by 0                                      #
#                                                            #
#                                                            #
# outputs:                                                        #
#     list of 4 vectors with values of indices in the C communities (names are those given in 'abundances')     #
#        - nbsp: number of species                                        #
#        - FRic: functional richness index                                    #
#        - FEve: functional evenness index                                    #
#        - FDiv: functional divergence index                                   #
#                                                            #
##########################################################################################################################

Fred <- function(traits, ind.plot = NULL) {
	
	if(is.null(ind.plot)) {
		ind.plot <- rep("all plots", dim(traits)[1])
		warnings("The argument `ind.plot` is null. Only one value will be compute for each metrics.") 
	}
	
	abundances <- table(ind.plot, 1:length(ind.plot))

	#loading required libraries
	#require(geometry)

	# T = number of traits
	T <- dim(traits)[2]

	# c = number of communities
	C <- dim(abundances)[1]

	# check coherence of number of species in 'traits' and 'abundances'
	if (dim(abundances)[2] != dim(traits)[1]) stop(" error : different number of individuals in 'traits' and in 'ind.plot' ")

	# check absence of NA in 'traits'
	if (length(which(is.na(traits) == T)) != 0) stop(" error : NA in 'traits' matrix ")

	# replacement of NA in 'abundances' by '0'
	abundances[which(is.na(abundances))] <- 0

	# definition of vector for results, with communities'names as given in 'abundances'
	nbsp <- rep(NA,C) ; names(nbsp) <- row.names(abundances)
	FRic <- rep(NA,C) ; names(FRic) <- row.names(abundances)
	FEve <- rep(NA,C) ; names(FEve) <- row.names(abundances)
	FDiv <- rep(NA,C) ; names(FDiv) <- row.names(abundances)

	# scaling and centering of each trait according to all species values
	traitsCS <- scale(traits, center = TRUE, scale = TRUE)


	############################################################################################################
	# loop to compute on each community the three Functional Diversity Indices, plus the number of species

	if(C>1) {pb <- txtProgressBar(1, C, style = 3)}

	for (i in 1:C)
	{
	 # selection of individuals present in the community
	 esppres <- which(abundances[i,]>0)

	 # number of individuals in the community
	 S <- length(esppres)
	 nbsp[i] <- S

	 # check if
	 if (S <= T) {stop(" number of individuals must be higher than number of traits ")}

	 # filter on 'traits' and 'abundances' to keep only values of individuals present in the community
	 tr <- traitsCS[esppres,]
	 ab <- as.matrix(abundances[i,esppres])

	 # scaling of abundances
	 abondrel <- ab/sum(ab)

	 # functional diversity indices
	   # FRIC
	    # use of convhulln function
	      # volume
	      FRic[i] <- round(convhulln(tr,"FA")$vol,6)

	      # identity of vertices
	      vert0 <- convhulln(tr,"Fx TO 'vert.txt'")
	      vert1 <- scan("vert.txt", quiet = T)
	 	    vert2 <- vert1+1
	      vertices <- vert2[-1]

	   # FEve
	    # computation of inter-species euclidian distances
	       distT <- dist(tr, method = "euclidian")
	    # computation of Minimum Spanning Tree and conversion of the 'mst' matrix into 'dist' class
	       linkmst <- mst(distT) ; mstvect <- as.dist(linkmst)
	    # computation of the pairwise cumulative relative abundances and conversion into 'dist' class
	       abond2 <- matrix(0,nrow = S,ncol = S)
	       for (q in 1:S)
	       for (r in 1:S)
	       abond2[q,r] <- abondrel[q]+abondrel[r]
	       abond2vect <- as.dist(abond2)
	    # computation of EW for the (S-1) segments pour relier les S points
	       EW <- rep(0,S-1)

	       flag <- 1
	       for (m in 1:((S-1)*S/2))
	       {if (mstvect[m] != 0) {EW[flag] <- distT[m]/(abond2vect[m]) ; flag <- flag+1}}
	    # computation of the PEW and comparison with 1/S-1, finally computation of FEve
	       minPEW <- rep(0,S-1) ; OdSmO <- 1/(S-1)
	       for (l in 1:(S-1))
	       minPEW[l] <- min( (EW[l]/sum(EW)) , OdSmO )
	       FEve[i] <- round( ( (sum(minPEW))- OdSmO) / (1-OdSmO ) ,6)


	   # FDiv

	    # traits values of vertices
	       trvertices <- tr[vertices,]

	    # coordinates of the center of gravity of the vertices (Gv)
	       baryv <- apply(trvertices,2,mean)

	    # euclidian dstances to Gv (dB)
	       distbaryv <- rep(0,S)
	       for (j in 1:S)
	       distbaryv[j] <- ( sum((tr[j,]-baryv)^2) )^0.5

	    # mean of dB values
	       meandB <- mean(distbaryv)

	    # deviations to mean of db
	       devdB <- distbaryv-meandB

	    # relative abundances-weighted mean deviation
	       abdev <- abondrel*devdB

	    # relative abundances-weighted mean of absolute deviations
	       ababsdev <- abondrel*abs(devdB)

	    # computation of FDiv
	       FDiv[i] <- round( (sum(abdev)+meandB) / (sum(ababsdev)+meandB) ,6)

		#Show the progress
		Sys.sleep(0.02)
		if(C>1) {setTxtProgressBar(pb, i)}

	} # end of i

	if(C>1) {close(pb)}
	unlink("vert.txt")
	res <- list(nbind = nbsp, FRic = FRic, FEve = FEve, FDiv = FDiv )
	invisible(res)

}# end of function

#Calcul pvalue for a list of index. This test equates to finding the quantile in exp in which obs would be found (under a one-tailed test)
Pval <- function (x, na.rm = TRUE) {

	if (!inherits(x, "listofindex")){
		x <- as.listofindex(x)
	}

	res <- list()

	for (i in seq(1, length(x), 2)){
		obs <- x[[i]]
		nullmodel <- x[[i+1]]

		if (is.vector(obs)){
			obs <- t(as.matrix(obs))
		}

		if (length(dim(obs)) != 2 ) {
			obs <- t(as.matrix(obs))
		}

		if (dim(obs)[1] == dim(obs)[2]) {
			warning("Observed matrix have the same number of rows and columns. The function is not able to detect automatically the correspondance between dimension of observed matrix and null model. You need to be sure that the null model is in the form of an array which the first and second dimension correspond respectively to the first and second dimension of the observed matrix and the third dimension correspond to the permutations")
			cond = c(1,2)
		}

		if (dim(obs)[1] != dim(obs)[2]) {

			if (class(nullmodel) == "list"){
				if (class(nullmodel[[1]]) == "list"){
					nullmodel <- array(unlist(nullmodel), dim = c(nrow(nullmodel[[1]][[1]]),ncol(nullmodel[[1]][[1]]), length(unlist(nullmodel))/nrow(nullmodel[[1]][[1]])/ncol(nullmodel[[1]][[1]])))
				}

				else {nullmodel <- array(unlist(nullmodel), dim = c(nrow(nullmodel[[1]]),ncol(nullmodel[[1]]), length(unlist(nullmodel))/nrow(nullmodel[[1]])/ncol(nullmodel[[1]])))}
			}

			if (class(obs) == "list"){
				obs <- matrix(obs[[1]], nrow = nrow(obs[[1]]), ncol = ncol(obs[[1]]))
			}

			if (!is.null(dim(obs))) {
				cond <- c(NA,NA)

				if (dim(obs)[1] == dim(nullmodel)[1]){
				cond[1] <- 1
				}

				if (dim(obs)[1] == dim(nullmodel)[2]){
				cond[1] <- 2
				}

				if (length(dim(nullmodel)) == 3){
					if (dim(obs)[1] == dim(nullmodel)[3]){
					cond[1] <- 3
					}
				}

				if (dim(obs)[2] == dim(nullmodel)[1]){
				cond[2] <- 1
				}

				if (dim(obs)[2] == dim(nullmodel)[2]){
				cond[2] <- 2
				}

				if (length(dim(nullmodel)) == 3){
					if (dim(obs)[2] == dim(nullmodel)[3]){
					cond[2] <- 3
					}
				}
			}
		}

		cond <- na.omit(cond)

		nullmodel <- aperm(nullmodel, c(cond, (1:3)[!1:3 %in% cond]))

		ni_tot <- names(x[[1]])[1]
		if(is.null(ni_tot)){
				ni_tot <- dimnames(x[[1]])[[2]][1]
		}

		res[[eval(names(x[i]))]] <- list()
		for(t in 1:length(eval(ni_tot))){

			ni <- names(x[[i]])[t]
			if(is.null(ni)){
				ni <- dimnames(x[[i]])[[2]][t]
			}

			res[[eval(names(x[i]))]][[eval(ni)]] <- list()

			for(g in 1:length(obs[,t])){
				res[[eval(names(x[i]))]][[eval(ni)]]$pval.inf[g] <- (1 + sum(obs[g,t] < nullmodel [g,t,], na.rm = T)) / (dim(nullmodel)[3] + 1)
				res[[eval(names(x[i]))]][[eval(ni)]]$pval.sup[g] <- (1 + sum(obs[g,t] > nullmodel [g,t,], na.rm = T)) / (dim(nullmodel)[3] + 1)
			}

			if(!is.null(dimnames(x[[1]])[[1]])){
				res[[eval(names(x[i]))]][[eval(ni)]]$pval.inf <- cbind(res[[eval(names(x[i]))]][[eval(ni)]]$pval.inf)
				rownames(res[[eval(names(x[i]))]][[eval(ni)]]$pval.inf) <- dimnames(x[[1]])[[1]]

				res[[eval(names(x[i]))]][[eval(ni)]]$pval.sup <- cbind(res[[eval(names(x[i]))]][[eval(ni)]]$pval.sup)
				rownames(res[[eval(names(x[i]))]][[eval(ni)]]$pval.sup) <- dimnames(x[[1]])[[1]]
			}

			if(!is.null(dimnames(x[[1]])[[2]][t])){
				colnames(res[[eval(names(x[i]))]][[eval(ni)]]$pval.sup) <- dimnames(x[[1]])[[2]][t]
				colnames(res[[eval(names(x[i]))]][[eval(ni)]]$pval.inf) <- dimnames(x[[1]])[[2]][t]
			}
		}
  }
  return(res)
}

#Function for colors from adegenet
funky.col <- colorRampPalette(c("#A6CEE3","#1F78B4","#B2DF8A",
               "#33A02C","#FB9A99","#E31A1C",
               "#FDBF6F","#FF7F00","#CAB2D6",
               "#6A3D9A","#FFFF99","#B15928"))


RandCom <- function(Ncom = 10, Nsp = 20, Nind.com = 100, sdlog = 1.5, min_value_traits = 80, max_value_traits = 200, cv_intra_sp = 1.5, cv_intra_com = 1.5, Int_Filter_Strength = 50, Ext_Filter_Strength = 50, Filter = "None"){

	Filter <- match.arg(Filter, c("None", "External", "Internal", "Both"))

	set <- letters
	SET <- LETTERS
	sp_Letters <- combn(set, round(Nsp/26/26)+1)
	com_Letters <- combn(SET, round(Ncom/26/26)+1)

	com <- paste("com", com_Letters[seq(1:Ncom)], sep ="_")
	sp <- paste("sp",  sp_Letters[seq(1:Nsp)], sep ="_")

	Nind <- Ncom * Nind.com
	trait_distri <- c("rnorm", "rlnorm") #For latter: may be an argument of the function
	Ntr <- length(trait_distri)

	Data <- data.frame(matrix(0, ncol=2+Ntr, nrow=Nind)); colnames(Data) = c("com","sp","trait1","trait2")

 #____
 # 1 - Assign a community to each individual
	### HYP 1: all communities have the same number of individuals
	Data$com = as.factor(rep(com, rep(Nind.com,length(com))))

 #____
 # 2 - Assign a species to each individual
	### HYP 1: species abudances follow a lognormal distribution within each comm
	### HYP 2: all species can be as abundant as each other within the regional pool (no rare/common species)
	### HYP 3: species abundance is not linked to any type of gradient -> could we influence the distribution by both lognormal and gradient?
	for(c in 1:Ncom){
		ex.sp = sample(sp, size = Nind.com, prob = rlnorm(Nsp, 0, sdlog), replace = T)
		Data$sp[((c-1)*Nind.com+1):(c*Nind.com)] = ex.sp
	}

 #____
 # 3 - Assign a trait value to each individual
  ## OPTION 1: random ##
	if(Filter=='None'){
		#trait1: normal distribution
		Data$trait1 = rnorm(Nind, (max_value_traits-min_value_traits), (max_value_traits-min_value_traits) * cv_intra_sp)
		#trait2: uniform distribution
		Data$trait2 <- runif(Nind, min_value_traits, max_value_traits)
	}

  ## OPTION 2: internal filtering ##
	if(Filter=='Internal'){
		# Parameter for the distance between species mean trait values:
		### HYP 1: the most extreme case (if Int_Filter_Strength==100) species have equally distributed mean values along the trait gradient
		Init_sp_mean <- round(sort(seq(min_value_traits, max_value_traits, length.out = Nsp)), 2)

		# Defining traits mean by species
		mean.sp <- sample(rnorm(Nsp, mean = Init_sp_mean, sd = (100-Int_Filter_Strength) * Nsp / (max_value_traits-min_value_traits)), replace=F)
		#mean.sp[mean.sp<0] = runif(sum(mean.sp<0), min_value_traits, max_value_traits) ## to avoid negative values!!!
			#--- I think it's not necessary if we use the absolute value for the variance

		# Draw the individual traits values depending on species attributes
		for(i in 1:Nind){
			#trait 1 : normal distribution
			Data$trait1[i] = rnorm(1, mean.sp[which(Data$sp[i] == sp)], abs(mean.sp[which(Data$sp[i] == sp)])*cv_intra_sp)
			#trait 2 : uniform distribution
			Data$trait2[i] = runif(1, mean.sp[which(Data$sp[i] == sp)]*(1-cv_intra_sp), abs(mean.sp[which(Data$sp[i] == sp)])*(1+cv_intra_sp))
		}
	}

  ## OPTION 3: external filtering ##
	if(Filter=='External'){

			  # Parameter for the distance between communities mean trait values:
			  ### HYP 1: the most extreme case (if Ext_Filter_Strength==100) communities have equally distributed mean values along the trait gradient
			  Init_com_mean <- round(sort(seq(min_value_traits, max_value_traits, length.out = Ncom)), 2)

			  # Defining traits mean by community
			  mean.com <- sample(rnorm(Ncom, mean = Init_com_mean, sd = (100-Ext_Filter_Strength) * Nsp / (max_value_traits-min_value_traits)), replace=F)
			  #mean.com[mean.com<0] = runif(sum(mean.com<0),min_value_traits, max_value_traits) ## to avoid negative values!!!

			  # Draw the individual traits depending on communities attributes
			  for(i in 1:Nind){
					#trait 1 : normal distribution
					Data$trait1[i] = rnorm(1, mean.com[which(Data$com[i] == com)], abs(mean.com[which(Data$com[i] == com)])*cv_intra_com)
					#trait 2 : uniform distribution
					Data$trait2[i] = runif(1, mean.com[which(Data$com[i] == com)]*(1-cv_intra_com), abs(mean.com[which(Data$com[i] == com)])*(1+cv_intra_com))
			  }
		}

  ## OPTION 4: external and internal filtering ##
	if(Filter=='Both'){

		# Parameter for the distance between species mean trait values:
		### HYP 1: the most extreme case (if Int_Filter_Strength==100) species have equally distributed mean values along the trait gradient
		Init_sp_mean <- round(sort(seq(min_value_traits, max_value_traits, length.out = Nsp)), 2)

		# Defining traits mean by species
		mean.sp <- sample(rnorm(Nsp, mean = Init_sp_mean, sd = (100-Int_Filter_Strength) * Nsp / (max_value_traits-min_value_traits)), replace = FALSE)

		# Parameter for the distance between communities mean trait values:
		### HYP 1: the most extreme case (if Ext_Filter_Strength==100) communities have equally distributed mean values along the trait gradient
		Init_com_mean <- round(sort(seq(min_value_traits, max_value_traits, length.out = Ncom)), 2)

		# Defining traits mean by community
		mean.com <- sample(rnorm(Ncom, mean = Init_com_mean, sd = (100-Ext_Filter_Strength) * Nsp / (max_value_traits-min_value_traits)), replace = FALSE)

		# Draw the individual traits values depending on species attributes
		for(i in 1:Nind){
			#trait 1 : normal distribution
			Data$trait1[i] = 1/2*(rnorm(1, mean.sp[which(Data$sp[i] == sp)], abs(mean.sp[which(Data$sp[i] == sp)])*cv_intra_sp) +
							 rnorm(1, mean.com[which(Data$com[i] == com)], abs(mean.com[which(Data$com[i] == com)])*cv_intra_com))# - mean(mean.com)
			#trait 2 : uniform distribution
			Data$trait2[i] = 1/2*(runif(1, mean.sp[which(Data$sp[i] == sp)]*(1-cv_intra_sp), abs(mean.sp[which(Data$sp[i] == sp)])*(1+cv_intra_sp)) +
							 runif(1, mean.com[which(Data$com[i] == com)]*(1-cv_intra_com), abs(mean.com[which(Data$com[i] == com)])*(1+cv_intra_com)))# - mean(mean.com)
		}
	}
	res <- list()
	res$data <- Data
	res$call <- match.call()
	return(res)
}


samplingSubsetData <- function(d = NULL, sampUnit = NULL, nperm = 9, type = "proportion", prop = seq(10, 100, by = 10), MinSample = 1, Size = NULL) {

	type <- match.arg(type, c("proportion", "count", "factorBySize", "propBySize"))

	subsets <- list()
	subsets_interm <- list()

	for(np in 1: nperm){ #start permutations
		subsets[[np]] <- list()

		## type propBySize
		if(type=="propBySize"){
			d.intern <- list()
			for(j in prop) {
				d.intern[[j]]<-list()
				for(i in levels(sampUnit)) {
					jj <- round(j*sum(sampUnit==i)/100)
					d.intern[[j]][[i]] <- d[sampUnit==i,][order(Size[sampUnit==i], decreasing = TRUE), ][1:jj,]
				}
				subsets_interm[[j]] <- lapply(d.intern[lapply(d.intern, length)>0], function(x) do.call(rbind, x))
			}
			subsets[[np]] <- subsets_interm[[j]]
		}

		else if (type == "proportion" | type == "count" | type == "factorBySize") {
			stock <- list()

			if (type == "proportion" | type == "factorBySize") {
				val_j <- prop
			}

			if (type == "count") {
				val_j <- c(1 : max(table(sampUnit)))
			}

		    for(j in val_j) {
				subsets_interm[[j]] <- list()
				stock[[j]] <- list()

				for(i in levels(sampUnit)) {
			        stock[[j]][[i]] <- c()

					## type proportion
			        if(type == "proportion"){
						jj <- round(j*sum(sampUnit==i)/100)
						if(jj==0){jj <- MinSample}
						if(jj!=0){
							#stock[[j]][[i]] <- c(sample((1 :nrow(d)) [sampUnit==i], jj))
							toSample <- (1 :nrow(d))[sampUnit==i]
							stock[[j]][[i]]  <- c(toSample[sample.int(length(toSample), size=jj)])
						}
						else(stock[[j]][[i]] <- NA)
			        }

					## type count
			        else if(type == "count"){
						stock[[j]][[i]] <- c()
						if(j <= sum(sampUnit == i)){
							stock[[j]][[i]] <- c(sample((1 :nrow(d)) [sampUnit==i], j))
						}
						else{stock[[j]][[i]] <- c((1 :nrow(d)) [sampUnit==i])}
			        }

			        else if(type=="factorBySize"){
						#calculate the rank of the factor using the SizeFactor
						if(length(Size) == nlevels(sampUnit)){
							rank_fact <- order(Size)
						}
						else {
							rank_fact <- order(tapply(Size, sampUnit, mean, na.rm = TRUE))
						}

						#Choose if we sample this level of factor thanks to the prop argument
						jj <- round(j*nlevels(sampUnit)/100)
						if(jj==0){jj <- MinSample}

						if(rank_fact [levels(sampUnit)==i] <= jj){
							stock[[j]][[i]] <- c((1 :nrow(d)) [sampUnit==i])
						}
					}
				}

				stock[[j]] <- na.omit(unlist(stock[[j]]))
				subsets_interm[[j]] <- d[stock[[j]], ]

				if(type == "proportion" | type=="factorBySize"){
					subsets[[np]] <- subsets_interm[lapply(subsets_interm, length)>0]
				}

				else if(type == "count"){
					subsets[[np]][[j]] <- subsets_interm[[j]]
				}
			}
		}
		print(paste(round(np/nperm, 2)*100, "%"))
	}
	return(subsets)
}
