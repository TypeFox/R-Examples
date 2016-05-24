# 		ASCA (ANOVA-simultaneous component analysis) tool, part of the 'MetStaT' package  
#		Copyright (C) 2012 Gooitzen Zwanenburg and Tim Dorscheidt
#		
#		This program is free software: you can redistribute it and/or modify
#		it under the terms of the GNU General Public License as published by
#		the Free Software Foundation, either version 3 of the License, or
#		(at your option) any later version.
#		
#		This program is distributed in the hope that it will be useful,
#		but WITHOUT ANY WARRANTY; without even the implied warranty of
#		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#		GNU General Public License for more details.
#		
#		You should have received a copy of the GNU General Public License
#		along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# 		Email: g.zwanenburg@uva.nl or tdorscheidt@gmail.com
###############################################################################

##### the first section contains helper functions

# method to determine the power set of an input set
ASCA.GetPowerSet <- function(input.set, exclude.empty.set = FALSE, exclude.complete.set = FALSE) {
	# through recursion, this function (and its nested function) will return all possible permutations of subsets of the inputset
	
	ASCA.PowerSetRecursion <- function(subset) {
		answer <- list()
		if (length(subset) > 0) {
			answer <- c(answer,list(subset))
			for (i in 1:length(subset)) {
				newset <- subset[-i]
				answer<-c(answer,ASCA.PowerSetRecursion(newset))
			}
		}
		answer
	}
	answer <- ASCA.PowerSetRecursion(input.set)
	if (exclude.complete.set) answer[[1]] <- NULL
	if (!exclude.empty.set) answer <- c(list(array()),answer)
	unique(answer)
}

# method to find all unique row patterns in the matrix, and return both these patterns and the row-indices per pattern
ASCA.GetRowRepeats <- function(mat) {
	# this function will find all unique row-patterns in the input, and return both these patterns (in order) and the row-indices where they occur 
	if (!is.matrix(mat)) {mat <- matrix(mat)} 
	result <- list()
	no.rows <- dim(mat)[1]
	no.cols <- dim(mat)[2]
	result$row.patterns <- matrix(nrow=0,ncol=no.cols)
	no.patterns.found <- 0
	result$indices.per.pattern <- list()
	
	ASCA.IsIdenticalIgnoreNames <- function(x, y) { # a new function to check for identical matrices was needed, one that ignores column and row names
		x.nameless <- c(x)
		y.nameless <- c(y)
		if (length(x.nameless)!=length(y.nameless)) {return(FALSE)}
		for (i in 1:length(x.nameless)) {
			if (x.nameless[i]!=y.nameless[i]) {return(FALSE)}
		}
		return(TRUE)
	}
	
	for (r in 1:no.rows) { # go through all the rows in the input matrix, and check whether that row-pattern was already discovered before
		pattern.number <- which(apply(result$row.patterns,1,ASCA.IsIdenticalIgnoreNames,y=mat[r,])==TRUE) # which pattern does this row match, if any
		if (length(pattern.number)==0) { # if the row does not match a previous pattern, then add it to the list of patterns
			result$row.patterns <- rbind(result$row.patterns,mat[r,])
			no.patterns.found <- no.patterns.found + 1
			result$indices.per.pattern[[no.patterns.found]] <- r 
		} else { # the row does match a previous pattern, therefore remember the index of this row as an occurence of that pattern in the matrix
			result$indices.per.pattern[[pattern.number]] <- c(result$indices.per.pattern[[pattern.number]],r)
		}
	}
	levels.order <- sort(as.numeric(apply(result$row.patterns,1,paste,collapse="")),index.return=TRUE)$ix # sort the patterns by numerical order
	result$indices.per.pattern <- result$indices.per.pattern[levels.order]
	result$row.patterns <- result$row.patterns[levels.order,,drop=FALSE]
	result
}

#### the second section contains functions instrinsic to the ASCA algorithm

# Main method of ASCA
ASCA.Calculate <- function(data, levels, equation.elements = "", scaling = FALSE, only.means.matrix = FALSE, use.previous.asca = NULL) {
	# first, a sub-function needs to be defined
	ASCA.GetEquationElement <- function(asca, evaluation, previous.asca) {
		# will calculate the means per unique combination of levels of the factor(s) involved
		# does not reduce this means-matrix with the means-matrix of the individual factors (is done by mother-function performASCA)
		s <- list()
		s$factors.evaluated <- evaluation
		if (!is.null(previous.asca)) {
			s$level.combinations <- previous.asca[[paste(evaluation,collapse="")]]$level.combinations
		} else {
			s$level.combinations <- ASCA.GetRowRepeats(asca$levels[,s$factors.evaluated,drop=FALSE])
		}
		s$means.matrix <- matrix(nrow=dim(asca$data)[1], ncol=dim(asca$data)[2])
		for (p in 1:dim(s$level.combinations$row.patterns)[1]) {
			mean.for.this.level.combination <- colMeans(asca$data[s$level.combinations$indices.per.pattern[[p]], ,drop=FALSE])
			for (i in s$level.combinations$indices.per.pattern[[p]]) {
				s$means.matrix[i,] <- mean.for.this.level.combination
			}
		}
		s
	}
	# this is the main functions for performing an ASCA
	# data must contain a matrix with the sample values
	# levels must contain a column for each factor, and then each row defining which level that sample is part of
	# equation.elements must contain the formula to be solved, but written as follows: 1,2,12 (which will solve for the first and second factor, and their interaction)
	# only.means.matrix can be false in case the calculation needs to be optimized and no SVD nor remainders need to be calculated (mainy used for permutation testing)
	# use.previous.asca can be used to provide an earlier performed ASCA with identical data (but potentially permutated rows in the data), as to skip certain costly steps in the algorithm (mainly used for permutation testing)
	# this function can handle as many factors and interactions as one would like, eg. 1,2,3,12,13,23,123
	s <- list()

	dataAdjusted <- MetStaT.ScalePip(data, center = TRUE, scale=scaling, quietly=TRUE)
	s$ssq.mean <- sum(rep( dataAdjusted$center.vector/dataAdjusted$scale.vector ,nrow(data))^2)
	s$ssq <- sum(data^2)
	s$data <- dataAdjusted$data
	
	if (!is.numeric(levels)) {stop("The supplied levels are not numeric.")}
	s$levels <- levels
	if (!only.means.matrix) {
		# perform svd on original data
		s$svd <- PCA.Calculate(s$data)
	}
	s$ee.names <- c()
	# if no formula was given (equation.elements is empty), then use the powerset of all factors as equation, ergo use all factors and all their interactions
	if (identical(equation.elements,"")) {equation.elements <- ASCA.GetPowerSet(c(1:dim(s$levels)[2]), exclude.empty.set = TRUE)}
	# if equation.elements is not empty, then split input at the commas to obtain the individual equation elements
	if (is.character(equation.elements))  equation.elements <- lapply(strsplit(strsplit(equation.elements,split=",")[[1]],split=""),as.numeric)
	# do some basic checks for correct equation elements, factors, levels, data-matrix sizes
	for (ee in equation.elements) {
		for (f in ee) if (f > dim(levels)[2] || f < 1) {stop(paste("Factor ",f," is beyond scope of study-design",sep=""))}
	}
	if (dim(data)[1] != dim(levels)[1]) {stop(paste("Number of rows in data (",dim(data)[1],") and study design (",dim(levels)[1],") do not match",sep=""))}
	# find the order in which to evaluate the elements of the equation, starting with single factors, then by complexity of interaction
	order.to.evaluate.ee <- sort(as.numeric(unlist(lapply(equation.elements,paste,collapse=""))),index.return=TRUE)$ix
	s$remainder <- s$data # the remainder matrix will eventually be reduced by all found means-matrices of the equation-elements
	for (ee in order.to.evaluate.ee) { # for each element in the equation, ordered by complexity, first derive its specific means-matrix, then (if need be) reduce it
		new.equation.element <- ASCA.GetEquationElement(s, equation.elements[[ee]], use.previous.asca) # this method will return the means-matrix for the specified element
		reductions <- ASCA.GetPowerSet(equation.elements[[ee]], exclude.empty.set = TRUE, exclude.complete.set = TRUE) # reductions of the just derived means-matrix are only necessary if the equation element was an interaction
		for (r in reductions) {
			new.equation.element$means.matrix <- new.equation.element$means.matrix - s[[c(paste(r,collapse=""))]]$means.matrix # reduce the means-matrix by all the sub-factors and sub-interactions that the equation elements holds
		}
		new.equation.element$ssq <- sum(new.equation.element$means.matrix^2)
		if (!only.means.matrix) {
			s$remainder <- s$remainder - new.equation.element$means.matrix # reduce the overall remainder by this means matrix
			new.equation.element$reduced.matrix <- s$remainder # store the current remainder in the new ASCA equation elements list
			# newEquationElement$pca <- prcomp(newEquationElement$meansMatrix)
			
			new.equation.element$svd <- PCA.Calculate(new.equation.element$means.matrix) # perform svd on the new equation element's means matrix
		}	
		ee.name <- paste(equation.elements[[ee]],collapse="")
		s$ee.names <- c(s$ee.names,ee.name)
		s[[ee.name]] <- new.equation.element
	}
	s$ssq.remainder <- sum(s$remainder^2)
	if (!only.means.matrix) ASCA.GetSummary(s)
	s
}

# permutation test for ASCA
ASCA.DoPermutationTest <- function(asca, perm=1000) {
	# test whether purmutated row orderings of the data in the asca model result in mean-matrices with lower ssq-values (sum of squares) then the mean-matrix of the asca-model
	
	# inner function returns the mean per column of the matrix M
	ASCA.GetColumnMeans <- function(M) {
		apply(M,2,mean)
	}
	
	ASCA.PermutateAndGetSSQ <- function(asca, ee, new.row.order) {
		# perform the ASCA again, but with permuated row order. Only perform ASCA for the factor involved, or, in case of an interaction, also the interaction's constituent factors.
		ee.to.evaluate <- ASCA.GetPowerSet(asca[[ee]]$factors.evaluated, exclude.empty.set = TRUE) # in case of an interaction, find its constituent factors (and potentially "smaller" interactions)
		perm.asca <- ASCA.Calculate(asca$data[new.row.order,], asca$levels, equation.elements = ee.to.evaluate, scaling = FALSE, only.means.matrix = TRUE, use.previous.asca = asca) # perform the asca with permutated data-rows
		sum(perm.asca[[ee]]$means.matrix^2) # return the sum of the square of the meansMatrix
	}
	
	# for a certain asca-model, and for a certain equation-element (factor or interaction) only certain permutations are allowed; find one such row-order permutation
	ASCA.GetPermutation <- function(asca, perm.ee) {
		free.factors <- unique(as.numeric(unlist(strsplit(asca$ee.names,split="")))) # get all factors
		perm.factors <- unique(as.numeric(unlist(strsplit(perm.ee,split="")))) # get the factors that need to remain fixed
		for (p in perm.factors) { free.factors <- free.factors[free.factors!=p] } # remove the fixed-factors from the set of all factors, which will now only contain 'free factors', so factors for which the corresponding row order is no longer relevant
		if (length(free.factors)==0) { # if no free factors are left, simply return a permutated row order
			base.row.order = sample(1:dim(asca$data)[1]) 
		} else {
			# What can we permutate now? The 'free factors' indicate which combination of rows exist with previously relevant factor-level information, we can now forget such contraints and permutate freely amongst them
			sets.within.permutations.are.allowed <- asca[[paste(free.factors,collapse="")]]$level.combinations$indices.per.pattern # the asca-model should contain information on which combinations of row-indices are free to permutate on the basis of the left free elements
			base.row.order <- 1:dim(asca$data)[1] # start with the unpermutated row-order
			for (s in sets.within.permutations.are.allowed) { # each of these sets contains a number of row-indices that only differ in the factor-level combinations which are now free, so they can be permutated
				permutated.set <- s[sample(1:length(s))] # permutate the order of indices in this set
				for (i in 1:length(permutated.set)) {
					base.row.order[permutated.set[i]] <- s[i] # store the permutated row-indices of this set
				}
			}
		}
		base.row.order # return the permutated row order
	}
	
	perm.ssq.per.ee <- matrix(ncol=length(asca$ee.names),nrow=perm+1)
	colnames(perm.ssq.per.ee) <- asca$ee.names
	current.col <- 1
	for (ee in asca$ee.names) {
		perm.ssq.per.ee[1,current.col] <- ASCA.PermutateAndGetSSQ(asca, ee, 1:dim(asca$data)[1])
		current.col <- current.col + 1
	}
	for (perm.nr in 2:(perm+1)) {
		current.col <- 1
		for (ee in asca$ee.names) {
			perm.row.order <- ASCA.GetPermutation(asca, ee)
			perm.ssq.per.ee[perm.nr,current.col] <- ASCA.PermutateAndGetSSQ(asca, ee, perm.row.order)
			current.col <- current.col + 1
		}
	}
	p = matrix(ncol=length(asca$ee.names),nrow=1)
	colnames(p) <- asca$ee.names
	for (i in 1:length(asca$ee.names)) {
		p[i] <- (1+sum(perm.ssq.per.ee[-1,i]>=perm.ssq.per.ee[1,i]))/perm # calculate max-limit p value
	}
	p
}

##### the third section contains all functions dealing with plotting of ASCA

# This method will plot all available plots for ASCA
ASCA.Plot <- function(asca) {
	ASCA.PlotScores(asca)
	ASCA.PlotLoadings(asca)
	for (ee in asca$ee.names) {
		ASCA.PlotScoresPerLevel(asca,ee)
		ASCA.PlotLoadings(asca,ee)
	}
}

# plot the scores for a specific factor/interaction
ASCA.PlotScores <- function(asca, ee = "", PCs = "1,2") {
	if (ee == "") {
		pr.object <- asca$svd
		main.title <- "Scoreplot whole data set"
	} else {
		pr.object <- asca[[ee]]$svd
		main.title <- paste("Scoreplot factor-combination",ee)
	}
	generate.plot <- TRUE
	if (generate.plot) {
		PCs <- as.numeric(strsplit(PCs,split=",")[[1]])
		# externalDevice <- NULL
		# if (dev.cur()=="1") {
		#   dev.new(width=8.27, height=11.69)
		# }
		# par(mfrow=c(2,1))
		plot(pr.object$t[,PCs[1]],pr.object$t[,PCs[2]],
				main=main.title,
				xlab=paste("PC",PCs[1]," (",formatC(pr.object$var.explained[PCs[1]] * 100,digits=2,format="f"),"%)",sep=""),
				ylab=paste("PC",PCs[2]," (",formatC(pr.object$var.explained[PCs[2]] * 100,digits=2,format="f"),"%)",sep=""));
	}
	
}

# plot the loadings of a specific factor/interaction
ASCA.PlotLoadings <- function(asca, ee= "", pcs = c(1,2)) {
	if (ee == "") {
		pr.object <- asca$svd
		main.title <- "Loadings whole data set"
	} else {
		pr.object <- asca[[ee]]$svd
		main.title <- paste("Loadings factor-combination",ee)
	}
	if (is.character(pcs)) eval(parse(text=paste(sep="","pcs <- c(",pcs,")")))
	
	list.of.pc.tuples <- MetStaT.GetPcTuples(length(pcs)) # get tuples of PC-indices that are plotted versus one another
	tuple.index <- 0
	prev.mfrow <- par("mfrow") # save current plot-window dimensions to restore this back to default after we're done
	prev.xpd <- par("xpd")# save current plot-window clipping parameter to restore this back to default after we're done
	prev.oma <- par("oma")# save current plot-figure margin parameters to restore this back to default after we're done
	prev.mar <- par("mar")# save current plot-area margin parameters to restore this back to default after we're done
	par(oma=c(0,0,0,0)) # reduce figure margins to maximize drawing area 
	par(mfrow=c(2,1)) # we'll be plotting two seperate loading-plots above one another
	par(xpd=NA) # diasable clipping
	for (tuple in list.of.pc.tuples) {
		tuple.index <- tuple.index + 1
		pc1 <- pcs[tuple[1]] # 1st element of tuple contains index of PC1 within PCs
		pc2 <- pcs[tuple[2]] # 2nd element of tuple contains index of PC2 within PCs
		
		# top plot
		par(mar=c(0,4.1,4.1,2.1)) # x-axis will only contain ticks, remove lower margins to the point where plots are touching
		main.title <- main.title
		ylab.title <- paste("Distance to PC",pc1,sep="") # y axis label
		color.to.Use <- pc1
		plot(1:dim(pr.object$v)[1],pr.object$v[,pc1],type="h",col=color.to.Use,
				main=main.title,
				axes=FALSE, # manually make axes
				xlab="",
				ylab=ylab.title,
				ylim= c(
						min(c(pr.object$v[,pc1],pr.object$v[,pc2])),
						max(c(pr.object$v[,pc1],pr.object$v[,pc2])) 
				) # limit y-axis to min and max of either pc
		)
		axis(1, labels=FALSE)
		axis(2, labels=TRUE)
		
		# bottom plot
		par(mar=c(5.1,4.1,0,2.1)) # remove upper margins to the point where plots are touching
		ylab.title <- paste("Distance to PC",pc2,sep="") # y axis label
		color.to.Use <- pc2
		plot(1:dim(pr.object$v)[1],pr.object$v[,pc2],type="h",col=color.to.Use,
				main=NULL,
				axes=FALSE, # manually make axes
				xlab="Variables",
				ylab=ylab.title,
				ylim= c(
						min(c(pr.object$v[,pc1],pr.object$v[,pc2])),
						max(c(pr.object$v[,pc1],pr.object$v[,pc2])) 
				) # limit y-axis to min and max of either pc
		)
		axis(1, labels=TRUE)
		axis(2, labels=TRUE)
	}
	par(mfrow=prev.mfrow) # restore plot window dimensions back to previous
	par(xpd=prev.xpd) # restore clipping parameter back to previous
	par(oma=prev.oma) # restore figure margin parameters back to previous
	par(mar=prev.mar) # restore plot margin parameters back to previous
}

# plot the scores and distinguish between the levels
ASCA.PlotScoresPerLevel <- function(asca, ee, pcs = "1,2") {
	pcs <- as.numeric(strsplit(pcs,split=",")[[1]])
	y <- (asca[[ee]]$means.matrix + asca$remainder) %*% asca[[ee]]$svd$v
	t.list.x <- list()
	t.list.y <- list()
	list.color.pattern <- list()
	color.per.variable <- rep(0,dim(asca$data)[1])
	pattern.per.variable <- rep(0,dim(asca$data)[1])
	kColOptions <- c(24,552,254,26,84,51,652,68,76,96,10,60,33,   245,147,12,26,164,181,52,512,344,201,111)
	kPointOptions <- 1:30
	for (p in 1:dim(asca[[ee]]$level.combinations$row.pattern)[1]) {
		if (length(asca[[ee]]$level.combinations$row.pattern[p,])==1) {
			list.color.pattern[[p]] <- c(kColOptions[p],kPointOptions[p])
		} else if (length(asca[[ee]]$level.combinations$row.pattern[p,])==2) {
			list.color.pattern[[p]] <- c(
					kColOptions[asca[[ee]]$level.combinations$row.pattern[p,1]],
					kPointOptions[asca[[ee]]$level.combinations$row.pattern[p,2]]
			)
		}else {
			list.color.pattern[[p]] <- c(
					kColOptions[asca[[ee]]$level.combinations$row.pattern[p,1]]%%9,
					floor(kPointOptions[asca[[ee]]$level.combinations$row.pattern[p,2]]/9)
			)
		}
		color.per.variable[asca[[ee]]$level.combinations$indices.per.pattern[[p]]] <- list.color.pattern[[p]][1]
		pattern.per.variable[asca[[ee]]$level.combinations$indices.per.pattern[[p]]] <- list.color.pattern[[p]][2]
		t.list.x[[p]] <- y[asca[[ee]]$level.combinations$indices.per.pattern[[p]],pcs[1]]
		t.list.y[[p]] <- y[asca[[ee]]$level.combinations$indices.per.pattern[[p]],pcs[2]]
	}
	
	legend.colors.patterns <- do.call(rbind,list.color.pattern)
	plot(asca[[ee]]$svd$t[,pcs[1]],asca[[ee]]$svd$t[,pcs[2]],
			xlim = range(c(min(unlist(t.list.x)),max(unlist(t.list.x)))), ylim = range(c(min(unlist(t.list.y)),max(unlist(t.list.y)))),
			main=paste("PC",pcs[1]," vs PC",pcs[2]," for factor-combination ",ee,sep=""),
			xlab=paste("PC",pcs[1]," (",formatC(asca[[ee]]$svd$var.explained[pcs[1]] * 100,digits=2,format="f"),"%)",sep=""),
			ylab=paste("PC",pcs[2]," (",formatC(asca[[ee]]$svd$var.explained[pcs[2]] * 100,digits=2,format="f"),"%)",sep=""),
			cex=1.5, lwd=3,
			col=colors()[color.per.variable], pch=pattern.per.variable)
	
	legend(x = "bottomright",
			apply(asca[[ee]]$level.combinations$row.patterns,1,paste,collapse=" "),
			cex=0.8, 
			col=colors()[legend.colors.patterns[,1]], pch=legend.colors.patterns[,2]);
	
	for (p in 1:length(t.list.x)) {
		points(t.list.x[[p]],t.list.y[[p]],col=colors()[list.color.pattern[[p]][1]],pch=list.color.pattern[[p]][2])
	}
#  TlistX
}

# give a table with summary statistics for ASCA
ASCA.GetSummary <- function(asca, quietly = FALSE) {
	summaryPCA <- matrix(NA, nrow = length(asca$ee.names) + 1, ncol=10)
	rownames(summaryPCA) <- c("data",asca$ee.names)
	colnames(summaryPCA) <- paste("PC",1:10,sep="")
	relevant.pcs <- which(asca$svd$var.explained[1:10]>0.01) # only display the first 10 PCs which are above 1% variance explained
	max.pc <- 0
	for (i in relevant.pcs) {
		max.pc <- max(max.pc, i)
		summaryPCA[1,i] <- asca$svd$var.explained[i]
	}
	ee.index <- 2
	for (ee in asca$ee.names) {
		relevant.pcs <- which(asca[[ee]]$svd$var.explained[1:10]>0.01) # only display the first 10 PCs which are above 1% variance explained
		for (i in relevant.pcs) {
			max.pc <- max(max.pc, i)
			summaryPCA[ee.index,i] <- asca[[ee]]$svd$var.explained[i]
		}
		ee.index <- ee.index + 1
	}
	summaryPCA <- summaryPCA[,1:max.pc,drop=FALSE]
	
	summarySSQ <- matrix(NA, nrow = 1, ncol = length(asca$ee.names) + 2)
	rownames(summarySSQ) <- "Contribution to ssq"
	colnames(summarySSQ) <- c("Overall means", asca$ee.names, "Residuals")
	summarySSQ[1,1] <- asca$ssq.mean/asca$ssq
	ee.index <- 2
	for (ee in asca$ee.names) {
		summarySSQ[1,ee.index] <- asca[[ee]]$ssq/asca$ssq
		ee.index <- ee.index + 1
	}
	summarySSQ[1,dim(summarySSQ)[2]] <- asca$ssq.remainder/asca$ssq
	
	summarySSQCentered <- matrix(NA, nrow = 1, ncol = length(asca$ee.names) + 1)
	rownames(summarySSQCentered) <- "Contribution to ssq"
	colnames(summarySSQCentered) <- c(asca$ee.names, "Residuals")
	ee.index <- 1
	for (ee in asca$ee.names) {
	    summarySSQCentered[1,ee.index] <- asca[[ee]]$ssq/(asca$ssq - asca$ssq.mean)
	    ee.index <- ee.index + 1
	}
	summarySSQCentered[1,dim(summarySSQCentered)[2]] <- asca$ssq.remainder/(asca$ssq - asca$ssq.mean)
    
	if (!quietly) {
		cat("Variance explained per principal component (if >1%):\n")
		for (r in 1:dim(summaryPCA)[1]) {
			if (r==1) { 									cat(paste("Whole data set \t",sep="")) 
			} else if (nchar(rownames(summaryPCA)[r])==1) {	cat(paste("Factor ",rownames(summaryPCA)[r],"     \t",sep=""))
			} else {										cat(paste("Interaction ",rownames(summaryPCA)[r],"\t",sep=""))}
			for (c in 1:dim(summaryPCA)[2]) {
				varExplained <- paste(colnames(summaryPCA)[c],": ",formatC(summaryPCA[r,c] * 100,digits=2,format="f"),"% ",sep="")
				varExplained <- paste(varExplained,paste(rep(" ",(13-nchar(varExplained))),collapse=""))
				cat(varExplained)
			}
			cat("\n")
		}
		cat("\n")
		cat("Percentage each effect contributes to the total sum of squares:\n")
		for (c in 1:dim(summarySSQ)[2]) {
			if (c==1) { 									cat(paste("Overall means  \t",sep=""))
			} else if (c==dim(summarySSQ)[2]) {				cat(paste("Residuals      \t",sep=""))
			} else if (nchar(colnames(summarySSQ)[c])==1) {	cat(paste("Factor ",colnames(summarySSQ)[c],"     \t",sep=""))
			} else {										cat(paste("Interaction ",colnames(summarySSQ)[c],"\t",sep=""))}
			cat(paste(formatC(summarySSQ[1,c] * 100,digits=2,format="f"),"%\n",sep=""))
		}
		cat("\n")
        
        
		cat("Percentage each effect contributes to the sum of squares of the centered data:\n")
		for (c in 1:dim(summarySSQCentered)[2]) {
		    if (c==dim(summarySSQCentered)[2]) {				cat(paste("Residuals      \t",sep=""))
		    } else if (nchar(colnames(summarySSQCentered)[c])==1) {	cat(paste("Factor ",colnames(summarySSQCentered)[c],"     \t",sep=""))
		    } else {										cat(paste("Interaction ",colnames(summarySSQCentered)[c],"\t",sep=""))}
		    cat(paste(formatC(summarySSQCentered[1,c] * 100,digits=2,format="f"),"%\n",sep=""))
		}
		cat("\n")
	}
	
	list(summary.pca = summaryPCA, summary.ssq = summarySSQ)
}
