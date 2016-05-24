## implement Tariff method and calculate accuracy metrics
## algorithm last update 2015.03.13
## wrapper last update 2015.12.06

#' Replicate Tariff methods
#'
#' This function implements Tariff method.
#'
#' @param causes.train character vector of causes, or the column name of cause in the training data
#' @param symps.train N.train by S matrix 
#' @param symps.test N.test by S matrix 
#' @param causes.table list of causes in the data
#' @param use.rank logical indicator for whether using ranks instead of scores
#' @param nboot.rank number of re-sampling for baseline rank comparison
#' @param use.sig logical indicator for whether using significant Tariff only
#' @param nboot.sig  number of re-sampling for testing significance.
#' @param use.top logical indicator for whether the tariff matrix should be cleaned to have only top symptoms
#' @param ntop number of top tariff kept for each cause
#' @param ... not used
#' 
#' @return \item{score}{matrix of score for each cause within each death}
#' \item{causes.train}{vector of most likely causes in training data} 
#' \item{causes.test}{vector of most likely causes in testing data} 
#' \item{csmf}{vector of CSMF}
#' \item{causes.table}{cause list used for output, i.e., list of existing causes in the training data}
#' @author Zehang Li, Tyler McCormick, Sam Clark
#' 
#' Maintainer: Zehang Li <lizehang@@uw.edu>
#' @references James, S. L., Flaxman, A. D., Murray, C. J., & Population Health Metrics Research Consortium. (2011). \emph{Performance of the Tariff Method: validation of a simple additive algorithm for analysis of verbal autopsies.} \emph{Population Health Metrics, 9(1), 1-16.}
#' @references Serina, P., Riley, I., Stewart, A., James, S. L., Flaxman, A. D., Lozano, R., ... & Ahuja, R. (2015). \emph{Improving performance of the Tariff Method for assigning causes of death to verbal autopsies.} \emph{BMC medicine, 13(1), 1.}
#' @references Tyler H. McCormick, Zehang R. Li, Clara Calvert, Amelia C.
#' Crampin, Kathleen Kahn and Samuel J. Clark(2016) \emph{Probabilistic
#' cause-of-death assignment using verbal autopsies},
#' \url{http://arxiv.org/abs/1411.3042} \emph{To appear, Journal of the American Statistical Association}
#' @keywords Tariff
#' @examples
#'\donttest{
#' data("RandomVA3")
#' test <- RandomVA3[1:200, ]
#' train <- RandomVA3[201:400, ]
#' allcauses <- unique(train$cause)
#' fit <- tariff(causes.train = "cause", symps.train = train, 
#' 				symps.test = test, causes.table = allcauses)
#' correct <- which(fit$causes.test[,2] == test$cause)
#' accuracy <- length(correct) / dim(test)[1]
#' }

tariff <- function(causes.train, symps.train, symps.test, causes.table = NULL,  use.rank = TRUE, nboot.rank = 1, use.sig = TRUE, nboot.sig = 500, use.top = FALSE, ntop = 40, ...){
	
	
	# if input cause is the column name
	if(class(causes.train) == "character" && length(causes.train) == 1){
		colindex <- match(causes.train, colnames(symps.train))
		colindex2 <- match(causes.train, colnames(symps.test))
    
		if(is.na(colindex)){
		  stop("Cannot find the cause-of-death column in training data")
		}
		causes.train <- symps.train[, colindex]
		symps.train <- symps.train[, -colindex]

		# also remove this from testing data if it is provided
		if(!is.na(colindex2)){
			causes.test <- symps.test[, colindex2]
			symps.test <- symps.test[, -colindex2]		
		}
	}

  if(is.null(causes.table)){
    causes.table <- unique(causes.train)
  }
  
	id.train <- symps.train[, 1]
	symps.train <- symps.train[, -1]
	id.test <- symps.test[, 1]
	symps.test <- symps.test[, -1]
	# make sure train and test has the same columns
	joint <- intersect(colnames(symps.train), colnames(symps.test))
	symps.test <- symps.test[, joint]
	symps.train <- symps.train[, joint]

	# function to convert InterVA input back to numeric
	# @para
	#	mat     : matrix containing "Y", "" and "."
	#   missing : impute value for "."
	toBinary <- function(mat, missing = 0){
		mat <- as.matrix(mat)
		mat2 <- matrix(0, dim(mat)[1], dim(mat)[2])
		mat2[which(toupper(mat) == "Y")] <- 1
		mat2[which(mat == ".")] <- missing
		return(mat2)
	} 

	# steps for re-sampling
	# function to count combinations from symps and causes
	# @para
	#	symps       : N by S matrix
	#   causes      : N matrix
	#   causelist   : C vector
	#   binary      : boolean indicating "Y" or 1
	# @return
	#    count matrix (C by S)
	count.combo <- function(symps, causes, causelist, binary){
		cond.count <- matrix(0, length(causelist), dim(symps)[2])
	 	for(i in 1:length(causelist)){
			cause <- causelist[i]
			list <- which(causes == cause)
			if(length(list) == 0){
				# cat("!")
				next;
			}
			cases <- symps[list, ,drop=FALSE]

			if(length(list) == 1){

			}
			# remove missing from calculation
			if(binary){
				count <- apply(cases, 2, function(x){length(which(x == 1))})
			}else{
				count <- apply(cases, 2, function(x){length(which(x == "Y"))})
			}
			cond.count[i, ] <-  count
		}
		return(cond.count)
	}

	# function to get tariff (not tariff score) from count matrix
	# @para
	#	X: C by S matrix
	# @return
	#   tariff: C by S matrix
	getTariff <- function(X){
		# find median by symptom
		med <- apply(X, 2, median)
		# find IQR by symptom, handle zero denominator
		iqr <- apply(X, 2, IQR)
		# if iqr is 0, replace with range
		rr <- apply(X, 2, function(x){max(x) - min(x)})	
		iqr[which(iqr == 0)] <- rr[which(iqr == 0)]
		# if range is also 0, means all the same value, 
		#	since it will be 0 anyway, just to avoid division by 0 here
		iqr[which(iqr == 0)] <- 0.05
		# rescale to Tariff now (still C by S matrix)
		# NOTICE < danger >:
		# R perform matrix - vector by column, so it has to be transposed first
		tariff <- (t(X) - med) / iqr

		tariff <- t(tariff)
		# Tariff paper suggests the rounding here, to avoid over-fitting?
		tariff <- round(2 * tariff) / 2
		return(tariff)
	}

	# function to convert score to rank compared only to Resampled training set
	# @para
	#	mat     : score matrix (C by N)
	#   all     : all score matrix (K by C)
	# return
	#	mat.rank: score rank matrix (C by N)
	toRank <- function(mat, all){
		cat("Calculating ranks\n")
		N <- dim(mat)[2]
		C <- dim(mat)[1]

		# notice rank from small to large, so take negative
		mat.rank <- lapply(seq(1:C), 
			function(kk){
				# cat("."); 
				return(
					sapply(mat[kk, ], function(x, y){rank(-c(x, y))[1]}, all[,kk])
			          )})
		out <- mat
		for(i in 1:C){
			out[i, ] <- mat.rank[[i]]
		}
		return(out)
	} 

	# find individual top n causes  
	# function to get the n-th largest element's index
	which.ordern <- function(x, n, isMax){
		return(order(x, decreasing = isMax)[n])
	}



	##
	## starts algorithm
	##
	
	# remove causes not in the training data
	nonexist <- which(causes.table %in% unique(causes.train) == FALSE)
	if(length(nonexist) > 0){
		causes.table2 <- causes.table[-nonexist]
	}else{
		causes.table2 <- causes.table
	}
	S <- dim(symps.train)[2]
	C <- length(causes.table2)
	N.train <- dim(symps.train)[1]
	N.test <- dim(symps.test)[1]
	symps.num <- toBinary(symps.train, missing = 0)

	##################################################################
	# first bootstrap step, removing insignificant cause-symptom combo 
	if(use.sig){
		cat("\nStart re-sampling for significant Tariff cells\n")
		all.tariff.boot <- array(0, dim = c(nboot.sig, C, S))
		# all.tariff.score.boot <- matrix(0, N.boot*N.train, C)
		for(i in 1 : nboot.sig){
			sample.boot <- sample(1:N.train, size = N.train, replace = TRUE)
			symps.boot <- symps.num[sample.boot, ]
			cause.boot <- causes.train[sample.boot]
			count.boot <- count.combo(symps.boot, cause.boot, causes.table2, binary=T)
			all.tariff.boot[i, , ] <- getTariff(count.boot)
			# if(i %% 10 == 0) cat(".")
		}
		# get the lower bound for tariff, C by S matrix
		lower <- apply(all.tariff.boot, c(2,3), function(x){quantile(x, 0.025)})
		# get the upper bound for tariff, C by S matrix
		upper <- apply(all.tariff.boot, c(2,3), function(x){quantile(x, 1-0.025)})
		# check if it covers zero
		cover <- sign(lower * upper)
		# define the indicator matrix of which tariff to be removed
		insig <- cover * 0
		insig[which(cover != -1)] <- 1
		if(sum(insig) == 0){
			warning("No Tariff is significant, remove bootstrapping step")
			insig <- matrix(1, C, S)
		}

	}else{
		insig <- matrix(1, C, S)
	}

	###########################################################################
	# calculate Tariff from training data
	# function to remove tail tariff values
	cleanTariff <- function(tariff, nonzero){
		for(i in 1:dim(tariff)[1]){
			order.tmp <- order(abs(tariff[i, ]), decreasing = TRUE)
			tariff[i, order.tmp[1:nonzero]] <- 0
		}
		return(tariff)
	}
	# calculate actual tariff and delete the insignificant combo
	X.train <- count.combo(symps.num, causes.train, causes.table2, binary=T)
	tariff <- getTariff(X.train) * insig

	if(use.top){
		tariff <- cleanTariff(tariff, ntop)
	}

	###########################################################################
	# second bootstrap step, getting Tariff score dist using uniform cause dist 
	# factor of each resampling draw. Truncate N.train to multiple of C
	if(use.rank){
		factor <- trunc(N.train / C)
		all.score.boot <- matrix(0, nboot.rank * factor * C, C)

		# calculate which deaths are by which cause
		index.by.cause <- lapply(causes.table2, function(k){which(causes.train == k)})
		for(i in 1:nboot.rank){
			sample.boot <- rep(0, C*factor)
			# re-sample stratified by cause
			for(j in 1:C){
				sample.boot[((j-1)*factor + 1):(j * factor)] <- sample(index.by.cause[[j]]
					, factor, replace = TRUE)
			}
			symps.boot <- symps.num[sample.boot, ]
			cause.boot <- causes.train[sample.boot]
			# count.boot <- count.combo(symps.boot, cause.boot, causes.table, binary=T)
			# tariff.boot <- getTariff(count.boot) * insig
			# score.boot <- tariff.boot %*% t(symps.boot)
			score.boot <- tariff %*% t(symps.boot)
			all.score.boot[((i-1)*factor*C + 1) : (i*factor*C), ] <- t(score.boot)
			# if(i %% 10 == 0) cat(".")
		}	
	}


	# convert InterVA input back into numeric form (N by S matrix)
	# test.blow <- alldata$symps.test[sample()]
	symps.num.test <- toBinary(symps.test, missing = 0) 
	# calculate tariff score (C by N matrix)
	score.num <- tariff %*% t(symps.num.test) 
	# find individual top cause (max score) 
	# use the cause names instead of the indexes
	causes.test <- causes.table2[apply(score.num, 2, which.max)]
	#
	score <- NULL

	if(use.rank){
	 	# find individual top cause (min rank) 
	 	score <- toRank(score.num, all.score.boot)
	 	causes.test <- causes.table2[apply(score, 2, which.min)]
	}else{
		score <- score.num
	}
	colnames(score) <- id.test
	rownames(score) <- causes.table2

	# find CSMF for testing set
	CSMF <- (table(c(causes.test, causes.table2)) - 1) / length(causes.test)
	if(use.rank){
		CSMF <- (table(c(causes.test, causes.table2)) - 1)/length(causes.test)
	}
	names(CSMF) <- causes.table2


	# might need to include a way to transform from causes.table2 back to causes.table here
	causes.train.out <- data.frame(ID = id.train)
	causes.train.out$cause <- as.character(causes.train)
	causes.test.out <- data.frame(ID = id.test)
	causes.test.out$cause <- as.character(causes.test)

	fit <- list(score = t(score),
				causes.train = causes.train.out,
				causes.test = causes.test.out,
				csmf = CSMF, 
				causes.table = causes.table2)
	class(fit) <- "tariff"
	return(fit)
}


