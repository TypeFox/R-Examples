#' Get individual COD probabilities from InSilicoVA Model Fits
#' 
#' This function calculates individual probabilities for each death and provide posterior credible intervals for each estimates. The default set up is to calculate the 95% C.I. when running the inSilicoVA model. If a different C.I. is desired after running the model, this function provides faster re-calculation without refitting the model.
#' 
#' 
#' @param object Fitted \code{"insilico"} object.
#' @param data data for the fitted \code{"insilico"} object. The first column of the data should be the ID that matches the \code{"insilico"} fitted model.
#' @param CI Credible interval for posterior estimates.
#' @param is.aggregate logical indicator for constructing aggregated distribution rather than individual distributions.
#' @param by list of column names to group by.
#' @param java_option Option to initialize java JVM. Default to ``-Xmx1g'', which sets the maximum heap size to be 1GB.
#' @param \dots Not used.
#' @return
#' \item{mean}{ individual mean COD distribution matrix.} 
#' \item{median}{ individual median COD distribution matrix.} 
#' \item{lower}{ individual lower bound for each COD probability.} 
#' \item{upper}{ individual upper bound for each COD probability.} 

#' @author Zehang Li, Tyler McCormick, Sam Clark
#' 
#' Maintainer: Zehang Li <lizehang@@uw.edu>
#' @seealso \code{\link{insilico}}, \code{\link{updateIndiv}}, \code{\link{plot.insilico}}
#' @references Tyler H. McCormick, Zehang R. Li, Clara Calvert, Amelia C.
#' Crampin, Kathleen Kahn and Samuel J. Clark Probabilistic cause-of-death
#' assignment using verbal autopsies, \emph{arXiv preprint arXiv:1411.3042}
#' \url{http://arxiv.org/abs/1411.3042} (2014)
#' @examples
#' \dontrun{
#' data(RandomVA1)
#' fit1<- insilico(RandomVA1, subpop = NULL,  
#'                 Nsim = 1000, burnin = 500, thin = 10 , seed = 1,
#'                 auto.length = FALSE)
#' summary(fit1, id = "d199")
#' 
#' # Calculate aggregated COD distributions
#' agg.csmf <- get.indiv(data = RandomVA2, fit1, CI = 0.95, 
#'                      is.aggregate = TRUE, by = NULL)
#' head(agg.csmf) 
#' 
#' agg.by.sex.age <- get.indiv(data = RandomVA2, fit1, CI = 0.95, 
#'                             is.aggregate = TRUE, by = list("sex", "age"))
#'							   head(agg.by.sex.age$mean) 
#' }
#' @export get.indiv
get.indiv <- function(object, data = NULL, CI = 0.95, is.aggregate = FALSE, by = NULL, java_option = "-Xmx1g", ...){
	if(is.null(java_option)) java_option = "-Xmx1g"
	options( java.parameters = java_option )
	id <- object$id
	obj <- .jnew("sampler/InsilicoSampler2")

	#-----------------------------------------------------------------#
	# if grouping is provided
	if(is.aggregate){
		if(is.null(data)){
			stop("Data not provided for grouping")
		}
		data <- data.frame(data)
		# check by variable
		if(is.null(by)){
			cat("No groups specified, aggregate for all deaths\n")
			data$sub <- "All"
			col.index <- which(colnames(data) == "sub")
		}else{
			if(class(by) != "list"){
				stop("Please specify by variable using list, see manual example for details.")
			}else{
				col.index <- match(by, colnames(data))
				if(length(which(is.na(col.index))) > 0){
					stop("Contains unknown grouping variables!")
				}
			}
		}
		# check id
		if(colnames(data)[1] != "ID"){
			colnames(data)[1] <- "ID"
			warning("The first column of data is assumed to be ID")
		}
		n.match <- length(which(!is.na(match(data[, 1], object$id)))) 
		if(length(n.match) == 0){
			stop("No data has the matching ID as in the fitted insilico object")
		}else{
			cat(paste(n.match, "deaths matched from data to the fitted object\n"))
		}

		# get grouping vector and names for all data
		datagroup.all <- data[, c(1, col.index)]
		if(length(col.index) > 1){
			datagroup.all$final.group<- apply(datagroup.all[, -1], 1, function(x){paste(x,collapse = ' ')})
		}else{
			datagroup.all$final.group <- datagroup.all[, -1]
		}
		# get the group name for subset in insilico fitted data
		index.tmp <- match(rownames(object$data), datagroup.all[, 1])
		if(length(which(!is.na(index.tmp))) == 0){
			stop("no matching ID found in the data")
		}
		datagroup <- datagroup.all$final.group[index.tmp]
		datagroup <- datagroup[!is.na(datagroup)]

		allgroups <- sort(unique(datagroup))
		datagroup <- match(datagroup, allgroups)
		datagroup.all$final.group.num <- match(datagroup.all$final.group, allgroups)
		datagroup.j <- .jarray(as.integer(datagroup), dispatch = TRUE)
		Ngroup <- as.integer(length(allgroups))
	}

	#---------------------------------------------------------------#

	# 3d array of csmf to be Nsub * Nitr * C
	if(is.null(object$subpop)){
		csmf <- array(0, dim = c(1, dim(object$csmf)[1], dim(object$csmf)[2]))
		csmf[1, , ] <- object$csmf
		subpop <- rep(as.integer(1), length(object$id))
	}else{
		csmf <- array(0, dim = c(length(object$csmf), dim(object$csmf[[1]])[1], dim(object$csmf[[1]])[2]))
		for(i in 1:length(object$csmf)){
			csmf[i, , ] <- object$csmf[[i]]
		}
		subpop <- as.integer(match(object$subpop, names(object$csmf)))
	}	


	if(object$external){
		csmf <- csmf[, , -object$external.causes, drop = FALSE]
		csmf <- apply(csmf, c(1, 2), function(x){
			if(sum(x) != 0){
				return(x/sum(x))
			}else{
				return(x)
			}})
		csmf <- aperm(csmf, c(2, 3, 1))
	}

	data.j <- .jarray(object$data, dispatch = TRUE)
	subpop.j <- .jarray(subpop, dispatch = TRUE)
	impossible.j <- .jarray(object$impossible.causes, dispatch = TRUE)
	# get condprob to be Nitr * S * C array
	if(object$updateCondProb == FALSE){
		condprob <- array(0, dim = c(dim(csmf)[2], dim(object$probbase)[1], dim(object$probbase)[2]))
		for(itr in 1:dim(csmf)[2]){
			condprob[itr, , ] <- object$probbase
		}
	}else{
		if(object$keepProbbase.level){
			condprob <- array(0, dim = c(dim(object$conditional.probs)[1], dim(object$probbase)[1], dim(object$probbase)[2]))
			for(i in 1:dim(condprob)[1]){
				#fix for interVA probbase
				object$probbase[object$probbase == "B -"] <- "B-"
				object$probbase[object$probbase == ""] <- "N"
				
				temp <- object$conditional.probs[i, match(object$probbase, colnames(object$conditional.probs))]
				condprob[i, , ] <- matrix(as.numeric(temp), 
											dim(object$probbase)[1], 
										    dim(object$probbase)[2])
			}
		}else{
			condprob <- object$conditional.probs
		}
	}

	# rJava package (0.9-8) currently has a bug with multi-dimension array
	# value passed to java contains 0's not in the original array
	# to bypass this open issue, condprob and csmf are transformed 
	# 	into vectors:
	# dim(csmf) = Nsub * Nitr * C 
	# dim(condprob) = Nitr * S * C
	# after transformation (i, j, k) -> (k-1)*d1*d2 + (j-1)*d1 + i
	Nsub <- as.integer(dim(csmf)[1])
	Nitr <- as.integer(dim(csmf)[2])
	C <- as.integer(dim(csmf)[3])
	S <- as.integer(dim(condprob)[2])

	condprob <- as.vector(condprob)
	csmf <- as.vector(csmf)
	condprob.j <- .jarray(condprob, dispatch = TRUE)
	csmf.j <- .jarray(csmf, dispatch = TRUE)
	#-------------------------------------------------------------------#
	if(!is.aggregate){
		cat("Calculating individual COD distributions...\n")
		# return rbind(mean, median, low, up), (4N) * C matrix
		indiv  <- .jcall(obj, "[[D", "IndivProb", 
					 data.j, impossible.j, csmf.j, subpop.j, condprob.j, 
					 (1 - CI)/2, 1-(1-CI)/2, 
					 Nsub, Nitr, C, S)
		indiv <- do.call(rbind, lapply(indiv, .jevalArray))
		
		# if not customized probbase, use longer cause names for output
		if(!object$is.customized){
			data("causetext", envir = environment())
			causetext<- get("causetext", envir  = environment())
			
			match.cause <- pmatch(causetext[, 1],  colnames(object$probbase))
			index.cause <- order(match.cause)[1:sum(!is.na(match.cause))]
			colnames(indiv) <- causetext[index.cause, 2]
		}else{
			colnames(indiv) <- colnames(object$probbase)
		}

		K <- dim(indiv)[1] / 4


		## add back all external cause death 41:51 in standard VA
		if(object$external){
			external.causes <- object$external.causes
			C0 <- dim(indiv)[2]
			ext.flag <- apply(object$indiv.prob[, external.causes], 1, sum)
			ext.probs <- object$indiv.prob[which(ext.flag == 1), ]

			indiv <- cbind(indiv[, 1:(external.causes[1] - 1)], 
				          matrix(0, dim(indiv)[1], length(external.causes)), 
				          indiv[, external.causes[1]:C0])
			colnames(indiv) <- colnames(object$indiv.prob)
			id.out <- c(id[match(rownames(object$data), id)], id[which(ext.flag == 1)])
		}else{
			id.out <- id
			ext.probs <- NULL
		}

		mean <- rbind(indiv[1:K, ], ext.probs)
		median <- rbind(indiv[(K+1):(2*K), ], ext.probs)
		lower <- rbind(indiv[(2*K+1):(3*K), ], ext.probs)
		upper <- rbind(indiv[(3*K+1):(4*K), ], ext.probs)

		rownames(mean) <- rownames(median) <- rownames(lower) <- rownames(upper) <- id.out
		return(list(mean = mean, median = median, 
					lower = lower, upper = upper))

	#-------------------------------------------------------------------#
	}else{
		cat("Aggregating individual COD distributions...\n")
		# return (4 x Ngroup) * (C + 1) matrix, last column is sample size
		indiv  <- .jcall(obj, "[[D", "AggIndivProb", 
					 data.j, impossible.j, csmf.j, subpop.j, condprob.j,
					 datagroup.j, Ngroup, (1 - CI)/2, 1-(1-CI)/2, 
					 Nsub, Nitr, C, S)
		indiv <- do.call(rbind, lapply(indiv, .jevalArray))
		data("causetext", envir = environment())
		causetext<- get("causetext", envir  = environment())
		
		weight <- indiv[, dim(indiv)[2]]
		indiv <- indiv[, -dim(indiv)[2]]

		match.cause <- pmatch(causetext[, 1],  colnames(object$probbase))
		index.cause <- order(match.cause)[1:sum(!is.na(match.cause))]
		colnames(indiv) <- causetext[index.cause, 2]

		K <- dim(indiv)[1] / 4


		## add back all external cause death 41:51 in standard VA
		if(object$external){
			external.causes <- object$external.causes
			C0 <- dim(indiv)[2]
			ext.flag <- apply(object$indiv.prob[, external.causes], 1, sum)
			ext.probs <- data.frame(object$indiv.prob[which(ext.flag == 1), external.causes])
			ext.probs$ID <- rownames(ext.probs)
			ext.probs <- merge(ext.probs, datagroup.all[, c("ID", "final.group.num")])
			probcolumns <- which(colnames(ext.probs) %in% c("ID", "final.group.num") == FALSE)

			ext.probs.agg <- matrix(0, length(allgroups), length(probcolumns))
			ext.weight <- rep(0, length(allgroups))

			for(i in 1:length(allgroups)){
				thisgroup <- which(ext.probs$final.group.num == i)
				if(length(thisgroup) > 0){
					ext.probs.agg[i, ] <- apply(ext.probs[thisgroup, probcolumns, drop=FALSE], 2, mean)	
					
					# get weight using the count instead of ratio
					# easier for later manipulation when the total count is divided
					ext.weight[i] <- length(thisgroup)
					ext.probs.agg[i, ] <- ext.probs.agg[i, ] * ext.weight[i] 			
				}
			}
			ext.probs.agg.4fold <- rbind(ext.probs.agg, ext.probs.agg,
										 ext.probs.agg, ext.probs.agg)
			# multiply group size
			indiv <- indiv * weight

			indiv <- cbind(indiv[, 1:(external.causes[1] - 1)], 
				          ext.probs.agg.4fold, 
				          indiv[, external.causes[1]:C0])

			# divide by full size
			indiv <- indiv / (weight + rep(ext.weight, 4))

			colnames(indiv) <- colnames(object$indiv.prob)
		}else{
			id.out <- id
			ext.probs <- NULL
		}

		mean <- indiv[1:K, ]
		median <- indiv[(K+1):(2*K), ]
		lower <- indiv[(2*K+1):(3*K), ]
		upper <- indiv[(3*K+1):(4*K), ]


		
		# if no grouping, return a matrix
		if(is.null(by)){
			out <- data.frame(cbind(mean, median, lower, upper))
			colnames(out) <- c("Mean", "Median", "Lower", "Upper")
			return(out)
		}else{
			rownames(mean) <- rownames(median) <- rownames(lower) <- rownames(upper) <- allgroups
			return(list(mean = t(mean), median = t(median), 
					lower = t(lower), upper = t(upper)))
		}
		
	}


}




#' Update individual COD probabilities from InSilicoVA Model Fits
#' 
#' This function updates individual probabilities for each death and provide posterior credible intervals for each estimates. 
#' 
#' 
#' @param object Fitted \code{"insilico"} object.
#' @param CI Credible interval for posterior estimates.
#' @param java_option Option to initialize java JVM. Default to ``-Xmx1g'', which sets the maximum heap size to be 1GB.
#' @param \dots Not used.
#' @return object Updated \code{"insilico"} object.
#' @author Zehang Li, Tyler McCormick, Sam Clark
#' 
#' Maintainer: Zehang Li <lizehang@@uw.edu>
#' @seealso \code{\link{insilico}}, \code{\link{get.indiv}}
#' @references Tyler H. McCormick, Zehang R. Li, Clara Calvert, Amelia C.
#' Crampin, Kathleen Kahn and Samuel J. Clark Probabilistic cause-of-death
#' assignment using verbal autopsies, \emph{arXiv preprint arXiv:1411.3042}
#' \url{http://arxiv.org/abs/1411.3042} (2014)
#' @examples
#' \dontrun{
#' data(RandomVA1)
#' fit1a<- insilico(RandomVA1, subpop = NULL,  
#'                 Nsim = 1000, burnin = 500, thin = 10 , seed = 1,
#'                 auto.length = FALSE)
#' summary(fit1a, id = "d199")
#' 
#' # The following script updates credible interval for individual 
#' fit1b <- updateIndiv(fit1a, CI = 0.95)
#' summary(fit1b, id = "d199")
#' }
#' @export updateIndiv 

updateIndiv <- function(object, CI = 0.95, java_option = "-Xmx1g", ...){
		
	 indiv  <- get.indiv(object, CI = CI, ...)
	 object$indiv.prob.lower <- indiv$lower
	 object$indiv.prob.upper <- indiv$upper
	 object$indiv.prob.median <- indiv$median
	 object$indiv.CI <- CI

	 return(object)
}