#' Obtain conditional probabilities from training data
#'
#' This is the function internally used in \code{insilico.train} function.
#' 
#' @param train Training data, it should be in the same format as the testing data
#' and contains one additional column (see \code{cause} below) specifying known
#' cause of death. The first column is also assumed to be death ID. 
#' @param gs the name of the column in \code{train} that contains cause of death.
#' @param gstable The list of causes of death used in training data.
#' @param thre a numerical value between 0 to 1. It specifies the maximum rate of
#' missing for any symptoms to be considered in the model. Default value is set to
#' 0.95, meaning if a symptom has more than 95\% missing in the training data, it
#' will be removed.
#' @param type Three types of learning conditional probabilities are provided: ``quantile''
#' or ``fixed''. Since InSilicoVA works with ranked conditional probabilities P(S|C), ``quantile''
#' means the rankings of the P(S|C) are obtained by matching the same quantile distributions
#' in the default InterVA P(S|C), and ``fixed'' means P(S|C) are matched to the closest values
#' in the default InterVA P(S|C) table. Empirically both types of rankings produce similar results. The third option ``empirical'' means no rankings are calculated, only the raw P(S|C) values are returned.
#' @param isNumeric Indicator if the input is already in numeric form. If the
#' input is coded numerically such that 1 for ``present'', 0 for ``absent'',
#' and -1 for ``missing'', this indicator could be set to True to avoid
#' conversion to standard InterVA format.
#'
#' @return 
#'  \item{cond.prob}{raw P(S|C) matrix}
#' 	\item{cond.prob.alpha}{ranked P(S|C) matrix}
#' 	\item{table.alpha}{list of ranks used}
#' 	\item{table.num}{list of median numerical values for each rank}
#' 	\item{symps.train}{training data after removing symptoms with too high missing rate.}
#' @export extract.prob
#'
extract.prob <- function(train, gs, gstable, thre = 0.95, type = c("quantile", "fixed", "empirical")[1], isNumeric = FALSE){

	## help functions to count vector contents
	if(isNumeric){
		yes <- 1
		no <- 0
		miss <- -1
	}else{
		yes <- "Y"
		no <- ""
		miss <- "."
	}

	countyes <- function(x){length(which(toupper(x) == yes))}
	countno <- function(x){length(which(x == no))}
	countmiss <- function(x){length(which(x == miss))}
	countna <- function(x){length(which(is.na(x)))}
	
	##
	## function to turn numeric matrix into level matrix
	##
	## 	cond.prob: S by C numeric matrix
	##	levels   : length L character vector (from small to large) 
	##				(e.g. E, D, C, B, A)
	##  cutoffs  : length L numeric vector (percentile cutoff, last = 1)
	##              (e.g. .1, .2, .3, .7, 1) -> E=(0, 0.1), D=(.1, .2), ..
	toLevel <- function(cond.prob, levels, cutoffs){
		cond.prob.vec <- as.vector(cond.prob)
		new <- rep(no , length(cond.prob.vec))
		cdf <- ecdf(cond.prob.vec)
		table <- table.median <- rep(0, length(cutoffs))
		for(i in length(cutoffs) : 1){
			table[i] <- quantile(cdf, cutoffs[i])	
			if(i > 1){
				table.median[i] <- quantile(cdf, (cutoffs[i] + cutoffs[i-1])/2)
			}else{
				table.median[i] <- quantile(cdf, (cutoffs[i])/2)
			}	
		}
		for(i in length(cutoffs) : 1){ 
			new[which(cond.prob.vec <= table[i])] <- levels[i]
		}
		new <- matrix(new, dim(cond.prob)[1], dim(cond.prob)[2])
		colnames(new) <- colnames(cond.prob)
		rownames(new) <- rownames(cond.prob)
		return(list(cond.prob = new, table.cut = table, table.median = table.median))
	}

	##
	##
	## step1: check missing items and update symptom matrix
    ##
    ##
	## first screening out symptoms with highest missing
	## remove missing rate > 0.95 
	train.id <- as.character(train[, 1])
	train <- train[, -1]
	train.gs <- as.character(train[, gs])
	if(is.null(gstable)){
		gstable <- unique(train.gs)
	}

	train <- train[, -which(colnames(train) == gs)]

	miss.train <- apply(train, 2, countmiss) / dim(train)[1]
	miss.too.many <- which(miss.train > thre)
	if(length(miss.too.many) > 0){
		train <- train[, -miss.too.many]
		warning(paste(length(miss.too.many), "symptoms deleted for missing rate over the pre-specified threshold", thre),
			immediate. = TRUE)
	}

	## check cause list
	gs.missing <- which(gstable %in% unique(train.gs) == FALSE)
	if(length(gs.missing) > 0){
		warning(paste("Causes not found in training data and deleted:",
			paste(gstable[gs.missing], collapse = ", ")), 
			immediate. = TRUE)
		gstable <- gstable[-gs.missing]
	}

	## 
	##
	## Get empirical prob base
	##
	##

	# get cond.prob, a S x C matrix of P(S|C)
	cond.prob <- NULL
	for(i in 1:length(gstable)){
		cause <- gstable[i]
		list <- which(train.gs == cause)
		
		cases <- train[list, ]
		## handle only case scenario
		if(length(list) == 1){
			cond <- rep(0, length(cases))
			cond[which(toupper(cases) == yes)] <- 1
			cond[which(cases == miss)] <-- NA
		}else{
			# remove missing from calculation
			count <- apply(cases, 2, countyes)
			cond <- count / (length(list) - apply(cases, 2, countmiss))						
		}
		cond.prob <- cbind(cond.prob, cond)
	}

	colnames(cond.prob)  <- gstable
	rownames(cond.prob)  <- colnames(train)

	S <- dim(cond.prob)[1]
	C <- dim(cond.prob)[2]

	## NA means for a particular cause and a particular symptom,
	## all cases are missing, 
	## but it's unreasonable to delete either, since they are useful
	## for other causes and symptoms
	bysymps <- apply(cond.prob, 1, countna)
	if(length(which(bysymps > 0)) > 0){
		cat(paste("\n", length(which(bysymps > 0)), "missing symptom-cause combinations in training data, imputed using average symptom prevalence\n"))
		## impute by Pr(s | c) by Prob(s | any c) 
		##   i.e. relatively no info provided from this s-c combination
		overall <- apply(train, 2, countyes) / dim(train)[1]
		for(i in 1:S){
			cond.prob[i, which(is.na(cond.prob[i,]))] <- overall[i]
		}		
	}

	## to avoid 0 and 1
	zeroes <- which(cond.prob == 0)
	ones <- which(cond.prob == 1)
	cond.prob[zeroes] <- runif(length(zeroes), min = 1e-6, max = 2e-6)
	cond.prob[ones] <- 1 - runif(length(ones), min = 1e-6, max = 2e-6)

	conditionals <- as.vector(cond.prob)

	# grab the InterVA table
	data("probbase", envir = environment())
	probbase <- get("probbase", envir  = environment())

	intertable <- probbase[2:246, 17:76]
	intertable[which(intertable == "B -")] <- "B-"
	intertable[which(intertable == "")] <- "N"
	levels <- c("I",  "A+", "A", "A-", "B+", "B" , "B-", "C+", "C", "C-", "D+", "D",  "D-", "E", "N")
	counts <- rep(0, length(levels))
	for(i in 1:length(levels)){
		counts[i] <- length(which(as.vector(intertable) == levels[i]))
	}
	freqs <- cumsum(counts / sum(counts))

	if(type == "quantile"){
		# get empirical results
		level.tmp <- toLevel(cond.prob, rev(levels), freqs)
		cond.prob.alpha <- level.tmp$cond.prob
		table.num <- rev(level.tmp$table.median)
	
	}else if(type == "fixed"){
		# get InterVA transformation results
		table.num <- c(1, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0001, 0.00001, 0)

		# get the middle point cut-off (first = 1)
		inter.cut <- c(1, table.num[-length(table.num)] + diff(table.num)/2 )

		# transform into the nearest InterVA values
		cond.prob.alpha <- matrix("N", dim(cond.prob)[1], dim(cond.prob)[2])
		for(i in 1:length(inter.cut)){
			cond.prob.alpha[which(cond.prob <= inter.cut[i])] <- levels[i]				
		}
		colnames(cond.prob.alpha) <- colnames(cond.prob) 
		rownames(cond.prob.alpha) <- rownames(cond.prob) 

	}else if(type == "empirical"){
		cond.prob.alpha <- NULL
		levels <- NULL
		table.num <- NULL
	}

	
	# the levels returned are from high to low
	
	return(list(cond.prob = cond.prob,
				cond.prob.alpha = cond.prob.alpha,
				table.alpha = levels, 
				table.num = table.num, 
				symps.train = train))
}
