## Created on Oct 1, 2015


#' Implement InSilicoVA methods with more flexible customization
#' 
#' This function implements InSilicoVA model. This is the lower level core function of InSilicoVA with more flexibility in customized input. For more detail of model specification, see the paper on \url{http://arxiv.org/abs/1411.3042} and the default function \code{\link{insilico}}.
#' 
#' 
#' @param data see \code{\link{insilico}}
#' @param isNumeric see \code{\link{insilico}}
#' @param updateCondProb see \code{\link{insilico}}
#' @param keepProbbase.level see \code{\link{insilico}}
#' @param CondProb see \code{\link{insilico}}
#' @param CondProbNum see \code{\link{insilico}}
#' @param datacheck see \code{\link{insilico}}
#' @param datacheck.missing see \code{\link{insilico}}
#' @param warning.write see \code{\link{insilico}}
#' @param external.sep see \code{\link{insilico}}
#' @param Nsim see \code{\link{insilico}}
#' @param thin see \code{\link{insilico}}
#' @param burnin see \code{\link{insilico}}
#' @param auto.length see \code{\link{insilico}}
#' @param conv.csmf see \code{\link{insilico}}
#' @param jump.scale see \code{\link{insilico}}
#' @param levels.prior see \code{\link{insilico}}
#' @param levels.strength see \code{\link{insilico}}
#' @param trunc.min see \code{\link{insilico}}
#' @param trunc.max see \code{\link{insilico}}
#' @param subpop see \code{\link{insilico}}
#' @param java_option see \code{\link{insilico}}
#' @param seed see \code{\link{insilico}}
#' @param phy.code see \code{\link{insilico}}
#' @param phy.cat see \code{\link{insilico}}
#' @param phy.unknown see \code{\link{insilico}}
#' @param phy.external see \code{\link{insilico}}
#' @param phy.debias see \code{\link{insilico}}
#' @param exclude.impossible.cause see \code{\link{insilico}}
#' @param customization.dev Logical indicator for customized variables
#' @param Probbase_by_symp.dev Not tested yet.
#' @param probbase.dev The customized conditional probabilities of symptoms given causes, which could be in a different format than InterVA default, but it should consist of \code{nlevel.dev} levels rather than numerical values.
#' @param table.dev The table of level names in \code{probbase.dev}. Default to be NULL
#' @param table.num.dev The corresponding prior numerical values for each level in \code{probbase.dev}, in the same order as \code{table.dev}. Default to be NULL
#' @param gstable.dev Table of gold standard causes for each death. Default to be NULL
#' @param nlevel.dev number of levels in \code{probbase.dev}. Default to be NULL
#' @param indiv.CI credible interval for individual probabilities
#' @param ... unused arguments

#' @return 
#' a insilico fit object, see see \code{\link{insilico}} for more detail.
#' @author Zehang Li, Tyler McCormick, Sam Clark
#' 
#' Maintainer: Zehang Li <lizehang@@uw.edu>
#' @seealso \code{\link{plot.insilico}}, \code{\link{summary.insilico}}
#' @references Tyler H. McCormick, Zehang R. Li, Clara Calvert, Amelia C.
#' Crampin, Kathleen Kahn and Samuel J. Clark(2014) \emph{Probabilistic
#' cause-of-death assignment using verbal autopsies},
#' \url{http://arxiv.org/abs/1411.3042} \cr \emph{Working paper no. 147, Center
#' for Statistics and the Social Sciences, University of Washington}
#' @keywords InSilicoVA
#' 
#' @export insilico.fit
insilico.fit <- function(data, isNumeric = FALSE, updateCondProb = TRUE, keepProbbase.level = TRUE,  CondProb = NULL, CondProbNum = NULL, datacheck = TRUE, datacheck.missing = TRUE, warning.write = FALSE, external.sep = TRUE, Nsim = 4000, thin = 10, burnin = 2000, auto.length = TRUE, conv.csmf = 0.02, jump.scale = 0.1, levels.prior = NULL, levels.strength = 1, trunc.min = 0.0001, trunc.max = 0.9999, subpop = NULL, java_option = "-Xmx1g", seed = 1, phy.code = NULL, phy.cat = NULL, phy.unknown = NULL, phy.external = NULL, phy.debias = NULL, exclude.impossible.cause = TRUE, customization.dev = FALSE, Probbase_by_symp.dev = FALSE, probbase.dev = NULL, table.dev = NULL, table.num.dev = NULL, gstable.dev = NULL, nlevel.dev = NULL, indiv.CI = NULL, ...){ 
  # handling changes throughout time
  args <- as.list(match.call())
  if(!is.null(args$length.sim)){
  	Nsim <- args$length.sim
  	cat("length.sim argument is replaced with Nsim argument, will remove in later versions.\n")
  }

##----------------------------------------------------------##
##       Helper functions                                   ##
##----------------------------------------------------------##

InterVA.table <- function(standard = TRUE, min = NULL, table.num.dev = NULL){
###########################################################
# function to return the interVA conversion table for alphabetic levels
# also change the smallest value from 0 to user input		
# @param:
#       standard: if TRUE, only need min and use interVA standard values
#		min: minimum level to replace 0 in interVA standard
#       table.num.dev: customized values
# @values:
# 		vector of interVA4 levels
	if(standard){
		if(is.null(min)){stop("Minimum level not specified")}
		return(c(1, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 
			  0.001, 0.0005, 0.0001, 0.00001, min))		
	}else{
		if(is.null(table.num.dev)){
			stop("Numerical level table not specified")
		}
		if(min(table.num.dev) == 0){
			table.num.dev[which.min(table.num.dev)] <- sort(table.num.dev, decreasing=FALSE)[2]/10			
		}
		return(sort(table.num.dev, decreasing = TRUE))
	}
}

scale.vec <- function(aaa, scale = NULL, scale.max = NULL, reverse = TRUE){
###########################################################
# function to change the scale of a given vector
# Note: need to reverse the order! so A has larger value, F has smaller value
# @param:
#	    aaa       : vector to be scaled
# 		scale     : sum of the vector after scaling
#       scale.max : max of the vector after scaling
#		reverse   : whether the order should be reversed
# @values:	
#  		scaled vector
	if(reverse){aaa <- max(aaa) + 1 - aaa}
	bbb <- aaa
	if(!is.null(scale)) return(bbb/sum(bbb) * scale)
	if(!is.null(scale.max)) return(bbb * scale.max / max(bbb))
}

scale.vec.inter <- function(aaa, scale = NULL, scale.max = NULL){
###########################################################
# function to map ordered input vector to InterVA alphabetic scale
# Note:  to avoid 0, change the last one to 0.000001 here
# @param:
#	    aaa       : vector to be scaled
# 		scale     : sum of the vector after scaling
#       scale.max : max of the vector after scaling
#		reverse   : whether the order should be reversed
# @values:	
#  
	dist <- InterVA.table(standard = TRUE, min = 0.000001)
	if(length(aaa) != length(dist)){stop("dimension of probbase prior not correct")}
	bbb <- dist[order(aaa)]
	if(!is.null(scale)) return(bbb/sum(bbb) * scale)
	if(!is.null(scale.max)) return(bbb * scale.max / max(bbb))
}

change.inter <- function(x, order = FALSE, standard = TRUE, table.dev = NULL, table.num.dev = NULL){
###########################################################
# function to translate alphebatic matrix into numeric matrix or order matrix
# @param:
# 	x      : alphabetic matrix 
#	order  : whether to change the matrix into order matrix
#   standard: whether to use the standard table
#   table.dev: new table of level names used high to low
#   table.num.dev: new table of numerical values correspond to table.dev
# @values:
#	numeric matrix by InterVA probbase rules, or the order matrix
	a <- dim(x)[1]
	b <- dim(x)[2]
	if(is.null(a)){
		y <- rep(0,length(x))
	}else{
		y <- matrix(0, a, b)
	}  	
	inter.table <- InterVA.table(standard = standard, table.num.dev = table.num.dev, min = 0)
	if(is.null(table.dev)){
		y[x == "I"] <- inter.table[1]
	    y[x == "A+"] <- inter.table[2]
	    y[x == "A"] <- inter.table[3]
	    y[x == "A-"] <- inter.table[4]
	    y[x == "B+"] <- inter.table[5]
	    y[x == "B"] <- inter.table[6]
	    y[x == "B-"] <- inter.table[7]
	    # historical error in probbase
	    y[x == "B -"] <- inter.table[7]

	    y[x == "C+"] <- inter.table[8]
	    y[x == "C"] <- inter.table[9]
	    y[x == "C-"] <- inter.table[10]
	    y[x == "D+"] <- inter.table[11]
	    y[x == "D"] <- inter.table[12]
	    y[x == "D-"] <- inter.table[13]
	    y[x == "E"] <- inter.table[14]
	    y[x == "N"] <- inter.table[15]
	    y[x == ""] <- inter.table[15]
	}else{
		if(length(table.dev) != length(table.num.dev)){
			stop("table.dev and table.num.dev have different length")
		}
		for(i in 1:length(table.dev)){
			y[x == table.dev[i]] <- inter.table[i] 
		}
	}


    if(order){
    	for(i in 1: length(inter.table)){
			y[y == inter.table[i] ] <- i
		}
    }
    if(!is.null(a)){
		y <- matrix(y, a, b)
	}  	
    return(y)   
}

cond.initiate <- function(probbase.order, expIni, Inter.ini, min, max){###########################################################
# Randomly initialize probbase from order matrix
# @param:
# 	probbase.order : order matrix
#   expIni         : initialize to exponential of uniform
#   Inter.ini      : initialize to InterVA probbase values
#   min            : minimum of probbase values
#	max            : maximum of probbase values
# @values:
#	new probbase matrix

		# move initial values away from 0 and 1
		if(min == 0) min <- 0.01
		if(max == 1) max <- 0.99		
		# find the number of levels
		nlevel <- max(probbase.order)
		# if the random levels are in log-linear fashion
		if(expIni == TRUE){
		randomlevels <- sort(exp(runif(nlevel, log(min), log(max))), decreasing = TRUE)	
		}else{
		randomlevels <- sort(runif(nlevel, min, max), decreasing = TRUE)		  
		}
		# if the random levels are initialized proportional to InterVA intervals
		if(Inter.ini){
			randomlevels <- InterVA.table(standard = TRUE, min = 0)
			randomlevels <- (randomlevels * (max - min)) + min
		}
		# change order matrix into actual values
		probrandom = probbase.order
		for(i in 1:nlevel){
			probrandom[which(probbase.order == i)] <- randomlevels[i]
		}

		return(probrandom)	
}

		
removeBad <- function(data, is.numeric, subpop){
###########################################################
# function to remove data with no sex/age indicators 
# @param:
#		data        : as in main function
#		is.numeric  : as in main function
#		subpop      : as in main function
# @values:
# 		a list of data and subpop after removing bad data
	if(!is.numeric){
		data.num <- matrix(0, dim(data)[1], dim(data)[2])		
		for(i in 2:dim(data)[2]){
			temp <- data[,i]
			temp.ind <- which(temp == "Y")
			data.num[temp.ind, i] <- 1
		}
	}else{
		data.num <- data
	}
	err <- NULL
	for(i in 1:dim(data.num)[1]){
		if(sum(data.num[i, 2:8]) < 1 ||
			sum(data.num[i, 9:10]) < 1 ||
			sum(data.num[i, 23:223]) < 1){
			err <- c(err, i)
		}
	}
	if(is.null(err)) return(list(data = data, subpop = subpop))
	if(!is.null(subpop)) return(list(data = data[-err, ], 
									 subpop = subpop[-err]))
	if(is.null(subpop)) return(list(data = data[-err, ], subpop = NULL))
}

datacheck.interVAJava <- function(data, obj, warning.write){
		
		data("probbase", envir = environment())
		probbase <- get("probbase", envir  = environment())

		# get text matrix
		dontask0 <- probbase[-1, 4:11]
		askif0 <- probbase[-1, 12]

		# get numerical matrix
		symps <- colnames(data[, -1])
		dontask <- match(as.vector(dontask0), symps)
		dontask[is.na(dontask)] <- 0
		dontask <- matrix(as.integer(dontask), dim(dontask0)[1], dim(dontask0)[2])
		askif <- as.integer(match(askif0, symps))
		askif[is.na(askif)] <- as.integer(0)
		askif <- matrix(askif, ncol = 1)

		data.num <- matrix(0, dim(data)[1], dim(data)[2] - 1)
		for(j in 2:dim(data)[2]){
			data.num[which(data[, j] == "Y"), j - 1] <- 1			
			data.num[which(data[, j] == "."), j - 1] <- -1
		}

		dontask.j <- .jarray(dontask, dispatch = TRUE)
		askif.j <- .jarray(askif, dispatch = TRUE)
		data.j <- .jarray(data.num, dispatch = TRUE)

		if(!warning.write){
			checked  <- .jcall(obj, "[[D", "Datacheck", dontask.j, askif.j, data.j)
		}else{
			id.j <- .jarray(as.character(data[, 1]), dispatch = TRUE)
			symps.j <- .jarray(colnames(data)[-1], dispatch = TRUE)
			checked  <- .jcall(obj, "[[D", "Datacheck", dontask.j, askif.j, data.j, id.j, symps.j, "warning_insilico.txt")
		}

		return(do.call(rbind, lapply(checked, .jevalArray)))
}


removeExt <- function(data, prob.orig, is.Numeric, subpop, subpop_order_list, external.causes, external.symps){
###########################################################
# function to remove external causes/symps and assign deterministic deaths
# @param:
#	   data
#      prob.orig: directly from InterVA4
# @values:
#	   data: after removing external death
#	   prob.orig: after deleting external symptoms
#	   exts: external death list		
	extSymps <- external.symps
	extCauses <- external.causes
	# extract subset of data for external symptoms
	N.all <- dim(data)[1]
	extData <- data[, extSymps + 1]
	if(is.Numeric){
		neg <- 0
		pos <- 1
	}else{
		neg <- ""
		pos <- "Y"
	}
	ext.where <- which(apply(extData, 1, function(x){
									length(which(x == pos)) }) > 0)
	extData <- as.matrix(extData[ext.where, ])
	ext.id <- data[ext.where, 1]
	ext.sub <- subpop[ext.where]
	# a smaller scale datacheck for external causes
	extData[which(extData[,3] == pos), 2] <- neg
	extData[which(extData[,4] == pos), c(2,3,5,6)] <- neg
	extData[which(extData[,5] == pos), c(2,3,4)] <- neg
	extData[which(extData[,6] == pos), c(8,10)] <- neg
	extData[which(extData[,8] == pos), c(2,3,4,5)] <- neg
	extData[which(extData[,9] == pos), c(2,3,4,8,10)] <- neg
	extData[which(extData[,10] == pos),c(2,4,5,6,8)] <- neg
	extData[which(extData[,11] == pos),c(3,4,8,9)] <- neg
	extData[which(extData[,12] == pos),c(4,7,8)] <- neg
	
	# initialize with all "unspecified ext causes"
	ext.cod <- rep(extCauses[11], length(ext.id))
	# begin checking symptoms
	# road traffic
	ext.cod[which(extData[,2] == pos)] <- extCauses[1]
	# non-road transport
	ext.cod[which(extData[,3] == pos)] <- extCauses[2]
	# accident fall
	ext.cod[which(extData[,4] == pos)] <- extCauses[3]
	# drowning
	ext.cod[which(extData[,5] == pos)] <- extCauses[4]
	# burns
	ext.cod[which(extData[,6] == pos)] <- extCauses[5]
	# assault
	ext.cod[which(extData[,7] == pos)] <- extCauses[10]
	# animal bite
	ext.cod[which(extData[,8] == pos)] <- extCauses[6]
	# force of nature
	ext.cod[which(extData[,9] == pos)] <- extCauses[7]
	# poison
	ext.cod[which(extData[,10] == pos)] <- extCauses[8]
	# inflict
	ext.cod[which(extData[,11] == pos)] <- extCauses[10]
	# suicide
	ext.cod[which(extData[,12] == pos)] <- extCauses[9]

	# delete death confirmed external
	data <- data[-ext.where, ]

	data <- data[, -(extSymps + 1)]
	# delete the causes from probbase
	prob.orig <- prob.orig[ -(extSymps), -(extCauses)]
	if(!is.null(subpop)){
		ext.csmf <- vector("list", length(subpop_order_list))
		for(i in 1:length(ext.csmf)){
			ext.csmf[[i]] <- rep(0, length(extCauses))
			ext.cod.temp <- ext.cod[which(ext.sub == subpop_order_list[i])]
			if(!is.null(ext.cod.temp)){
				for(j in 1:length(extCauses)){
					ext.csmf[[i]][j] <- length(which(ext.cod.temp == extCauses[j]))
				}
				ext.csmf[[i]] <- ext.csmf[[i]]/length(which(subpop == subpop_order_list[i]))		
			}
		}
	}else{
		ext.csmf <- rep(0, length(extCauses))
		for(i in 1:length(extCauses)){
			ext.csmf[i] <- length(which(ext.cod == extCauses[i]))
		}
		ext.csmf <- ext.csmf/N.all		
	}
	return(list(data = data, 
				subpop = subpop[-ext.where],
				prob.orig = prob.orig, 
				ext.sub  = ext.sub,
				ext.id = ext.id, 
				ext.cod = ext.cod,
				ext.csmf = ext.csmf))
}

ParseResult <- function(N_sub.j, C.j, S.j, N_level.j, pool.j, fit){
###########################################################
# function to parse results from Java into correct place
# @param:
#	   various java arguments
# 	   Java output
# @values:
#		list of variables parsed   
   	counter <- 1
    csmf.sub <- NULL
    p.hat <- NULL
    probbase.gibbs <- NULL
    levels.gibbs <- NULL

    # extract N_thin
    N_thin = fit[1]
    counter = counter + 1

    # extract CSMF
    if(N_sub.j > 1){
        csmf.sub <- vector("list", N_sub.j)
        for(sub in 1:N_sub.j){
            csmf.sub[[sub]] <- fit[counter : (counter + C.j * N_thin - 1)]
            csmf.sub[[sub]] <- matrix(csmf.sub[[sub]], nrow = N_thin, ncol = C.j, 
            						  byrow = TRUE)
            counter <- counter + C.j * N_thin
        }
    }else{
        p.hat <- fit[counter : (counter + C.j * N_thin- 1)]
        p.hat <- matrix(p.hat, nrow = N_thin, ncol = C.j, byrow = TRUE)
        counter <- counter + C.j * N_thin
    }
    
    # extract individual probabilities
    p.indiv <- fit[counter : (counter + N.j * C.j - 1)]
    p.indiv <- matrix(p.indiv, nrow = N.j, ncol = C.j, byrow = TRUE)
    counter <- counter + N.j * C.j

    # extract probbase 
    if(pool.j != 0){
        probbase.gibbs <- fit[counter:(counter + S.j * C.j * N_thin - 1)]
        # array(..., dim =c(a,b,c))
        # what it does it for each c, fill by column
        # i.e. c(x[1,1,1], x[2, 1, 1], x[1, 2, 1], ...) 
        # i.e. in Java, loop in the order of C.j -> S.j -> N_thin
        probbase.gibbs <- array(probbase.gibbs, dim = c(N_thin, S.j, C.j))
        counter <- counter + S.j * C.j * N_thin
    }else{
        levels.gibbs <- fit[counter:(counter + N_level.j * N_thin - 1)]
        levels.gibbs <- matrix(levels.gibbs, nrow = N_thin, ncol = N_level.j, 					   byrow = TRUE)
        counter <- counter + N_level.j * N_thin
    }

    # find last time configurations
    mu_last <- matrix(fit[counter:(counter + N_sub.j * C.j - 1)], 
                     nrow = N_sub.j, ncol = C.j, byrow = TRUE)
    counter <- counter + N_sub.j * C.j
    sigma2_last <-fit[counter:(counter + N_sub.j - 1)]
    counter <- counter + N_sub.j  
    theta_last <-  matrix(fit[counter:(counter + N_sub.j * C.j - 1)], 
                         nrow = N_sub.j, ncol = C.j, byrow = TRUE)
    counter <- counter + N_sub.j * C.j     	
    

    out <- list(csmf.sub = csmf.sub, p.hat = p.hat, p.indiv = p.indiv, 
                probbase.gibbs = probbase.gibbs, 
                levels.gibbs = levels.gibbs, 
                mu.last = mu_last, sigma2.last = sigma2_last, 
                theta.last = theta_last)

    return(out)
}

##----------------------------------------------------------##
##       Helper functions all loaded                        ##
##----------------------------------------------------------##

	if(is.null(java_option)) java_option = "-Xmx1g"
	options( java.parameters = java_option )
	obj <- .jnew("sampler/InsilicoSampler2")

	time0 <- Sys.time()
	# method <- tolower(method)
	# restrict to only normal model for now
	method <- "normal"
	alpha.scale = NULL

	##----------------------------------------------------------##
	## check model-specific parameters are provided
	if(method == "dirichlet"){
		if(is.null(alpha.scale)){stop("No alpha provided for Dirichlet prior")}
	}else if(method == "normal"){
		if(is.null(jump.scale)){stop("No jump scale provided for Normal prior")}
	}else{
		stop("invalid method")
	}

	##----------------------------------------------------------##
	if(is.null(Nsim) || is.null(thin) || is.null(burnin)){
		stop("Length of chain/thinning/burn-in not specified")
	}
	
	##----------------------------------------------------------##
	if(keepProbbase.level && Probbase_by_symp.dev){
		stop("keepProbbase.level and Probbase_by_symp.dev cannot be set to TRUE simultaneously.")
	}

	##----------------------------------------------------------##
	if(!is.null(nlevel.dev)){
		nlevel <- nlevel.dev
	}else{
		nlevel <- 15
	}

	##----------------------------------------------------------##
	## initialize key data dependencies
	##----------------------------------------------------------##	
	data("probbase", envir = environment())
	probbase<- get("probbase", envir  = environment())
	data("causetext", envir = environment())
	causetext<- get("causetext", envir  = environment())
	
	
	##----------------------------------------------------------##
	## Extract sub-populations		
	##----------------------------------------------------------##
	# get subpopulation if it's a columnname
  	if(class(subpop) == "list" || length(subpop) == 1){
  		col.index <- match(subpop, colnames(data))
  		if(length(which(is.na(col.index))) > 0){
  			stop("error: invalid sub-population name specification")
  		}
  		if(length(col.index) == 1){
  			subpop <- data[, col.index]
  		}else{
  			subpop <- data[, col.index]
  			subpop <- apply(subpop, 1, function(x){paste(x, collapse = ' ')})
  		}
  	}
  	# the subpop either from specification above, or vector input
  	# make sure they are not factors but characters
  	if(!is.null(subpop)){
  		subpop <- as.character(subpop)
  	} 
  	if(length(unique(subpop)) == 1){
  			subpop <- NULL
  			warning("Only one level in subpopulation, running the dataset as one population")
  	}

	##----------------------------------------------------------##
	## without developer customization		
	##----------------------------------------------------------##
	if(!customization.dev){
		# get interVA probbase
	  	prob.orig <- probbase[2:246,17:76]
	  	
	  	
	  	if(dim(data)[2] != dim(probbase)[1] ){
	  		correct_names <- probbase[2:246, 2]
	  		exist <- correct_names %in% colnames(data)
	  		if(length(which(exist == FALSE)) > 0){
		        stop(paste("error: invalid data input format. Symptom(s) not found:", correct_names[!exist]))
	  		}else{
	  			data <- data[, c(1, match(correct_names, colnames(data)))]
	  		}
	    }

	    ## check the column names and give warning
	    data("RandomVA1", envir = environment())
	    RandomVA1 <- get("RandomVA1", envir  = environment())
	    valabels <- colnames(RandomVA1)
	    vacauses <- causetext[4:63,2]
	    external.causes = seq(41, 51)
	    external.symps = seq(211, 222)
	    
	    count.changelabel = 0
	    for(i in 1:246){
	        if(tolower(colnames(data)[i]) != tolower(valabels)[i]){
	            warning(paste("Input columne '", colnames(data)[i], "' does not match InterVA standard: '", 
	                    valabels[i], "'", sep = ""),
	                    call. = FALSE, immediate. = TRUE)
	            count.changelabel = count.changelabel + 1
	        }         
	    }
	    if(count.changelabel > 0){
	        warning(paste(count.changelabel, "column names changed in input. \n If the change in undesirable, please change in the input to match standard InterVA4 input format.\n"), call. = FALSE, immediate. = TRUE)
	        colnames(data) <- valabels
	    }

	  	if(!is.null(CondProb)){
	  		prob.orig <- CondProb
	  		exclude.impossible.cause <- FALSE
	  	}
	  	if(!is.null(CondProbNum)){
	  		prob.orig <- CondProbNum 
	  		updateCondProb <- FALSE			
	  	}
	 
		##-----------------------------------------------------##
		## remove bad data happens before taking into missing
		## (i.e. data without age/sex or has no real symptoms)
		tmp <- removeBad(data, isNumeric, subpop)
	  	data <- tmp[[1]]
	  	subpop <- tmp[[2]]
  	  	if(is.null(subpop)){
	  		subpop_order_list <- NULL
	  	}else{
		  	subpop_order_list <- sort(unique(subpop))
	  	}

  	}else{
		##----------------------------------------------------------##
		## with developer customization		
		##----------------------------------------------------------##	
		# only match columns exactly as in probbase
	  	correct_names <- rownames(probbase.dev)
	  	exist <- correct_names %in% colnames(data)
	  	if(sum(exist) == 0){
		    stop("error: invalid data input, no matching symptoms found")
	  	}else{
	  		data <- data[, c(1, match(correct_names, colnames(data)))]
	  	}
	    
	    prob.orig <- probbase.dev	
  		valabels <- colnames(data)
	    vacauses <- gstable.dev
	    external.causes <- NULL
  	}

  	# standardize to Upper case
	data <- data.frame(lapply(data, as.character), 
					   stringsAsFactors=FALSE)
	for(j in 2:dim(data)[2]){
		data[is.na(data[, j]), j] <- ""
		data[, j] <- toupper(data[, j])
	}
  	##----------------------------------------------------------##
  	## remove external causes
  	if(external.sep){
  		externals <- removeExt(data,prob.orig, isNumeric, subpop, subpop_order_list, external.causes, external.symps)
  		data <- externals$data
  		subpop <- externals$subpop
  		prob.orig <- externals$prob.orig
  	}

  	if(datacheck.missing && datacheck){
		cat("Performing data consistency check...\n")
		checked <- datacheck.interVAJava(data, obj, warning.write)
		cat("Data check finished.\n")
		## update in data with missing, nothing is updated into missing, so only changing Y and N
		for(i in 1:(dim(data)[2]-1)){
			data[which(checked[, i] == 1), i+1] <- "Y"
			data[which(checked[, i] == -1), i+1] <- ""
		}	
  	}
	##----------------------------------------------------------##
   	## check the missing list
   	## this step is after removing bad data and before data-checking
   	missing.all <- NULL
	## add the all missing items that are not in missing list
	for(i in 1:(dim(data)[2]-1)){
		if(length(which(data[,i+1] == ".")) >= 1 * dim(data)[1] && 
		   !((i+1) %in% missing.all)){
			missing.all <- sort(c(missing.all, i))
		}
	}
	if(length(missing.all) > 0){
		warning(paste(length(missing.all), "symptom missing completely and added to missing list", 
			"\nList of missing symptoms: \n", 
			paste( probbase[missing.all + 1, 2], collapse = ", ")), 
		call. = FALSE, immediate. = TRUE)		
	}
	## remove all missing symptoms from both data and probbase
	if(!is.null(missing.all)){	
		data <- data[, -(missing.all + 1)]
		prob.orig <- prob.orig[-missing.all, ]	
	}	

	##----------------------------------------------------------##
   	## check interVA rules, ignoring missing at this step,
   	## since missing could be rewritten 
   	if((!datacheck.missing) && datacheck){
		cat("Performing data consistency check...\n")
		checked <- datacheck.interVAJava(data, obj, warning.write)
		cat("Data check finished.\n")
		## update in data with missing, nothing is updated into missing, so only changing Y and N
		for(i in 1:(dim(data)[2]-1)){
			data[which(checked[, i] == 1), i+1] <- "Y"
			data[which(checked[, i] == 0), i+1] <- ""
		}	
	}


	##----------------------------------------------------------##	
	# initiate numerical matrix "cond.prob.true"
	# obtained from "prob.orig", which is could be one of the following: 
	#	1. level matrix: if CondProbNum is NULL
	#   2. already a numerical matrix: if otherwise

	##----------------------------------------------------------##
	## without developer customization		
	##----------------------------------------------------------##
	if(!customization.dev){
	  	## translate original probbase into InterVA interpreted values
	  	if(is.null(CondProbNum)){
	  		prob.order <- change.inter(prob.orig, order = TRUE, standard = TRUE)
	  		cond.prob.true <- change.inter(prob.orig, order = FALSE, standard = TRUE)
	 	}else{
	 		cond.prob.true <- prob.orig
	 		prob.order <- matrix(1, dim(prob.orig)[1], dim(prob.orig)[2])
	 	}
	}else{
		##----------------------------------------------------------##
		## with developer customization		
		##----------------------------------------------------------##
	 	# TODO: make sure this is ok 
	 	# if(!is.null(prob.order.dev)){
	  	# 		prob.order <- prob.order.dev
	  	# 		cond.prob.true <- prob.orig
	  	# 	}else{
	  	if(updateCondProb){
	 		prob.order <- change.inter(prob.orig, order = TRUE, standard = FALSE, table.dev = table.dev, table.num.dev = table.num.dev)
	  		cond.prob.true <- change.inter(prob.orig, order = FALSE, standard = FALSE, table.dev = table.dev, table.num.dev = table.num.dev)
	  	}else{
	  		cond.prob.true <- prob.orig
	 		prob.order <- matrix(1, dim(prob.orig)[1], dim(prob.orig)[2])
	  	}
 	 } 	
	##----------------------------------------------------------##
  	# get data dimensions
  	## number of data
  	N <- dim(data)[1]
	if(N <= 1) stop("Not enough sample")
	## number of symptoms (need to delete ID column)
	S <- dim(data)[2] - 1
	if(S != dim(cond.prob.true)[1]) stop("Length of symptoms is not right")
	## number of causes
	C <- dim(cond.prob.true)[2]

	## external causes have been removed 
	if(external.sep){
		vacauses.current <- vacauses[-external.causes]
	}else{
		vacauses.current <- vacauses
	}

	##----------------------------------------------------------##
	## check impossible pairs of symptoms and causes
  	## check only first demographic symptoms (7 age + 2 gender)
  	## also the value saved is the index (starting from 1)
  	if(exclude.impossible.cause && (!customization.dev)){
	  	impossible <- NULL
	  	demog.set <- c("elder", "midage", "adult", "child", "under5", "infant", "neonate", "male", "female")
	  	demog.index <- match(demog.set, colnames(data)[-1])
	  	demog.index <- demog.index[!is.na(demog.index)]
  		for(ss in demog.index){
			for(cc in 1:C){
				if(cond.prob.true[ss, cc] == 0){
					impossible <- rbind(impossible, c(as.integer(cc), as.integer(ss)))
				}
			}
		}
		impossible <- as.matrix(impossible)	
	}else{
		# java checks if impossible has 2 columns
		impossible <- matrix(as.integer(0), 1, 3)
	}


	
	## physician debias data
	if(!is.null(phy.debias)){
	   phy.code <- phy.debias$code.debias
	}
	
	## physician codes
	if(!is.null(phy.code)){

		# also remove external causes too if needed
		if(external.sep){
			external.match <- match(phy.cat[, 1], vacauses[external.causes])
			phy.cat <- phy.cat[-which(!is.na(external.match)), ]
			phy.code <- phy.code[, -which(colnames(phy.code) == phy.external)]
		}
		# check every cause in phy.cat is an actual cause
		testmatch <- match(phy.cat[, 1], vacauses.current)
		if(length(which(is.na(testmatch))) > 0){
			stop("Cause name incorrect in physician category matching")
		}
		if(length(unique(phy.cat[,1])) != length(phy.cat[,1])){
			stop("Repeated cause names in physician category matching")
		}

		if(!(phy.unknown %in% colnames(phy.code))){
			stop("Cannot find Unknown category in phy.code")
		}		
		# causes in the physician coded categories
		cause.phy <- unique(phy.cat[, 2])
		if(phy.unknown %in% cause.phy){
			stop("Unknown cause exist in phy.cat matrix! Please remove it")
		}
		# Added varible for java:
		#   number of categories provided by physician + unknown
		C.phy <- length(cause.phy) + 1
		# when the two dimensions disagree
		if(dim(phy.code)[2] != C.phy + 1){
			# exist cause not appeared in phy.code, fine, just add empty column
			morecause <- which(!(cause.phy %in% colnames(phy.code)))
			if(length(morecause) > 0){
				empty <- matrix(0, dim(phy.code)[1], length(morecause))
				colnames(empty) <- cause.phy[morecause]
				phy.code <- cbind(phy.code, empty)
			}
			# exist cause not appeared in phy.cat, error
			if(dim(phy.code)[2] != C.phy + 1){
				stop("List of physician coded causes in phy.code and phy.cat do not match")
			}
		}
		# make unknown cause the first column
		matchunknown <- which(colnames(phy.code) == phy.unknown)
		restcolumns <- match(cause.phy, colnames(phy.code))
		phy.code <- phy.code[, c(1, matchunknown, restcolumns)]

		# Added varible for java:
		#   ordered set of broader cause
		vacauses.broader <- phy.cat[match(vacauses.current, phy.cat[,1]), 2]
		vacauses.broader <- match(vacauses.broader, cause.phy)
		# match ID
		matchid <- match(phy.code[, 1], data[, 1])
		# Added varible for java:
		#    assigned cause
		assignment <- matrix(0, N, C.phy)	
		# first column is unknown	
		assignment[, 1] <- 1
		# remove unmatched physician coding
		if(length(which(is.na(matchid))) > 0){
			phy.code <- phy.code[-which(is.na(matchid)), ]
			matchid <- matchid[-which(is.na(matchid))]
		}
		assignment[matchid, ] <- as.matrix(phy.code[, -1])
		
		#normalize assignment
		for(index in 1:dim(assignment)[1]){
			assignment[index, ] <- assignment[index, ] / sum(assignment[index, ])
		}
		
		if(external.sep){
			cat(paste(length(matchid), 
				"deaths found known physician coding after removing deaths from external causes.\n"))			
		}else{
			cat(paste(length(matchid), 
 				"deaths found known physician coding.\n"))
		}
	}else{
		# if no physician coding, everything is unknown
		C.phy <- 1
		assignment <- matrix(0, N, C.phy)	
		assignment[, 1] <- 1
		vacauses.broader <- 1:length(vacauses.current)
	}
	
	##---------------------------------------------------------------##
	## Specify the prior for truncated beta distribution
	prior.b.cond = trunc(1.5 * N)			
	#prior.b.cond = N

	if(keepProbbase.level){
		# if update only table of interpretation, need stronger prior
		levelcount <- table(cond.prob.true)
		# levels.strength should be multiplied here to avoid truncating to 0
		prior.b.cond <- trunc(prior.b.cond * median(levelcount) * levels.strength) 
	}else{
		prior.b.cond <- trunc(prior.b.cond * levels.strength)
	}

	if(is.null(levels.prior)){
			levels.prior <- scale.vec.inter(seq(1,nlevel), 
					scale.max = prior.b.cond * 0.99)		
	}
	##---------------------------------------------------------------##
	## get sub-population information
	if(!is.null(subpop)){
		if(length(subpop) != N) {
			stop("Sub-population size not match")
		}
		subpop.numeric <- rep(0, length(subpop))
		sublist <- vector("list", length(subpop_order_list))
		subbelong <- rep(0, N)
		for(i in 1:length(subpop_order_list)){
			sublist[[i]] <- which(subpop == subpop_order_list[i])
			subbelong[which(subpop == subpop_order_list[i])] <- i
			subpop.numeric[sublist[[i]]] <- i-1
		}
		names(sublist)<- subpop_order_list
		N.sub <- length(sublist)
	}
	##---------------------------------------------------------------##
  	## initiate indicator matrix	
	indic <- matrix(0, N, S)
	id <- data[, 1]
	## convert "Y" to 1 in indicator matrix
	## if numeric data, just use data metrix without ID
	if(isNumeric){
		indic <- data[, -1]
	}else{		
		for(i in 1:S){
			temp <- data[,i+1]
			temp.ind <- which(temp == "Y")
			indic[temp.ind, i] <- 1
		}
	}
	## if data contains missing ".", translate to -1
	if(isNumeric){
		contains.missing <- (length(which(data == "-1"))> 0)
	}else{
		contains.missing <- (length(which(data == ".")) > 0)
	}
	# if there are missing and data not numeric, initiate missing
	if(contains.missing && !(isNumeric)){
		for(i in 1:S){
			temp <- data[,i+1]
			temp.ind <- which(temp == ".")
			indic[temp.ind, i] <- -1
		}
	}
	## if data contains subpopulation, 
	## for complete missing within subpopulation,
	##	   change missing from -1 to -2 to separate.
	if(!is.null(subpop)){
		for(j in 1:N.sub){
			for(i in 1:S){
				temp <- indic[which(subbelong == j), i]
				if(length(which(temp == -1)) == length(temp)){
					indic[which(subbelong == j), i] <- -2	
				}
			}			
		}
	}

	##---------------------------------------------------------------##
	## parameter initialization
	# csmf.prior <- rep(1/C, C)
	Sys_Prior <- as.numeric(change.inter(probbase[1,17:76], order = FALSE), standard = TRUE)
	# Number of indicators + 13 description variables. A_group:14-16;B_group:17:76;D_group:77:81
	D <- length(Sys_Prior)
	csmf.prior <- Sys_Prior/sum(Sys_Prior)
	##---------------------------------------------------------------##
	##  initialize prior and adjust for external causes
	if(method == "dirichlet"){
		# determine alpha for CSMF
  		alpha <- csmf.prior * C * alpha.scale
  	}else if(method == "normal"){
  		if(!is.null(csmf.prior)){
	  		Z <- csmf.prior[1]/exp(1)
	  		mu <- log(csmf.prior / Z)
	  	}else{
	  		mu <- rep(1, C)
	  	}
	  	sigma2 <- 1
  	}
	if(external.sep){
		csmf.prior <- csmf.prior[-external.causes]
		if(method == "normal" && length(mu) != C){ 
			mu <- mu[-external.causes]
		}else if(method == "dirichlet" && length(alpha) != C){
			alpha <- alpha[-external.causes]
		}
	}
	##---------------------------------------------------------------##
	## initialize probbase and start java
	
	if(updateCondProb){
		cond.prob <- cond.initiate(prob.order, expIni = TRUE, Inter.ini = TRUE, min = trunc.min, max = trunc.max)
    }else if(is.null(CondProbNum)){
    	cond.prob <- change.inter(prob.orig, standard = TRUE)
    }else{
    	cond.prob <- prob.orig
    }



    N.j <- as.integer(N)
    S.j <- as.integer(S)
    C.j <- as.integer(C)
    probbase.j <- .jarray(as.matrix(cond.prob), dispatch=TRUE)
    
    # at this stage, there are a few possibilities:
    #  customization.dev & updateCondProb: need to check existing levels
    #  !updateCondProb: do nothing
    #  !is.null(CondProb): need to check existing levels
    #  
    if((customization.dev && updateCondProb) || !is.null(CondProb)){
    	
		# get new numerical levels
    	if(customization.dev){
	    	dist <- InterVA.table(standard = FALSE, table.num.dev = table.num.dev)
    	}else{
    		dist <- InterVA.table(standard = TRUE, min = 0)
    	}
    	# check existence of each level's index in prob.order
    	level.exist <- seq(1:length(dist)) %in% unique(as.vector(prob.order))
    	
    	# update order matrix
    	prob.order.new <- prob.order
    	for(i in 1:length(dist)){
    		if(!level.exist[i]){
    			# I is 1, N is 15
    			# if a level is missing, all level after should minus 1
    			prob.order.new[which(prob.order > i)] <- prob.order.new[which(prob.order > i)] - 1
    		}
    	}
    	prob.order <- prob.order.new
    	probbase_order.j <- .jarray(as.matrix(prob.order), dispatch = TRUE)
    	# update level vector
    	dist <- dist[level.exist]
    	levels.prior <- levels.prior[level.exist]
    	N_level.j <- as.integer(sum(level.exist))
    }else{
    	level.exist <- rep(TRUE, nlevel)
    	probbase_order.j <- .jarray(as.matrix(prob.order), dispatch = TRUE)
    	N_level.j <- as.integer(nlevel)
    	dist <- InterVA.table(standard = TRUE, min = 0)
	}
   
    level_values.j <- .jarray(dist, dispatch = TRUE)
    prior_a.j <- .jarray(levels.prior , dispatch = TRUE)
    prior_b.j <- prior.b.cond
    jumprange.j <- jump.scale
    trunc_min.j <- trunc.min
    trunc_max.j <- trunc.max
    indic.j <- .jarray(as.matrix(indic), dispatch=TRUE)
    contains_missing.j <- as.integer(contains.missing)
    # update Oct 1, 2015: pool variable redefined:
    # 					  0 - pool to table; 1 - by cause; 2 - by symptom
    pool.j <- as.integer(!keepProbbase.level) + as.integer(Probbase_by_symp.dev)
    seed.j <- as.integer(seed)
    N_gibbs.j <- as.integer(Nsim)
    burn.j <- as.integer(burnin)
    thin.j <- as.integer(thin)
    mu.j <- .jarray(mu, dispatch = TRUE)
    sigma2.j <- sigma2 
    isUnix <-  .Platform$OS.type == "unix"

    assignment.j <- .jarray(as.matrix(assignment), dispatch = TRUE)
    C.phy.j <- as.integer(C.phy)
    vacauses.broader.j <- .jarray(vacauses.broader+0.0, dispatch = TRUE)

    impossible.j <- .jarray(as.matrix(impossible), dispatch = TRUE)

    if(is.null(subpop)){
		N_sub.j <- as.integer(1)
		subpop.j <- .jarray(as.integer(rep(0, N)), dispatch = TRUE)
	}else{
	    N_sub.j <- as.integer(N.sub)
        subpop.j <- .jarray(as.integer(subpop.numeric), dispatch = TRUE)
	}

	dimensions <- c(N.j, S.j, C.j, N_sub.j, N_level.j)
	dimensions.j <- .jarray(dimensions, dispatch = TRUE)

    isAdded <- FALSE
    mu.last.j <- .jarray(matrix(0, N_sub.j, C), dispatch = TRUE)
    sigma2.last.j <- .jarray(rep(0, N_sub.j), dispatch = TRUE)
	theta.last.j <- .jarray(matrix(0, N_sub.j, C), dispatch = TRUE)
	keepProb.j <- !updateCondProb

	ins  <- .jcall(obj, "[D", "Fit", 
		dimensions.j, 
		probbase.j, probbase_order.j, level_values.j, 
		prior_a.j, prior_b.j, jumprange.j, trunc_min.j, trunc_max.j, 
		indic.j, subpop.j, contains_missing.j, pool.j, 
		seed.j, N_gibbs.j, burn.j, thin.j, 
		mu.j, sigma2.j, isUnix, keepProb.j, 
		isAdded, mu.last.j, sigma2.last.j, theta.last.j, 
		C.phy.j, vacauses.broader.j, assignment.j, impossible.j) 
    # one dimensional array is straightforward
    fit <- ins
    # fit <-  .jeval(ins, .jevalArray))
##---------------------------------------------------------------------------------##	
##
## save data from java output to proper format
##    
    results <- ParseResult(N_sub.j, C.j, S.j, N_level.j, pool.j, fit)
	# check convergence
    conv <- tryCatch({
					   if(!is.null(results$csmf.sub)){
							csmf.diag(results$csmf.sub, conv.csmf, 
								              test = "heidel", verbose = FALSE) 	
				    	}else{
				    		csmf.diag(results$p.hat, conv.csmf, 
								              test = "heidel", verbose = FALSE) 
				    	}
				}, error = function(condition) {
				    print("error checking diagnostics")
				    FALSE
				})
    # check convergence and if the chain needs to run longer
    if(auto.length){	
    	# if not converge, run again, max number of runs = 3
    	add <- 1
    	while(!conv && add < 3){
    		mu.last.j <- .jarray(as.matrix(results$mu.last), dispatch = TRUE)
    		sigma2.last.j <- .jarray(results$sigma2.last, dispatch = TRUE)
    		theta.last.j <- .jarray(as.matrix(results$theta.last), dispatch = TRUE)
    		# same length as previous chain if added the first time
    		# double the length	if the second time
    		Nsim <- Nsim * 2 
    		burnin <- Nsim / 2
    		N_gibbs.j <- as.integer(trunc(N_gibbs.j * (2^(add-1))))
			burn.j <- as.integer(0)
			keepProb.j <- !updateCondProb

			cat(paste("Not all causes with CSMF >", conv.csmf, "are convergent.\n"))
    		cat(paste("Increase chain length with another", N_gibbs.j, "iterations\n"))
    		obj <- .jnew("sampler/InsilicoSampler2")
    		ins  <- .jcall(obj, "[D", "Fit", 
						dimensions.j, 
						probbase.j, probbase_order.j, level_values.j, 
						prior_a.j, prior_b.j, jumprange.j, trunc_min.j, trunc_max.j, 
						indic.j, subpop.j, contains_missing.j, pool.j, 
						seed.j, N_gibbs.j, burn.j, thin.j, 
						mu.j, sigma2.j, isUnix, keepProb.j, 
						TRUE, mu.last.j, sigma2.last.j, theta.last.j, 
						C.phy.j, vacauses.broader.j, assignment.j, impossible.j)
    		# one dimensional array is straightforward
    		fit.add <- ins
    		# fit.add <-  t(sapply(ins, .jevalArray))

    		# remove the first row of fit)
    		results <- ParseResult(N_sub.j, C.j, S.j, N_level.j, pool.j, fit.add)
    		add = add + 1
    		# check convergence
	    	if(!is.null(results$csmf.sub)){
				conv <- csmf.diag(results$csmf.sub, conv.csmf, 
					              test = "heidel", verbose = FALSE) 	
	    	}else{
	    		conv <- csmf.diag(results$p.hat, conv.csmf, 
					              test = "heidel", verbose = FALSE) 
	    	}
    	}
    }
    ## if still not convergent 
    if(!conv){
    	cat(paste("Not all causes with CSMF >", conv.csmf, "are convergent.\n",
    			  "Please check using csmf.diag() for more information.\n"))
    }
    csmf.sub  <- results$csmf.sub 
    p.hat  <- results$p.hat
    p.indiv  <- results$p.indiv 
    probbase.gibbs  <- results$probbase.gibbs 
    levels.gibbs  <- results$levels.gibbs
##---------------------------------------------------------------------------------##
## To make output consistent for further analysis,
## 		add back external results
	##
	## if separated external causes
	##
    if(external.sep){
    	# get starting and ending index
    	ext1 <- external.causes[1]
    	ext2 <- external.causes[2]
    	##
    	## if with subgroup
    	##
    	if(!is.null(subpop)){
    		# set p.hat to NULL, and set up csmf.sub.all as a list of p.hat 
    		p.hat <- NULL
    		csmf.sub.all <- vector("list", N_sub.j)
    		names(csmf.sub.all) <- subpop_order_list
    		# iterate over all subpopulation
    		for(j in 1:length(csmf.sub)){
    			# initialize the csmf matrix 
    			csmf.sub.all[[j]] <- matrix(0, 
    				dim(csmf.sub[[j]])[1], 
    				C.j + length(external.causes))
    			# rescale the non-external CSMF once the external causes are added 
    			rescale <- length(sublist[[j]]) / (length(sublist[[j]]) + length(which(externals$ext.sub == subpop_order_list[j])))
    			temp <- csmf.sub[[j]] * rescale
    					
    			# combine the rescaled non-external CSMF with the external CSMF	
    			csmf.sub.all[[j]]  <- cbind(temp[, 1:(ext1 - 1)], 
				    						matrix(externals$ext.csmf[[j]], 
				    						  	   dim(temp)[1],
				    						  	   length(external.causes), 
				    						  	   byrow = TRUE), 
				    						temp[, ext1:C.j])	
    		}
    		csmf.sub <- csmf.sub.all
    	##
    	## if no subgroup
    	##
    	}else{
    		p.hat <- p.hat * N/(N + length(externals$ext.id))
    		temp <- p.hat[, ext1:C.j]
    		extra <- matrix(externals$ext.csmf, dim(p.hat)[1], length(external.causes), byrow = TRUE)
    		p.hat <- cbind(p.hat[, 1:(ext1 - 1)], extra, temp)
    	}

    	p.indiv <- cbind(p.indiv[, 1:(ext1 - 1)], 
    						  matrix(0, dim(p.indiv)[1], length(external.causes)), 
    						  p.indiv[, (ext1:C.j)])
    	
    	p.indiv.ext <- matrix(0, nrow = length(externals$ext.id), ncol = C.j + length(external.causes) )
    	for(i in 1:length(externals$ext.id)){p.indiv.ext[i, externals$ext.cod[i]] <- 1}
    	p.indiv <- rbind(p.indiv, p.indiv.ext) 
    	id <- c(id, externals$ext.id)
    	subpop <- c(subpop, externals$ext.sub)
   }
##---------------------------------------------------------------##    	
## add column names to outcome
if(pool.j != 0){
	##======= CHECK THIS CHUNK ===========##
	if(external.sep){
		# remove column "ID"
		valabels <- valabels[-1]
		# remove external
		valabels <- valabels[-external.symps]
		vacauses.ext <- vacauses[-external.causes]
	    # remove missing, sequential deleting since that's 
	    #   how we obtained the indices before
		valabels <- valabels[-missing.all]
		dimnames(probbase.gibbs)[[2]] <- valabels
		dimnames(probbase.gibbs)[[3]] <- vacauses.ext		
	}else{
		dimnames(probbase.gibbs)[[2]] <- valabels
		dimnames(probbase.gibbs)[[3]] <- vacauses
	}
	##======= CHECK THIS CHUNK ===========##
}else{
	if(customization.dev){
		colnames(levels.gibbs) <- rev(table.dev[level.exist])
	}else{
		colnames(levels.gibbs) <- c("I", "A+", "A", "A-", "B+", "B", "B-", "C+", "C", "C-", "D+", "D", "D-", "E", "N")
	}
	probbase.gibbs <- levels.gibbs
} 
if(!is.null(subpop)){
	  for(j in 1:length(csmf.sub)){colnames(csmf.sub[[j]]) <- vacauses}
	  p.hat <- csmf.sub
}else{
	  colnames(p.hat) <- vacauses
}
colnames(p.indiv) <- vacauses
rownames(p.indiv) <- id
##---------------------------------------------------------------------------------##    	
if(!updateCondProb){
    	probbase.gibbs <- NULL
}

cleandata <- as.matrix(indic)
# because the external deaths might be appended to the end
rownames(cleandata) <- id[1:dim(cleandata)[1]]

out <- list(
		id = id,
		data = cleandata,
	    indiv.prob = p.indiv, 
		csmf = p.hat,
		conditional.probs = probbase.gibbs,
		probbase = prob.orig,
		missing.symptoms = missing.all,
		external = external.sep, 
		external.causes = external.causes,
		impossible.causes = impossible,
	
		updateCondProb = updateCondProb, 
		keepProbbase.level = keepProbbase.level, 
		datacheck = datacheck,
		Nsim = Nsim, 
		thin = thin, 
		burnin = burnin, 
	
		jump.scale = jump.scale, 
		levels.prior = levels.prior, 
		levels.strength = levels.strength, 
		trunc.min = trunc.min, 
		trunc.max = trunc.max, 
		subpop = subpop, 
		indiv.CI = indiv.CI, 
		is.customized = customization.dev)

# get also individual probabilities
if(!is.null(indiv.CI)){
	indiv <- get.indiv(data = NULL, object = out, indiv.CI)
	out$indiv.prob.median <- indiv$median
	out$indiv.prob.upper <- indiv$upper
	out$indiv.prob.lower <- indiv$lower			
}

class(out) <- "insilico"
return(out)  	
} 
