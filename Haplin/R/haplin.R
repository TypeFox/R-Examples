haplin <- function(filename, data, pedIndex,
markers = "ALL", n.vars = 0, sep = " ", allele.sep = ";", na.strings = "NA",
design = "triad", use.missing = FALSE, xchrom = FALSE, maternal = FALSE, test.maternal = FALSE, poo = FALSE, scoretest = "no",
ccvar = NULL, strata = NULL, sex = NULL, comb.sex = "double",
reference = "reciprocal", response = "free", threshold = 0.01, max.haplos = NULL, haplo.file = NULL,
resampling = "no", max.EM.iter = 50, data.out = "no", verbose = TRUE, printout = TRUE
)
{
##
## filename: (character string) filename 
## n.vars is the number of columns before genetic data columns, sep = " " is the column separator, allele.sep = ";" separates alleles within marker, na.strings = "NA" identifies missing values
## ccvar is the position of the case-control variable, covar is the position of the environment covariate
## 
# PAABEGYNT 13/1-04 KL. 22.46
#
#
## EXTRACT AND CHECK ARGUMENTS
.info <- f.catch(match.call(), formals())
#
## SET PARAMETERS, FOR SIMPLICITY
design <- .info$model$design
xchrom <- .info$model$xchrom
maternal <- .info$model$maternal
response <- .info$haplos$response
max.EM.iter <- .info$control$max.EM.iter
#
# VET IKKE MED DISSE:
use.missing <- .info$model$use.missing
verbose <- .info$control$verbose
resampling <- .info$control$resampling
data.out <- .info$control$data.out
printout <- .info$control$printout
#
## INSTALL MASS (FOR THE mvrnorm FUNCTION):
#require(MASS)
#
## START 
if(verbose)cat("\n## HAPLIN, VERSION 5.5 ##\n")
#
## LOAD DATA, EITHER FROM FILE OR FROM DATA OBJECT
.data.read <- f.get.data(data, pedIndex, .info)
#
## DATA OUT, IF REQUESTED
if(data.out == "basic"){
	.tmp <- f.data(data.read = .data.read, quick = T)
	return(.tmp)
	#.data <- .tmp$data
	#.info <- .tmp$info
}
#
## PREPARE DATA, RETURN HERE
.tmp <- f.data(data.read = .data.read)
.data <- .tmp$data
.info <- .tmp$info
#
## DATA OUT, IF REQUESTED
if(data.out == "prelim"){
	.data.out <- f.prep.dataout(info = .info, data = .data)
	return(.data.out)
}
#
## DISSE BLIR NAA DEFINERT BARE I f.data OG MAA DERFOR GJORES EKSPLISITTE:
.ntri.seq <- .info$data$ntri.seq
ref.cat <- .info$haplos$ref.cat
#
## DECIDE TEST SEQUENCE (NEED ONLY SPECIFY ARGUMENTS THAT CHANGE OVER TEST SEQUENCE)
.testseq.standard <- list(list(maternal = F, response = "simple", x = F), list(maternal = maternal, response = response, x = T))
.messages.standard <- c("\nUsing EM to estimate model with no effect:\n", "\nUsing EM to estimate full model:\n")
#
.testseq <- .testseq.standard
.messages <- .messages.standard
#
if(!.info$model$test.maternal & .info$model$scoretest == "only"){
	## AS A STANDARD SEQUENCE, BUT WHERE THE LAST ESTIMATION IS DONE WITH ONLY ONE ITERATION, SIMPLY TO OBTAIN X (AND NUMBER OF PARAMETERS) IN FULL MODEL
	.testseq <- list(list(maternal = F, response = "simple", x = F), list(maternal = maternal, response = response, x = T, max.EM.iter = 0, suppressEmWarnings = T))
	###.res.oneiter <- f.EM.missing(data = .data, response = .info$haplos$response, maternal = maternal, verbose = verbose, ref.cat = ref.cat, design = design, max.EM.iter = 0, x = T, suppressEmWarnings = T)
	.messages <- c("\nUsing EM to estimate model with no effect:\n", "\nComputing design matrix for full model:\n")
		## MERK: EGENTLIG NOK BARE AA FAA UT X. *OG* ET ESTIMAT AV coefficients SIDEN LENGDEN AV DET FAKTISK BRUKES... (BURDE EGENTLIG IKKE VARE NODVENDIG)
		## DESSUTEN: BURDE HER KUNNE KJORE BAADE mult OG free. MEN: free GAAR IKKE NODVENDIGVIS FOR cc
}
#
if(.info$model$test.maternal){
	.testseq <- list(list(maternal = F, response = "simple", x = F), list(maternal = F, response = response, x = T), list(maternal = maternal, response = response, x = T))
	.messages <- c("\nUsing EM to estimate model with no effect:\n", "\nEstimating intermediate model, with only child effects:\n", "\nUsing EM to estimate full model:\n")
}
#
.test.response <- F
if(.test.response){
	.testseq <- list(list(maternal = F, response = "simple", x = F), list(maternal = maternal, response = "mult", x = T), list(maternal = maternal, response = "free", x = T))
	cat("LKJLKJL ADVARSEL MAA KJORES MED REF.CAT!?")
	.messages <- c("\nUsing EM to estimate model with no effect:\n", "\nEstimating intermediate model, with only child effects:\n", "\nUsing EM to estimate full model:\n")
}


##################################################################################################
#
## START ESTIMATION, INCLUDING EM (AND RESAMPLING, IF REQUESTED)
#
##################################################################################################
#
##
## ESTIMATE ALL MODELS IN .testseq
#
## "DEFAULT" ARGUMENTS TO EM
.EM.args <- list(data = .data, maternal = maternal, response = response, max.EM.iter = max.EM.iter, x = T, verbose = verbose, suppressEmWarnings = F, info = .info)
#
.res.list <- vector(length(.testseq), mode = "list")
for (i in seq(along = .res.list)){
	if(verbose)cat(.messages[i])
	#
	## REPLACE DEFAULT ARGMENTS WITH THOSE SPECIFIED IN .testseq[[i]]
	.EM.args.tmp <- .EM.args
	.EM.args.tmp[names(.testseq[[i]])] <- .testseq[[i]] 
	#
	###.res.list[[i]] <- f.EM.missing(data = .data, response = .testseq[[i]]$response, maternal = .testseq[[i]]$maternal, verbose = verbose, ref.cat = ref.cat, design = design, xchrom = xchrom, max.EM.iter = max.EM.iter, x = .testseq[[i]]$x)
	## DO THE ACTUAL ESTIMATION
	.res.list[[i]] <- do.call("f.EM.missing", .EM.args.tmp)
	if(verbose)cat("\nDone\n")
}
#
## EXTRACT CONVERGENCE DIAGNOSTICS
.info$estimation$iter.used <- sapply(.res.list, function(x) attr(x, "iter.used"))
.info$estimation$EM.conv <- sapply(.res.list, function(x) attr(x, "EM.conv"))


.res.0 <- .res.list[[1]]
.res <- .res.list[[length(.res.list)]]
#
## JUST CHECK THAT THE TRIAD ACCOUNTING IS CORRECT
if(scoretest == "only"){
	if(abs(.res.0$ntri - .ntri.seq[4]) > 1e-6) warning("There may be a problem with the data summary")
}else if(max(abs(.res$ntri - .ntri.seq[4]), abs(.res.0$ntri - .ntri.seq[4])) > 1e-6) warning("There may be a problem with the data summary")
#
#
if(resampling == "jackknife"){
	if(design != "triad") stop("Jackknifing has not been tested with case-control data....", call. = F)
	if(.info$model$scoretest == "only") stop('Jackknifing has not been tested when scoretest == "only"', call. = F)
	if(xchrom) stop("Jackknifing not tested with xchrom data!", call. = F)
	if(!is.null(.info$variables$covar))stop("Jackknifing not tested with covariates!", call. = F)
	if(verbose) cat("\nStarting jackknife\n")
	.res.resamp <- f.jackknife(data = .data, maternal = maternal, ref.cat = ref.cat, verbose = F, use.EM = T, max.EM.iter = max.EM.iter, info = .info)
	attr(.res, "cov.resamp") <- .res.resamp$cov
}
#
## DATA OUT, IF REQUESTED
if(data.out %in% c("null", "full")){
	if(data.out == "null") .restemp <- .res.0
	if(data.out == "full") .restemp <- .res
	.data.out <- f.prep.dataout(info = .info, data = .data, res = .restemp)
	return(.data.out)
}
#
## COMPUTE VARIANCE-COVARIANCE MATRICES, LIKELIHOOD RATIO TESTS AND SCORE TESTS
if(T){
	#
	## INITIALIZE LISTS FOR VAR.COVAR, SCORE AND LR
	.l <- length(.res.list)
	.o.var.covar.list <- vector(.l, mode = "list")## ER DET BRUK FOR DISSE???
	.score.list <- vector(.l, mode = "list")## BRUKES IKKE!!
	.scoretest.list <- vector(.l, mode = "list")
	.var.cov.list <- vector(.l, mode = "list")## ER DET BRUK FOR DISSE???
	.lratio.list <- vector(.l, mode = "list")
	.npars.0 <- rep(NA, .l)
	#
	if(.l > 2) for(i in 1:(.l-1)){
		#
		### PAIRWISE COMPARISON, STEP BY STEP
		### USED ONLY IF THERE ARE INTERMEDIATE STEPS IN .res.list, FOR INSTANCE WHEN TESTING OVERALL MATERNAL EFFECTS
		#
		## COMPUTE VAR.COVAR AND SCORE USING PREDICTION FROM i'th AND X FROM (i+1)'th
		.o.var.covar.list[[i]] <- f.var.covar(pred = .res.list[[i]]$pred, X = .res.list[[i+1]]$result$x, data = .data, info = .info)
		.npars.0[i] <- length(.res.list[[i]]$result$coefficients)
		.var.cov.list[[i]] <- .o.var.covar.list[[i]][["var.covar"]]
		.score.list[[i]] <- .o.var.covar.list[[i]][["score"]] ## BRUKES IKKE!!
		#
		## DO A SCORE TEST USING PREDICTION FROM i'th AND X FROM (i+1)'th
		.scoretest.list[[i]] <- f.scoretest(o.var.covar.0 = .o.var.covar.list[[i]], npars.0 = .npars.0[i])
		#
		## DO AN LR TEST, (i+1)'th AGAINST i'th
		.lratio.list[[i]] <- f.like.ratio(res.0 = .res.list[[i]], res = .res.list[[i+1]], data = .data, info = .info)
	}
	#
	### COMPARE FIRST WITH LAST, SHOULD ALWAYS BE DONE
	### SAVED AS LAST RESULT IN ALL LISTS
	#
	## COMPUTE VAR.COVAR AND SCORE USING PREDICTION FROM FIRST AND X FROM LAST
	.o.var.covar.list[[.l]] <- f.var.covar(pred = .res.list[[1]]$pred, X = .res.list[[.l]]$result$x, data = .data, info = .info)
	.npars.0[.l] <- length(.res.list[[1]]$result$coefficients) ## NB!
	.var.cov.list[[.l]] <- .o.var.covar.list[[.l]][["var.covar"]]
	.score.list[[.l]] <- .o.var.covar.list[[.l]][["score"]] ## BRUKES IKKE!!
	#
	## DO A SCORE TEST USING PREDICTION FROM FIRST AND X FROM LAST
	.scoretest.list[[.l]] <- f.scoretest(o.var.covar.0 = .o.var.covar.list[[.l]], npars.0 = .npars.0[.l])
	#
	## DO AN LR TEST, LAST AGAINS FIRST
	.lratio.list[[.l]] <- f.like.ratio(res.0 = .res.list[[1]], res = .res.list[[.l]], data = .data, info = .info)
}


# MERK: NAA ER DET INGEN OPSJON FOR test.maternal = T HVIS scoretest = F!



.tmp.0 <- .o.var.covar.list[[.l]]
.var.cov.0 <- .var.cov.list[[.l]] # NOTE THAT THIS IS COMPUTED FOR THE FULL MODEL, BUT ASSUMING H0. IT IS USED, FOR INSTANCE, IN f.suest


.score.ut <- .scoretest.list[[.l]]

#
## COMPUTE "EXACT" VAR-COVAR, TAKING EM INTO ACCOUNT
.temp1 <- f.var.covar(pred = .res$pred, X = .res$result$x, data = .data, info = .info)
.var.cov <- .temp1$var.covar
# .var.cov <- f.var.covar(pred = .res$pred, X = .res$result$x, data = .data, info = .info)[["var.covar"]]
attr(.res, "cov.correct") <- .var.cov


.lratio.ut <- f.like.ratio(res.0 = .res.0, res = .res, data = .data, info = .info)
.var.cov.ut <- list(var.cov.0 = .var.cov.0, var.cov = .var.cov)

if(verbose) cat("\nEstimation finished, preparing output...  ")
#
#
.out <- list(result = .res, design = design, alleles = .info$haplos$alleles, selected.haplotypes = .info$haplos$selected.haplotypes, resampling = resampling, orig.call = sys.call(), date = date(), reference.method = .info$haplos$reference.method, rows.dropped = .info$data$rows.dropped, HWE.res = .info$check$HWE.res, ntri.seq = .ntri.seq, loglike = .lratio.ut, score = .score.ut, var.cov = .var.cov.ut, info = .info, temp = list(o.var.covar.list = .o.var.covar.list, npars.0 = .npars.0), temp1 = .data)
class(.out) <- "haplin"	#
if(verbose) cat("Done\n")
#
if(printout){
## PRINTOUT AND PLOTTING, FOR THE CONVENIENCE OF THE USER:
	if(scoretest != "only"){
		plot(.out, reference = .info$haplos$reference.method)#
	}
#
	if(verbose) cat("\n#################################\n")
	#
	print(summary(.out, reference = .info$haplos$reference.method)) #
}

#
##
return(invisible(.out))
}
