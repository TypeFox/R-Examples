f.suest <- function(reslist, main = "", ad.hoc = 0, debug = T, diag.plot = T){
##
## reslist IS A LIST OF RESULTS FROM haplin, RUN ON THE SAME DATA FILE
## (CURRENTLY USES ONLY OVERLAPPING LINES FROM EACH ESTIMATION)
##
#
## KEEP AND RESET SEED VALUE TO FACILITATE SIMULATIONS
.ex.seed <- exists(".Random.seed")
if(.ex.seed) .oldseed <- .Random.seed
set.seed(24)
if(.ex.seed) on.exit(.Random.seed <<- .oldseed)
#
##

.l <- length(reslist)
### if(.l < 2) stop("f.suest needs more than one estimation result!")
#
##  EXTRACT SCORE MATRICES
if(ad.hoc == 0) .score <- lapply(reslist, function(x) x$score$score)
if(ad.hoc > 0) .score <- lapply(reslist, function(x) x$temp$o.var.covar.list[[ad.hoc]]$score)
#
## EXTRACT NUMBER OF NUISANCE PARAMETERS (BETTER WAYS?)
if(ad.hoc == 0) .npars <- sapply(reslist, function(x) sum(x$info$haplos$selected.haplotypes))
if(ad.hoc > 0) .npars <- sapply(reslist, function(x) x$temp$npars.0[ad.hoc])
if(debug) cat("ad hoc!\n")
if(reslist[[1]]$info$model$design %in% c("cc", "cc.triad")) .npars <- .npars + 1
if(reslist[[1]]$info$model$xchrom & is.null(reslist[[1]]$info$variables$sel.sex)) .npars <- .npars + 1
if(debug) cat("...dettte ogsaa!\n")
#
## JOIN SCORES INTO A SINGLE MATRIX
.lines <- lapply(.score, function(x) dimnames(x)[["orig.lines"]])
.lines.common <- .lines[[1]]
if(.l >= 2) for (i in 2:.l){
	.lines.common <- intersect(.lines.common, .lines[[i]])
}
.lines.account <- c(sapply(.lines, length), common = length(.lines.common))
#
## 
.scoremat <- vector(.l, mode = "list")
.kill <- vector(.l, mode = "list")
.d <- rep(NA, .l)
for (i in 1:.l){
	.scoremat[[i]] <- .score[[i]][.lines.common, -(1:.npars[i]), drop = F]
	.kill[[i]] <- which(apply(.scoremat[[i]], 2, function(x) length(unique(x)) == 1))
	if(length(.kill[[i]]) > 0){
		.scoremat[[i]] <- .scoremat[[i]][, -.kill[[i]], drop = F]
		if(debug) cat("########\n#########KILL########\n###########\n")
	}
	.d[i] <- dim(.scoremat[[i]])[2] # NUMBER OF PARAMETERS TO BE TESTED
}
.scoremat <- do.call("cbind", .scoremat)
#
## EMPIRICAL COVARIANCE MATRIX FOR THE (TOTAL) SCORE:
.vcov.score.emp <- cov(.scoremat) * (dim(.scoremat)[1])
#
## LIST OF EXPECTED VARIANCE-COVARIANCE MATRICES
.vcov.list <- vector(.l, mode = "list")
for (i in 1:.l){
	if(F & length(.kill[[i]]) > 0){## LITT I TVIL HER: MEN .kill BOR VEL HAANDTERES SOM DE ANDRE NUISANCE-PARAMETRENE?
		## SVARET SER MEGET LIKT UT I NOEN SMAA TESTER
		.tempinf <- solve(reslist[[i]]$temp)
		.tempinf <- .tempinf[-(.npars[i] + .kill[[i]]), -(.npars[i] + .kill[[i]]), drop = F]
		.temp <- solve(.tempinf)
		.vcov.list[[i]] <- .temp[-(1:.npars[i]),-(1:.npars[i]), drop = F]
	}else{
		if(ad.hoc == 0).temp <- reslist[[i]]$var.cov$var.cov.0
		if(ad.hoc > 0) .temp <- reslist[[i]]$temp$o.var.covar.list[[ad.hoc]]$var.cov
		.vcov.list[[i]] <- .temp[-c(1:.npars[i], .npars[i] + .kill[[i]]),-c(1:.npars[i], .npars[i] + .kill[[i]]), drop = F]
	}
}


.f.start.stop <- function(d){
## COMPUTE STARTPOINTS AND ENDPOINTS IN AN EXTRACTION SEQUENCE
## NO, SIMPLER JUST TO MAKE THE ACTUAL SEQUENCE:
	.dc1 <- cumsum(d) # END POINTS
	.dc0 <- c(0, .dc1[-length(.dc1)]) + 1 # START POINTS
	.ut <- list(length(d), mode = "list")
	for (i in seq(along = .dc0)){
		.ut[[i]] <- .dc0[i]:.dc1[i]
	}
	return(.ut)
}


.f.insert <- function(matfull, matsublist){
## FUNCTION FOR INSERTING A LIST OF SUBMATRICES (BLOCK DIAGONALS) INTO A LARGER MATRIX
	.d <- sapply(matsublist, function(x)dim(x)[1])
	if(dim(matfull)[1] != sum(.d))stop()
	#
	.st <- .f.start.stop(.d)
	for(i in seq(along = matsublist)){
		matfull[.st[[i]], .st[[i]]] <- matsublist[[i]]
	}
	return(matfull)
}


.f.extract <- function(matfull, d){
## FUNCTION FOR EXTRACTING A BLOCK-DIAGONAL FROM A LARGER MATRIX
	if(sum(d) != dim(matfull)[1]) stop()
	.ut <- vector(length(d), mode = "list")
	.st <- .f.start.stop(d)
	for (i in seq(along = d)){
		.ut[[i]] <- matfull[.st[[i]], .st[[i]], drop = F]
	}
	return(.ut)
}



.vcov.score.emp.list <- .f.extract(.vcov.score.emp, .d)

.G <- .vcov.score.emp
.G[,] <- 0
.Glist <- vector(.l, mode = "list")
.err <- F

for (i in seq(along = .d)){
	.tmp <- .vcov.score.emp.list[[i]]
	.tmp <- try(solve(.tmp), silent = T) # KRYSS FINGRE!
	if(inherits(.tmp, "try-error")) {.err <- T; next}
	if(dim(.tmp)[1] == 1) .tmp <- sqrt(.tmp)
	else .tmp <- mroot(.tmp) # KRYSS FINGRE!
###	.vc <- .inf.list[[i]]
	.vc <- solve(.vcov.list[[i]])
	if(dim(.vc)[1] == 1) .vc <- sqrt(.vc)
	else .vc <- mroot(.vc)
	.fact <- try(.tmp %*% t(.vc), silent = T) ## NBNB! HER BLIR DET FEILMELDING HVIS .vc IKKE ER AV FULL RANK NAAR MAN BRUKER mroot I LINJEN OVER!
	## SKJER SJELDENT, MEN SKJEDDE I KJORING AV suest PAA JACOBSSONS "tempres/o.run.hap_slide.3.cc.RData", GEN "COL5A1"
	if(inherits(.fact, "try-error")) {.err <- T; next}
	.Glist[[i]] <- .fact
}

.f.test <- function(scoresum, vcov.list, matr = F){
	###.scoresum <- apply(score, 2, sum)
	.d <- sapply(vcov.list, function(x)dim(x)[1])
	.st <- .f.start.stop(.d)
	#
	if(matr){
		.sc.test <- vector(length(.d), mode = "list")
		.pv.chisq <- vector(length(.d), mode = "list")
		#
		for(i in seq(along = .d)){
			.tmpsc <- scoresum[,.st[[i]], drop = F]
			.tmp <- .tmpsc %*% vcov.list[[i]]
			###.sc.test[[i]] <- rowSums(.tmp * .tmpsc) # TRENGER IKKE OPPBEVARE DENNE
			.sc.test <- rowSums(.tmp * .tmpsc)
			.pv.chisq[[i]] <- pchisq(.sc.test, df = .d[i], lower.tail = F)
		}
		.pv.chisq <- do.call("cbind", .pv.chisq)
		return(.pv.chisq)
	}
	else{
	
		.sc.test <- rep(NA, length(.d))
		.pv.chisq <- rep(NA, length(.d))
		#
		for(i in seq(along = .d)){
			.sc.test[i] <- scoresum[.st[[i]]] %*% vcov.list[[i]] %*% scoresum[.st[[i]]]
		}
		.pv.chisq <- pchisq(.sc.test, df = .d, lower.tail = F)
		return(.pv.chisq)
	}
}

if(!.err){
	.G <- .f.insert(.G, .Glist)

	.vcov.score.emp.just <-  t(.G) %*% .vcov.score.emp %*% .G

	#print(round(cov2cor(.vcov.tot) / cov2cor(.vcov.tot.just), 3))
	#plot(jitter(cov2cor(.vcov.tot)), jitter(cov2cor(.vcov.tot.just)))
	
	if(dim(.vcov.score.emp.just)[1] == 1) .mroot <- sqrt(.vcov.score.emp.just)
	else .mroot <- mroot(.vcov.score.emp.just)
	
	.normsim <- mvrnorm(n = 100000, mu = rep(0, dim(.mroot)[2]), Sigma = diag(rep(1, dim(.mroot)[2])), empirical = T)
	.normsim <- .normsim %*% t(.mroot)
	
	
	
	if(F){# TESTING
		### .tull <- cov(.normsim)
		### cat("gugg\n")
		### print(round(.vcov.tot / .vcov.tot.just, 3))
		### print(round(.tull / .vcov.tot.just, 3))
	}
	
	### .normtest <- apply(.normsim, 1, .f.test, vcov.list = .vcov.list)
	
	.normtest <- .f.test(.normsim, vcov.list = .vcov.list, matr = T)
}

.obs.pval <- .f.test(colSums(.scoremat), vcov.list = .vcov.list)

## DEMONSTRATES THAT EXCLUDING A FEW LINES HAS A FAIRLY MARKED EFFECT RIGHT HERE IN THIS EXAMPLE
tullsc <- NULL
for (i in seq(along = .score)){
	.tulltmp <- .score[[i]][,-c(1:.npars[i]), drop = F]
	if(length(.kill[[i]]) > 0){
		.tulltmp <- .tulltmp[, -.kill[[i]], drop = F]
	}
	tullsc <- c(tullsc, colSums(.tulltmp))
}
.obs.pval.alt <- .f.test(tullsc, vcov.list = .vcov.list)

if(!.err){
	.min.normtest <- do.call("pmin", as.dframe(.normtest))
	
	.obs.pval.corr <- sum(.min.normtest < min(.obs.pval))/length(.min.normtest)
	.obs.pval.alt.corr <- sum(.min.normtest < min(.obs.pval.alt))/length(.min.normtest)
	
	if(diag.plot){
		.ant <- dim(.normtest)[2]
		.dimplot <- ceiling(sqrt(.ant))
		#
		if(.ant == 1){
			par(mfrow = c(1,1)) 
			hist(.normtest[1:1000,], main = main, nclass = 50)
		}else
		if(.ant <= 30) {
			par(mfrow = c(1,1))
			pairs(.normtest[1:200,], main = main)
			par(mfrow = rep(min(.dimplot, 5), 2))# cat("Forhindrer crash naar mange plot!\n")
			apply(.normtest[1:1000,], 2, hist, main = main, nclass = 50)
		}
	}
}else{## IF FAILS, USE BONFERRONI
	.obs.pval.corr <- min(1, min(.obs.pval)*length(.obs.pval))
	.obs.pval.alt.corr <- min(1, min(.obs.pval.alt)*length(.obs.pval.alt))
}

return(list(pval.obs = .obs.pval, pval.alt = .obs.pval.alt, pval.obs.corr = .obs.pval.corr, pval.alt.corr = .obs.pval.alt.corr, lines.account = .lines.account, bonferroni = .err, kill = .kill))




}
