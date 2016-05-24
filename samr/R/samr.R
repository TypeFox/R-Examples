# for robustness of code, we define constants to be used
#   later
samr.const.quantitative.response <- "Quantitative"
samr.const.twoclass.unpaired.response <- "Two class unpaired"
samr.const.survival.response <- "Survival"
samr.const.multiclass.response <- "Multiclass"
samr.const.oneclass.response <- "One class"
samr.const.twoclass.paired.response <- "Two class paired"
samr.const.twoclass.unpaired.timecourse.response <- "Two class unpaired timecourse"
samr.const.oneclass.timecourse.response <- "One class timecourse"
samr.const.twoclass.paired.timecourse.response <- "Two class paired timecourse"
samr.const.patterndiscovery.response <- "Pattern discovery"
samr.const.red.color <- 3
samr.const.green.color <- 10
samr.const.black.color <- 1
# Note: the table samr.xl.data.types created by the
# the Excel calling program has exactly these names as well
samr <- function(data, resp.type = c("Quantitative", 
	"Two class unpaired", "Survival", "Multiclass", "One class", 
	"Two class paired", "Two class unpaired timecourse", "One class timecourse", 
	"Two class paired timecourse", "Pattern discovery"), assay.type = c("array", 
	"seq"), s0 = NULL, s0.perc = NULL, nperms = 100, center.arrays = FALSE, 
	testStatistic = c("standard", "wilcoxon"), time.summary.type = c("slope", 
		"signed.area"), regression.method = c("standard", "ranks"), 
	return.x = FALSE,
	knn.neighbors = 10, random.seed = NULL, 
	nresamp = 20, nresamp.perm = NULL, xl.mode = c("regular", 
		"firsttime", "next20", "lasttime"), xl.time = NULL, xl.prevfit = NULL) {
	##SAM method. copyright june 2000: Goss, Tibshirani and
	#   Chu.
	## coded by r tibshirani;
	## y is response measure: 1,2 for two twoclass groups,
	## y=1,1,1 for onesample problem,
	## -1, 1,2-,2 for paired groups, with -1 paired with 1 etc
	## or survival time, or categorical for unordered groups
	#   1,2,3,..
	## quantitative for continuous ordered y#
	## timecourse data is also handled. In addition, pattern
	#   discovery is available- in which the eigengenes
	##  are used as a quantitative outcome
	##
	##
	## s0 is the exchangeability factor; you can specify
	## s0 as an actual value
	## s0.perc, the percentile of sd values to use for s0
	## or if both s0 and s0.perc are null (the default), then
	#   s0 is automatically estimated
	## returns
	##  evo= expected order statistics (length p=# of genes)
	## tt=numer/(sd+s0) test statistics on original data (and
	#   the ingredients)
	## ttstar0= p by nperms matrix of test statistics on
	#   permted data
	## ttstar= ttstar0 with columns sorted (largest values in
	#   row 1)
	## also returns permuted values: foldchange.star, ystar,
	#   sdstar, censoring.statusstar (for survival data)
	## in xl.mode firsttime or next 20, the function also
	#   returns x , which data$x, except for time-course data,
	##   where it is computed from in this function
	# from time summaries of data$x. However, this quantity is
	#   deleted the last time the function is called,
	#   as it is very large and not needed further
	this.call = match.call()
	resp.type.arg = match.arg(resp.type)
	assay.type = match.arg(assay.type)
	xl.mode = match.arg(xl.mode)
	if (!is.null(random.seed)) {
		set.seed(random.seed)
	}
	if (is.null(nresamp.perm)) {
		nresamp.perm = nresamp
	}
	nresamp.perm = min(nresamp, nresamp.perm)
	if (xl.mode == "regular" | xl.mode == "firsttime") {
		# initialize some things (harmlessly), just so that xl.mode
		#   will work correctly
		x = NULL
		xresamp = NULL
		ttstar0 = NULL
		evo = NULL
		ystar = NULL
		sdstar.keep = NULL
		censoring.status = NULL
		#censoring.status.star.keep = NULL	# Jun commented this line
		sdstar = NULL
		pi0 = NULL
		stand.contrasts = NULL
		stand.contrasts.star = NULL
		stand.contrasts.95 = NULL
		foldchange = NULL
		foldchange.star = NULL
		perms = NULL
		permsy = NULL
		eigengene = NULL
		eigengene.number = NULL
		##
		testStatistic <- match.arg(testStatistic)
		time.summary.type <- match.arg(time.summary.type)
		regression.method <- match.arg(regression.method)
		x = data$x
		y = data$y
		argy = y
		if (!is.null(data$eigengene.number)) {
			eigengene.number = data$eigengene.number
		}
		# impute missing data
		if (sum(is.na(x)) > 0) {
			require(impute)
			x = impute.knn(x, k = knn.neighbors)
			if (!is.matrix(x)) {
				x = x$data
			}
		}
		are.blocks.specified = FALSE
		# check that resp.type ok for seq data
		cond = (resp.type == "One class") | (resp.type == "Two class unpaired timecourse") | 
			(resp.type == "One class unpaired timecourse") | 
			(resp.type == "Two class paired timecourse") | (resp.type == 
			"Pattern discovery")
		if (assay.type == "seq" & cond) {
			stop(paste("Resp.type=", resp.type, " not allowed when assay.type='seq'"))
		}
		if (assay.type == "seq" & min(x) < 0) {
			stop(paste("Negative values not allowed when assay.type='seq'"))
		}
		
		## Jun added starts
		## check whether x are counts
		if (assay.type == "seq" & (sum(x%%1 != 0) != 0)) {
			stop("Non-integer values not alled when assay.type='seq'")
		}
		## Jun added ends
		
		# center columns of  array data if requested.
		# should not be allowed for seq data
		if (assay.type == "seq" & center.arrays) {
			stop(paste("Centering  not allowed when assay.type='seq'"))
		}
		if (assay.type == "seq" & regression.method == "ranks") {
			stop(paste("regression.method==ranks not allowed when assay.type='seq'"))
		}
		if (center.arrays) {
			x <- scale(x, center = apply(x, 2, median), scale = FALSE)
		}
		depth = scaling.factors = rep(NA, ncol(x))
		scaling.factors = (prod(depth)^(1/length(depth)))/depth
		# estimate seq depth and create resampled 3way array for
		#   seq data
		if (assay.type == "seq") {
			cat("Estimating sequencing depths...", fill = T)  ## Jun added this line
			depth = samr.estimate.depth(x)
			cat("Resampling to get new data matrices...", fill = T)  ## Jun added this line
			xresamp = resample(x, depth, nresamp = nresamp)
		}
		scaling.factors = (prod(depth)^(1/length(depth)))/depth
		# check if there are blocks for 2 class unpaired case
		if (resp.type == samr.const.twoclass.unpaired.response) {
			if (substring(y[1], 2, 6) == "Block" | substring(y[1], 
				2, 6) == "block") {
				junk = parse.block.labels.for.2classes(y)
				y = junk$y
				blocky = junk$blocky
				are.blocks.specified = TRUE
			}
		}
		# make sure 1,2, -1,1,, etc are non-character values coming
		#   from Excel
		if (resp.type == samr.const.twoclass.unpaired.response | 
			resp.type == samr.const.twoclass.paired.response | 
			resp.type == samr.const.oneclass.response | resp.type == 
			samr.const.quantitative.response | resp.type == samr.const.multiclass.response) {
			y = as.numeric(y)
		}
		# parse and summarize, if timecourse data
		sd.internal = NULL
		if (resp.type == samr.const.twoclass.unpaired.timecourse.response | 
			resp.type == samr.const.twoclass.paired.timecourse.response | 
			resp.type == samr.const.oneclass.timecourse.response) {
			junk = parse.time.labels.and.summarize.data(x, y, 
				resp.type, time.summary.type)
			y = junk$y
			x = junk$x
			sd.internal = sqrt(rowMeans(junk$sd^2))
			if (min(table(y)) == 1) {
				cat("", fill = T)
				cat("Warning: only one timecourse in one or more classes;\nSAM plot and FDRs will be unreliable; only gene scores are informative", 
				  fill = T)
			}
		}
		# if the data is timecourse, we have already summarized the
		#   time aspect.
		# Thus we change the resp.type to the appropriate
		#   non-time-course type. Note that the original value
		#  of resp.type was saved above in resp.type.arg
		if (resp.type == samr.const.twoclass.unpaired.timecourse.response) {
			resp.type = samr.const.twoclass.unpaired.response
		}
		if (resp.type == samr.const.twoclass.paired.timecourse.response) {
			resp.type = samr.const.twoclass.paired.response
		}
		if (resp.type == samr.const.oneclass.timecourse.response) {
			resp.type = samr.const.oneclass.response
		}
		stand.contrasts = NULL
		stand.contrasts.95 = NULL
		if (resp.type == samr.const.survival.response) {
			censoring.status = data$censoring.status
		}
		# do a thorough error  checking of the response data
		check.format(y, resp.type = resp.type, censoring.status = censoring.status)
		# transform to ranks if appropriate
		if (resp.type == samr.const.quantitative.response & regression.method == 
			"ranks") {
			y = rank(y)
			x = t(apply(x, 1, rank))
		}
		n <- nrow(x)
		ny <- length(y)
		sd <- NULL
		numer <- NULL
		# initial computation to get sd
		if (resp.type == samr.const.twoclass.unpaired.response & 
			testStatistic == "standard" & assay.type == "array") {
			init.fit <- ttest.func(x, y, sd = sd.internal)
			numer <- init.fit$numer
			sd <- init.fit$sd
		}
		if (resp.type == samr.const.twoclass.unpaired.response & 
			testStatistic == "wilcoxon" & assay.type == "array") {
			init.fit <- wilcoxon.func(x, y)
			numer <- init.fit$numer
			sd <- init.fit$sd
		}
		if (resp.type == samr.const.oneclass.response & assay.type == 
			"array") {
			init.fit <- onesample.ttest.func(x, y, sd = sd.internal)
			numer <- init.fit$numer
			sd <- init.fit$sd
		}
		if (resp.type == samr.const.twoclass.paired.response & 
			assay.type == "array") {
			init.fit <- paired.ttest.func(x, y, sd = sd.internal)
			numer <- init.fit$numer
			sd <- init.fit$sd
		}
		if (resp.type == samr.const.survival.response & assay.type == 
			"array") {
			init.fit <- cox.func(x, y, censoring.status)
			numer <- init.fit$numer
			sd <- init.fit$sd
		}
		if (resp.type == samr.const.multiclass.response & assay.type == 
			"array") {
			init.fit <- multiclass.func(x, y)
			numer <- init.fit$numer
			sd <- init.fit$sd
		}
		if (resp.type == samr.const.quantitative.response & assay.type == 
			"array") {
			init.fit <- quantitative.func(x, y)
			numer <- init.fit$numer
			sd <- init.fit$sd
		}
		if (resp.type == samr.const.patterndiscovery.response & 
			assay.type == "array") {
			init.fit <- patterndiscovery.func(x)
			numer <- init.fit$numer
			sd <- init.fit$sd
		}
		# for wilcoxon or rank regression or patterndiscovery , we
		#   set s0 to the 5th percentile of the sd values
		#   (automatic
		# estimation is not possible as the values of sd are too
		#   coarse)
		# also if dataset is small (< 500 genes), we don't attempt
		#   automatic
		# estimation of s0
		if ((resp.type == samr.const.quantitative.response & 
			(testStatistic == "wilcoxon" | regression.method == 
				"ranks" & assay.type == "array") | resp.type == 
			samr.const.patterndiscovery.response) | resp.type == 
			samr.const.twoclass.unpaired.response & assay.type == 
			"array" & testStatistic == "wilcoxon" | (nrow(x) < 
			500) & is.null(s0) & is.null(s0.perc)) {
			s0 = quantile(sd, 0.05)
			s0.perc = 0.05
		}
		# estimate s0 if necessary
		if (is.null(s0) & assay.type == "array") {
			if (!is.null(s0.perc)) {
				if ((s0.perc != -1 & s0.perc < 0) | s0.perc > 
				  100) {
				  stop("Illegal value for s0.perc: must be between 0 and 100, or equal\nto (-1) (meaning that s0 should be set to zero)")
				}
				if (s0.perc == -1) {
				  s0 = 0
				}
				if (s0.perc >= 0) {
				  s0 <- quantile(init.fit$sd, s0.perc/100)
				}
			}
			if (is.null(s0.perc)) {
				s0 = est.s0(init.fit$tt, init.fit$sd)$s0.hat
				s0.perc = 100 * sum(init.fit$sd < s0)/length(init.fit$sd)
			}
		}
		if (assay.type == "seq") {
			s0 = 0
			s0.perc = 0
		}
		# compute test statistics on original data
		##################
		# array type data
		if (resp.type == samr.const.twoclass.unpaired.response & 
			testStatistic == "standard" & assay.type == "array") {
			tt <- ttest.func(x, y, s0 = s0, sd = sd.internal)$tt
		}
		if (resp.type == samr.const.twoclass.unpaired.response & 
			testStatistic == "wilcoxon" & assay.type == "array") {
			tt <- wilcoxon.func(x, y, s0 = s0)$tt
		}
		if (resp.type == samr.const.oneclass.response & assay.type == 
			"array") {
			tt <- onesample.ttest.func(x, y, s0 = s0, sd = sd.internal)$tt
		}
		if (resp.type == samr.const.twoclass.paired.response & 
			assay.type == "array") {
			tt <- paired.ttest.func(x, y, s0 = s0, sd = sd.internal)$tt
		}
		if (resp.type == samr.const.survival.response & assay.type == 
			"array") {
			tt <- cox.func(x, y, censoring.status, s0 = s0)$tt
		}
		if (resp.type == samr.const.multiclass.response & assay.type == 
			"array") {
			junk2 <- multiclass.func(x, y, s0 = s0)
			tt = junk2$tt
			stand.contrasts = junk2$stand.contrasts
		}
		if (resp.type == samr.const.quantitative.response & assay.type == 
			"array") {
			tt <- quantitative.func(x, y, s0 = s0)$tt
		}
		if (resp.type == samr.const.patterndiscovery.response & 
			assay.type == "array") {
			junk <- patterndiscovery.func(x, s0 = s0, eigengene.number = eigengene.number)
			tt <- junk$tt
			eigengene = junk$eigengene
		}
		#seq data
		if (resp.type == samr.const.twoclass.unpaired.response & 
			assay.type == "seq") {
			junk = wilcoxon.unpaired.seq.func(xresamp, y)
			tt = junk$tt
			numer = junk$numer
			sd = junk$sd
		}
		if (resp.type == samr.const.twoclass.paired.response & 
			assay.type == "seq") {
			junk <- wilcoxon.paired.seq.func(xresamp, y)
			tt = junk$tt
			numer = junk$numer
			sd = junk$sd
		}
		if (resp.type == samr.const.quantitative.response & assay.type == 
			"seq") {
			junk <- quantitative.seq.func(xresamp, y)
			tt = junk$tt
			numer = junk$numer
			sd = junk$sd
		}
		if (resp.type == samr.const.survival.response & assay.type == 
			"seq") {
			junk <- cox.seq.func(xresamp, y, censoring.status)
			tt = junk$tt
			numer = junk$numer
			sd = junk$sd
		}
		if (resp.type == samr.const.multiclass.response & assay.type == 
			"seq") {
			junk2 <- multiclass.seq.func(xresamp, y)
			tt = junk2$tt
			numer = junk2$numer
			sd = junk2$sd
			stand.contrasts = junk2$stand.contrasts
		}
		###########
		# construct matrix of permutations
		if (resp.type == samr.const.quantitative.response | resp.type == 
			samr.const.multiclass.response | resp.type == samr.const.survival.response) {
			junk <- getperms(y, nperms)
			perms = junk$perms
			all.perms.flag = junk$all.perms.flag
			nperms.act = junk$nperms.act
		}
		if (resp.type == samr.const.twoclass.unpaired.response) {
			if (are.blocks.specified) {
				junk = compute.block.perms(y, blocky, nperms)
				permsy = matrix(junk$permsy, ncol = length(y))
				all.perms.flag = junk$all.perms.flag
				nperms.act = junk$nperms.act
			}
			else {
				junk <- getperms(y, nperms)
				permsy = matrix(y[junk$perms], ncol = length(y))
				all.perms.flag = junk$all.perms.flag
				nperms.act = junk$nperms.act
			}
		}
		if (resp.type == samr.const.oneclass.response) {
			if ((length(y) * log(2)) < log(nperms)) {
				allii = 0:((2^length(y)) - 1)
				nperms.act = 2^length(y)
				all.perms.flag = 1
			}
			else {
				nperms.act = nperms
				all.perms.flag = 0
			}
			permsy = matrix(NA, nrow = nperms.act, ncol = length(y))
			if (all.perms.flag == 1) {
				k = 0
				for (i in allii) {
				  junk = integer.base.b(i, b = 2)
				  if (length(junk) < length(y)) {
					junk = c(rep(0, length(y) - length(junk)), 
					  junk)
				  }
				  k = k + 1
				  permsy[k, ] = y * (2 * junk - 1)
				}
			}
			else {
				for (i in 1:nperms.act) {
				  permsy[i, ] = sample(c(-1, 1), size = length(y), 
					replace = TRUE)
				}
			}
		}
		if (resp.type == samr.const.twoclass.paired.response) {
			junk = compute.block.perms(y, abs(y), nperms)
			permsy = junk$permsy
			all.perms.flag = junk$all.perms.flag
			nperms.act = junk$nperms.act
		}
		if (resp.type == samr.const.patterndiscovery.response) {
			nperms.act = nperms
			perms = NULL
			permsy = NULL
			all.perms.flag = FALSE
		}
		# compute test statistics on permuted  data
		# sdstar.keep <- matrix(0,ncol=nperms.act,nrow=nrow(x)) ##
		#   Jun commented this line
		## Jun added starts
		sdstar.keep <- NULL
		if (assay.type != "seq") {
			sdstar.keep <- matrix(0, ncol = nperms.act, nrow = nrow(x))
		}
		## Jun added ends
		
		## Jun commented the following 4 lines
#		censoring.status.star.keep <- NULL
#		if (resp.type == samr.const.survival.response) {
#			censoring.status.star.keep <- matrix(0, ncol = nperms.act, 
#				nrow = length(y))
#		}
		
		ttstar <- matrix(0, nrow = nrow(x), ncol = nperms.act)
		foldchange.star = NULL
		if (resp.type == samr.const.twoclass.unpaired.response | 
			resp.type == samr.const.twoclass.paired.response) {
			foldchange.star <- matrix(0, nrow = nrow(x), ncol = nperms.act)
		}
		if (resp.type == samr.const.multiclass.response) 
		{
			stand.contrasts.star = array(NA, c(nrow(x), length(table(y)), 
				nperms.act))
		}
		# end of if(xltime=='regular' etc
	}
	if (xl.mode == "next20" | xl.mode == "lasttime") {
		# get stuff from prevfit
		evo = xl.prevfit$evo
		tt = xl.prevfit$tt
		numer = xl.prevfit$numer
		eigengene = xl.prevfit$eigengene
		eigengene.number = xl.prevfit$eigengene.number
		sd = xl.prevfit$sd - xl.prevfit$s0
		sd.internal = xl.prevfit$sd.internal
		ttstar = xl.prevfit$ttstar
		ttstar0 = xl.prevfit$ttstar0
		n = xl.prevfit$n
		pi0 = xl.prevfit$pi0
		foldchange = xl.prevfit$foldchange
		y = xl.prevfit$y
		x = xl.prevfit$x
		xresamp = xl.prevfit$xresamp
		censoring.status = xl.prevfit$censoring.status
		argy = xl.prevfit$argy
		testStatistic = xl.prevfit$testStatistic
		foldchange.star = xl.prevfit$foldchange.star
		s0 = xl.prevfit$s0
		s0.perc = xl.prevfit$s0.perc
		# ystar= xl.prevfit$ystar
		resp.type = xl.prevfit$resp.type
		resp.type.arg = xl.prevfit$resp.type.arg
		#censoring.status.star.keep = xl.prevfit$censoring.status.star.keep	# Jun commented this line
		assay.type = xl.prevfit$assay.type
		# sdstar= xl.prevfit$sdstar
		sdstar.keep = xl.prevfit$sdstar.keep
		resp.type = xl.prevfit$resp.type
		stand.contrasts = xl.prevfit$stand.contrasts
		stand.contrasts.star = xl.prevfit$stand.contrasts.star
		stand.contrasts.95 = xl.prevfit$stand.contrasts.95
		perms = xl.prevfit$perms
		permsy = xl.prevfit$permsy
		nperms = xl.prevfit$nperms
		nperms.act = xl.prevfit$nperms.act
		all.perms.flag = xl.prevfit$all.perms.flag
		depth = xl.prevfit$depth
		scaling.factors = xl.prevfit$scaling.factors
		nresamp = xl.prevfit$nresamp
		nresamp.perm = xl.prevfit$nresamp.perm
	}
	if (xl.mode == "regular") {
		first = 1
		last = nperms.act
	}
	if (xl.mode == "firsttime") {
		first = 1
		last = 1
	}
	if (xl.mode == "next20") {
		first = xl.time
		last = min(xl.time + 19, nperms.act - 1)
	}
	if (xl.mode == "lasttime") {
		first = nperms.act
		last = nperms.act
	}
	for (b in first:last) {
		cat(c("perm=", b), fill = TRUE)
		if (assay.type == "array") {
			xstar <- x
		}
		if (assay.type == "seq") {
			xstar <- xresamp[, , 1:nresamp.perm]
		}
		if (resp.type == samr.const.oneclass.response) {
			ystar = permsy[b, ]
			if (testStatistic == "standard") {
				ttstar[, b] <- onesample.ttest.func(xstar, ystar, 
				  s0 = s0, sd = sd.internal)$tt
			}
		}
		if (resp.type == samr.const.twoclass.paired.response) {
			ystar = permsy[b, ]
			if (assay.type == "array") {
				ttstar[, b] <- paired.ttest.func(xstar, ystar, 
				  s0 = s0, sd = sd.internal)$tt
				foldchange.star[, b] = foldchange.paired(xstar, 
				  ystar, data$logged2)
			}
			if (assay.type == "seq") {
				ttstar[, b] <- wilcoxon.paired.seq.func(xstar, 
				  ystar)$tt
				foldchange.star[, b] <- foldchange.seq.twoclass.paired(x, 
				  ystar, depth)  ## Jun added this line
			}
		}
		if (resp.type == samr.const.twoclass.unpaired.response) {
			ystar = permsy[b, ]
			if (assay.type == "array") {
				if (testStatistic == "standard") {
				  junk <- ttest.func(xstar, ystar, s0 = s0, sd = sd.internal)
				}
				if (testStatistic == "wilcoxon") {
				  junk <- wilcoxon.func(xstar, ystar, s0 = s0)
				}
				ttstar[, b] <- junk$tt
				sdstar.keep[, b] <- junk$sd
				foldchange.star[, b] = foldchange.twoclass(xstar, 
				  ystar, data$logged2)
			}
			if (assay.type == "seq") {
				ttstar[, b] <- wilcoxon.unpaired.seq.func(xstar, 
				  ystar)$tt
				foldchange.star[, b] <- foldchange.seq.twoclass.unpaired(x, 
				  ystar, depth)  ## Jun added this line
			}
		}
		if (resp.type == samr.const.survival.response) {
			o <- perms[b, ]
			if (assay.type == "array") {
				ttstar[, b] <- cox.func(xstar, y[o], censoring.status = censoring.status[o], 
				  s0 = s0)$tt
			}
			if (assay.type == "seq") {
				ttstar[, b] <- cox.seq.func(xstar, y[o], censoring.status = censoring.status[o])$tt
				#censoring.status.star.keep[, b] <- censoring.status[o]	# Jun commented this line
			}
		}
		if (resp.type == samr.const.multiclass.response) {
			ystar = y[perms[b, ]]
			if (assay.type == "array") {
				junk <- multiclass.func(xstar, ystar, s0 = s0)
				ttstar[, b] <- junk$tt
				sdstar.keep[, b] <- junk$sd
				stand.contrasts.star[, , b] = junk$stand.contrasts
			}
			if (assay.type == "seq") {
				junk <- multiclass.seq.func(xstar, ystar)
				ttstar[, b] <- junk$tt
				stand.contrasts.star[, , b] <- junk$stand.contrasts
			}
		}
		if (resp.type == samr.const.quantitative.response) {
			ystar = y[perms[b, ]]
			if (assay.type == "array") {
				junk <- quantitative.func(xstar, ystar, s0 = s0)
				ttstar[, b] <- junk$tt
				sdstar.keep[, b] <- junk$sd
			}
			if (assay.type == "seq") {
				junk <- quantitative.seq.func(xstar, ystar)
				ttstar[, b] <- junk$tt
			}
		}
		if (resp.type == samr.const.patterndiscovery.response) {
			xstar = permute.rows(x)
			junk <- patterndiscovery.func(xstar, s0 = s0, eigengene.number = eigengene.number)
			ttstar[, b] <- junk$tt
			sdstar.keep[, b] <- junk$sd
		}
		# end of for b in first:last
	}
	# sort columns of statistics from permuted samples, and
	#   compute expected order statistics
	if (xl.mode == "regular" | xl.mode == "lasttime") {
		ttstar0 <- ttstar
		for (j in 1:ncol(ttstar)) {
			ttstar[, j] <- -1 * sort(-1 * ttstar[, j])
		}
		for (i in 1:nrow(ttstar)) {
			ttstar[i, ] <- sort(ttstar[i, ])
		}
		evo <- apply(ttstar, 1, mean)
		evo <- evo[length(evo):1]
		#censoring.statusstar <- censoring.status.star.keep	## Jun commented this line
		sdstar <- sdstar.keep
		# estimation of pi0= prop of null genes
		pi0 = 1
		if (resp.type != samr.const.multiclass.response) {
			qq <- quantile(ttstar, c(0.25, 0.75))
		}
		if (resp.type == samr.const.multiclass.response) {
			qq <- quantile(ttstar, c(0, 0.5))
		}
		pi0 <- sum(tt > qq[1] & tt < qq[2])/(0.5 * length(tt))
		# compute fold changes, when applicable
		foldchange = NULL
		if (resp.type == samr.const.twoclass.unpaired.response & 
			assay.type == "array") {
			foldchange = foldchange.twoclass(x, y, data$logged2)
		}
		if (resp.type == samr.const.twoclass.paired.response & 
			assay.type == "array") {
			foldchange = foldchange.paired(x, y, data$logged2)
		}
		if (resp.type == samr.const.oneclass.response & assay.type == 
			"array") {
		}
		stand.contrasts.95 = NULL
		if (resp.type == samr.const.multiclass.response)
		{
			stand.contrasts.95 = quantile(stand.contrasts.star, 
				c(0.025, 0.975))
		}
		#if(assay.type=='seq'){foldchange=rep(NA,length(tt))} ##
		#   Jun commented this line
		# the last time through, unless otherwise specified, we
		#   delete x, since it is very big and is not needed
		#   further
		## Jun added starts
		if (resp.type == samr.const.twoclass.unpaired.response & 
			assay.type == "seq") {
			foldchange <- foldchange.seq.twoclass.unpaired(x, 
				y, depth)
		}
		if (resp.type == samr.const.twoclass.paired.response & 
			assay.type == "seq") {
			foldchange <- foldchange.seq.twoclass.paired(x, y, 
				depth)
		}
		## Jun added ends
		if (return.x == FALSE) {
			x = NULL
		}
	}
	
	return(list(n=n,
	x=x,
	xresamp=xresamp,
	y=y,
	argy=argy, 
	censoring.status= censoring.status,
	testStatistic=testStatistic,
	nperms=nperms,
	nperms.act=nperms.act,
	tt=tt,
	numer=numer,
	sd=sd+s0,
	sd.internal=sd.internal,
	s0=s0,
	s0.perc=s0.perc,
	evo=evo,
	perms=perms,
	permsy=permsy,
	nresamp=nresamp,
	nresamp.perm=nresamp.perm,
	all.perms.flag=all.perms.flag,
	ttstar=ttstar,
	ttstar0=ttstar0,
	eigengene=eigengene,
	eigengene.number=eigengene.number,
	pi0=pi0,
	foldchange=foldchange, 
	foldchange.star=foldchange.star,
	sdstar.keep=sdstar.keep,
	#censoring.status.star.keep=censoring.status.star.keep,	## Jun commented this line
	resp.type=resp.type,
	resp.type.arg=resp.type.arg,
	assay.type=assay.type,
	stand.contrasts=stand.contrasts,
	stand.contrasts.star=stand.contrasts.star,
	stand.contrasts.95=stand.contrasts.95,
	depth=depth,
	#scaling.factors=scaling.factors,	## Jun commented this line
	call=this.call))
} 
