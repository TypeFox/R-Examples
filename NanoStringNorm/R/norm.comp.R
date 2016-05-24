# The NanoStringNorm package is copyright (c) 2012 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

norm.comp <- function(x, anno, replicates = NULL,  CodeCount.methods = c('none', 'sum', 'geo.mean'), Background.methods = c('none','mean', 'mean.2sd','max'), SampleContent.methods = c('none','housekeeping.sum', 'housekeeping.geo.mean', 'total.sum','top.mean', 'top.geo.mean', 'low.cv.geo.mean'), OtherNorm.methods = c('none','quantile','zscore', 'rank.normal', 'vsn'),  histogram = FALSE, verbose = TRUE, icc.method = "mixed") { 

	if (!requireNamespace("lme4")) {
		stop("norm.comp: lme4 is required");
		}

	# get correct list item from xls or NSN output
	if (class(x) == 'NanoString') {
		x <- x[[1]];
		header <- x[[2]];
		} 
	else if (class(x) == 'NanoStringNorm') {
		x <- x$normalized.data;
		}

	if (colnames(x)[1] == "CodeClass") colnames(x)[1] = "Code.Class";

	# set the replicates to be the sample names if missing
	if (is.null(replicates)) {
		replicates <- colnames(x[!colnames(x) %in% c("Name", "Code.Class", "Accession")]);
		}

	# check the samples in the replicates
	x.sample.names <- colnames(x)[!colnames(x) %in% c("Name", "Code.Class", "Accession")];

	if (length(x.sample.names) != length(replicates)){
		stop("norm.comp: Replicates has a different length than Samples in x");
		}
	else if (length(duplicated) == 0) {
		stop("norm.comp: You did not specify any replicates.  All values are unique. ");
		}

	CodeCount.method.results <- NULL;
	Background.method.results <- NULL;
	SampleContent.method.results <- NULL;
	OtherNorm.method.results <- NULL;
	cv.pos.results <- NULL;
	cv.hk.results <- NULL;
	cv.end.results <- NULL;
	icc.anova.results <- NULL;
	icc.mixed.results <- NULL;
	#corr.group.results <- NULL;
	#cv.group.results <- NULL;

	# function for coefficient of variation with scaling of negative values
	get.cv <- function(x) {
		if(all(is.na(x))) return(NA);
		x <- na.omit(x);
		if(any(x < 0)) x <- x-min(x);
		sd(x)/mean(x) * 100;
		}

	get.mean.corr <- function(x, ignore.pc = .1, method = "pearson") {
		#print(ignore.pc)
		#if (all(x>0)) x[x < quantile(x, ignore.pc)] <- NA;
		cor.matrix <- cor(x, method = method, use = "pairwise.complete.obs");
		mean(cor.matrix[upper.tri(cor.matrix)]);
		}

	# function for intra class correlation
	get.icc <- function(x, replicates, method) {

		# which samples are actual replicates
		actual.replicates <- is.duplicated(replicates);
		
		# check missing, no groups left, only one group
		if (all(is.na(x))) return(NA);
		if (!anyDuplicated(replicates[actual.replicates & !is.na(x)])) return(NA);
		if (nlevels(as.factor(replicates[actual.replicates & !is.na(x)])) == 1) return(NA);

		# see package multilvel 
		if (method == "anova") {
			anova.data <- data.frame(x = x[actual.replicates], replicates = replicates[actual.replicates]);
			anova.fit <- summary(aov(x ~ as.factor(replicates), data = anova.data));
			ms_between <- anova.fit[[1]][1,3];
			ms_within  <- anova.fit[[1]][2,3];
			group.size <- (anova.fit[[1]][2,1] + (anova.fit[[1]][1,1] + 1))/(anova.fit[[1]][1,1] + 1);
			icc <- (ms_between - ms_within)/(ms_between + ((group.size - 1) * ms_within));
			#icc <- anova.fit[1,3]/(anova.fit[1,3]+anova.fit[2,3]);

			# low df models sometimes have problems 
			# returning negative variances
			# just drop these
			#if(is.na(icc)) browser()
			if (icc <= 0 | (ms_between - ms_within) == 0) {
				icc <- NA;
				}

			}
		if (method == "mixed") {
			#expr = lme(x ~ random = ~ 1|replicates)
			mixed.fit <- tryCatch(
				expr = lme4::lmer(x ~ 1|replicates),
				error = function(e) { return ( c(NULL)); }
				);

			if (is.null(mixed.fit)) {
				icc <- NA;
				}
			else {
				s2_between <- as.numeric(lme4::VarCorr(mixed.fit)$replicates[1,1]);
				s2_within <- as.numeric(attr(lme4::VarCorr(mixed.fit), "sc"))^2;
				icc <- s2_between/(s2_between + s2_within);
				}

			# low df models sometimes have problems 
			# returning negative variances which get rounded to zero
			# just drop these
			if (icc == 0) {
				icc <- NA;
				}
			}

		return(icc);
		}

	get.accuracy <- function(x, group, n.reps) {
		# split the data to get accuracy
		}

	# function to find duplicate *pairs* as logical vectors
	is.duplicated <- function(x) {
		duplicated(x) | duplicated(x, fromLast = TRUE);
		}
	
	is.unique <- function(x) {
		!(duplicated(x) | duplicated(x, fromLast = TRUE));
		}

	# loop over each normalization
	for (CodeCount.method in CodeCount.methods) {
		for (Background.method in Background.methods){
			for (SampleContent.method in SampleContent.methods) {
				for (OtherNorm.method in OtherNorm.methods) {
					if (OtherNorm.method %in% c("zscore", "rank.normal") & (CodeCount.method != "none" | SampleContent.method != "none" )) next;

					if (verbose == TRUE) {
						print(paste(CodeCount.method, Background.method, SampleContent.method, OtherNorm.method, sep="_"));
						}

					# normalize wrapped in tryCatch to avoid truncating output if step fails
					data.nsn <- tryCatch(
						expr = NanoStringNorm(
							x = x,
							CodeCount = CodeCount.method,
							Background = Background.method,
							SampleContent = SampleContent.method,
							OtherNorm = OtherNorm.method,
							take.log = TRUE,
							verbose = FALSE,
							return.matrix.of.endogenous.probes = FALSE
							),
						error = function(e) {return(NA)}
						)

					if(!all(is.na(data.nsn))) {
						#data.nsn <- data.nsn[[1]][grepl("[Ee]ndogenous", data.nsn[[1]]$Code.Class),-c(1:3)];
						data.nsn <- data.nsn[[1]][,-c(1:3)];

						if (OtherNorm.method %in% c("rank.normal","zscore")) {
							# calc missing
							missing.pc <- apply(data.nsn, 1, function(x) sum(is.na(x))/length(x));
							}
						else {
							# calc missing
							missing.pc <- apply(data.nsn, 1, function(x) sum(x==0)/length(x));
							# set missing to NA
							data.nsn[data.nsn == 0] <- NA;
							}

						# remove genes with lots of missing
						missing.gt90pc <- missing.pc > .90;

						if (all(missing.gt90pc)) {
							stop("norm.comp: all genes have greater than 90% missing.");
							}
						else if (any(missing.gt90pc)) {
							data.nsn[missing.gt90pc,] <- NA;
							}

						# calculate mean cv for pos, hk, end
						cv <- apply(
							X = data.nsn,
							MARGIN = 1,
							FUN = get.cv
							);

						is.90pc.missing <- apply(
							X = data.nsn,
							MARGIN = 1,
							FUN = function(x) sum(is.na(x)) >= round((0.9 * length(x)),0)
							);

						cv.pos <- tryCatch(
							expr = mean(cv[x$Name %in% c("POS_A(128)","POS_B(32)","POS_C(8)", "POS_D(2)", "POS_A", "POS_B", "POS_C", "POS_D")], na.rm = TRUE),
							error = function(e) { return ( c(NA)); }
							);

						cv.hk <- tryCatch(
							expr = mean(cv[x$Code.Class %in% c("Control", "Housekeeping")], na.rm = TRUE),
							error = function(e) { return ( c(NA)); }
							);

						cv.end <- tryCatch(
							expr = mean(cv[grepl("[Ee]ndogenous", x$Code.Class) & !is.90pc.missing], na.rm = TRUE),
							error = function(e) { return ( c(NA)); }
							);

						# only try icc if more than 1 group	
						if (nlevels(as.factor(replicates[is.duplicated(replicates)])) > 1) {

							icc.anova <- NA;
							icc.mixed <- NA;

							if (any(grepl("anova", icc.method))) {
								# fit model with and compare residual error
								icc.anova <- apply(
									X = data.nsn[grepl("[Ee]ndogenous", x$Code.Class),],
									MARGIN = 1,
									FUN = get.icc,
									replicates = replicates,
									method = "anova"
									);

								icc.anova <- median(icc.anova, na.rm = TRUE);
							}

							if (any(grepl("mixed", icc.method))) {
								icc.mixed <- apply(
									X = data.nsn[grepl("[Ee]ndogenous", x$Code.Class),],
									MARGIN = 1,
									FUN = get.icc,
									replicates = replicates,
									method = "mixed"
									);

								icc.mixed <- median(icc.mixed, na.rm = TRUE);

								#x.long <- matrix(as.vector(as.matrix(data.nsn[grepl("[Ee]ndogenous", x$Code.Class),])),ncol = 1);
								#rep.log <- rep(replicates, each = nrow(data.nsn[grepl("[Ee]ndogenous", x$Code.Class),]));

								#icc.mixed(x.long, rep.long, "mixed")
							}

	#						corr.group <- NA;
	#						cv.group <- NA;
							}
						else {

	#						corr.group <- get.mean.corr(data.nsn[,is.duplicated(replicates)], ignore.pc = .1);

	#						cv.group <- apply(
	#							X = data.nsn[,is.duplicated(replicates)],
	#							MARGIN = 1,
	#							FUN = get.cv
	#							);

	#						cv.group <- mean(cv.group, na.rm = TRUE); 

							icc.anova <- NA;
							icc.mixed <- NA;
							}
						}
					else {
						cv.pos <- NA;
						cv.hk <- NA;
						cv.end <- NA;
						icc.anova <- NA;
						icc.mixed <- NA;
						}

					# concatenate results
					CodeCount.method.results <- c(CodeCount.method.results, CodeCount.method);
					Background.method.results <- c(Background.method.results, Background.method);
					SampleContent.method.results <- c(SampleContent.method.results,SampleContent.method);
					OtherNorm.method.results <- c(OtherNorm.method.results,OtherNorm.method);

					cv.pos.results <- c(cv.pos.results, round(cv.pos,1));
					cv.hk.results <- c(cv.hk.results, round(cv.hk,1));
					cv.end.results <- c(cv.end.results, round(cv.end,1));
			
					icc.anova.results <- c(icc.anova.results, round(icc.anova,3));
					icc.mixed.results <- c(icc.mixed.results, round(icc.mixed,3));
					#corr.group.results <- c(corr.group.results, round(corr.group,1));
					#cv.group.results <- c(cv.group.results, round(cv.group,1));
					}
				}
			}
		}

	norm.comp.results <- data.frame(
		method = paste(CodeCount.method.results, Background.method.results, SampleContent.method.results, OtherNorm.method.results, sep="_"),
		CodeCount.method = CodeCount.method.results,
		Background.method = Background.method.results,
		SampleContent.method = SampleContent.method.results,
		OtherNorm.method = OtherNorm.method.results,
		cv.pos.results = cv.pos.results,
		cv.hk.results = cv.hk.results,
		cv.end.results = cv.end.results,
		cv.bio2tech.ratio = round(cv.end.results/cv.pos.results,2),
		#corr.group.results = corr.group.results,
		#cv.group.results = cv.group.results,
		icc.anova.results = icc.anova.results,
		icc.mixed.results = icc.mixed.results,
		stringsAsFactors = FALSE
		);

#	if (nlevels(as.factor(replicates[is.duplicated(replicates)])) > 1) {
#		norm.comp.results$corr.group.results <- NULL;
#		norm.comp.results$cv.group.results <- NULL;
#		}
#	else {
#		norm.comp.results$icc.mixed.results <- NULL;
#		norm.comp.results$icc.anova.results <- NULL;
#		}

	if (histogram == TRUE) {
		ns.green.rgb  <- rgb(193, 215, 66, maxColorValue = 255);
		ns.orange.rgb <- rgb(228, 108, 0, maxColorValue = 255);
		op.default <- par(cex = 1, cex.lab = 1.5, cex.axis = 1, cex.main = 1.5, las = 1, pch = 20, col.lab = 'grey30', col.axis = 'grey30', col.main = 'grey30', mar = c(5.1,5.1,4.1,2.1));
		hist(
			norm.comp.results$icc.mixed.results,
			breaks = 35,
			xlab = "ICC",
			freq = FALSE,
			main = "Replicate Intra Class Correlation",
			col = 'grey85',
			border = 'grey50',
			xlim = c(min(norm.comp.results$icc.mixed.results),1)
			);

		density.icc <- density(na.omit(norm.comp.results$icc.mixed.results));
		lines(density.icc, lwd = 4, col = ns.green.rgb);

		box(col = 'grey60', lwd = 3);

		}

	return(norm.comp.results);
	}
