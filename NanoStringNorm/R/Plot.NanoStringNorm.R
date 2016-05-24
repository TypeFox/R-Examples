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

Plot.NanoStringNorm <- function(x, plot.type = 'RNA.estimates', samples = NA, genes = NA, label.best.guess = TRUE, label.ids = list(), label.as.legend = TRUE, label.n = 10, title = TRUE, col = NA) {

	# check that the object being plotted is the right class
	if ( class(x) != 'NanoStringNorm' ) { stop('In order to plot the input object needs to be of class NanoStringNorm.  Try changing return.matrix.of.endogenous.probe to FALSE.'); }

	# check input plot.type
	if (!any(plot.type %in% c('all', 'mean.sd', 'cv', 'norm.factors', 'missing', 'volcano','batch.effects','RNA.estimates','positive.controls'))) {
		stop('Plot: Unrecognized plot.type.');
		} 

	# replace plot type all 
	if ('all' %in% plot.type) { plot.type <- c('mean.sd', 'cv', 'norm.factors', 'missing', 'volcano','batch.effects','RNA.estimates','positive.controls'); }

	if (all(is.na(col)) | length(col) != 2 ) {
		# setup default plotting parameters
	#	ns.green.rgb  <- rgb(193, 215, 46, maxColorValue = 255); # original
	#	ns.orange.rgb <- rgb(228, 108, 11, maxColorValue = 255);

		ns.green.rgb  <- rgb(193, 215, 66, maxColorValue = 255); # better greyscale 
		ns.orange.rgb <- rgb(228, 108, 0, maxColorValue = 255);

		col[1] <- ns.green.rgb;
		col[2] <- ns.orange.rgb;
		}

	op.default <- par(cex = 1, cex.lab = 1.5, cex.axis = 1, cex.main = 1.5, las = 1, pch = 20, col.lab = 'grey30', col.axis = 'grey30', col.main = 'grey30', mar = c(5.1,5.1,4.1,2.1));

	##############################################
	### gene. plot the mean vs sd of each gene ###
	##############################################

	if ('mean.sd' %in% plot.type) {
		# define controls
		is.control <- !grepl('endogenous', x$normalized.data$Code.Class,ignore.case=TRUE);
		is.positive.control <- grepl('pos', x$normalized.data$Code.Class,ignore.case=TRUE);
		is.negative.control <- grepl('neg', x$normalized.data$Code.Class,ignore.case=TRUE);
		is.hk.control <- grepl("control|housekeeping", x$normalized.data$Code.Class,ignore.case=TRUE);

		# setup plotting environment
		plot(
			x = x$gene.summary.stats.norm[grepl('Endogenous', x$normalized.data$Code.Class),'Mean'],
			y = x$gene.summary.stats.norm[grepl('Endogenous', x$normalized.data$Code.Class),'SD'],
			ylab = 'SD',
			xlab = 'Mean',
			main = if (title == TRUE) 'Gene: Mean vs Standard Deviation' else NA,
			xlim = c(0,max(x$gene.summary.stats.norm$Mean)),
			ylim = c(0,max(x$gene.summary.stats.norm$SD)),
			col = col[1],
			);

		# add the data points
		points(
			x = x$gene.summary.stats.norm[is.negative.control,'Mean'],
			y = x$gene.summary.stats.norm[is.negative.control,'SD'],
			col = col[2],
			pch = 17
			);
		points(
			x = x$gene.summary.stats.norm[is.positive.control,'Mean'],
			y = x$gene.summary.stats.norm[is.positive.control,'SD'],
			col = col[2],
			pch = 18
			);
		points(
			x = x$gene.summary.stats.norm[is.hk.control,'Mean'],
			y = x$gene.summary.stats.norm[is.hk.control,'SD'],
			col = col[2]
			);

		# add the lowess best fit line
		lines(lowess(x = x$gene.summary.stats.norm$Mean, y = x$gene.summary.stats.norm$SD), lwd = 4, col = 'grey60');
		box(col = 'grey60', lwd = 3);

		# add the legend
		legend(
			x = 'topleft',
			legend = c('Endogenous','Positive','Negative','Housekeeping'),
			col = c(col[1], col[2], col[2], col[2]),
			pch = c(16, 18,17,16),
			text.col = 'grey30',
#			lwd = 4,
			pt.cex = 1.5,
			bty = 'n'
			);


		# which genes are highly expressed and low sd (top 25% of mean and bottom 3% sd)
		is.high.mean <- x$gene.summary.stats.norm$Mean > quantile(x$gene.summary.stats.norm$Mean, .75) &  !is.control;
		sd.threshold <- quantile(x$gene.summary.stats.norm[is.high.mean,'SD'],.03);
		is.high.mean.and.low.sd <- is.high.mean & x$gene.summary.stats.norm$SD < sd.threshold;

		# which genes are outliers i.e. mean or sd is 3 sdev from mean
		SD.outliers <- (x$gene.summary.stats.norm$SD - mean(x$gene.summary.stats.norm$SD)) / sd(x$gene.summary.stats.norm$SD) > 3;

		# add best guess labels
		if (label.best.guess == TRUE) {
# & !'genes' %in% names(label.ids)) { 

			# high mean low sd
			if (any(is.high.mean.and.low.sd)) {
				text(
					x = x$gene.summary.stats.norm$Mean[is.high.mean.and.low.sd], 
					y = x$gene.summary.stats.norm$SD[is.high.mean.and.low.sd], 
					labels = gsub('hsa-','', rownames(x$gene.summary.stats.norm)[is.high.mean.and.low.sd]), 
					pos = c(1,2,3,4),
					col = 'grey30',
					cex = .7
					);
				}

			# outlier in sd
			if (any(SD.outliers)) {
				text(
					x = x$gene.summary.stats.norm$Mean[SD.outliers], 
					y = x$gene.summary.stats.norm$SD[SD.outliers], 
					labels = gsub('hsa-','',rownames(x$gene.summary.stats.norm)[SD.outliers]), 
					pos = c(1,2,3,4),
					col = 'grey30',
					cex = .7
					);
				}
			}

		# add explicit labels
		if ('genes' %in% names(label.ids)) {
			to.label <- rownames(x$gene.summary.stats.norm) %in% label.ids$genes;
			text(
				x = x$gene.summary.stats.norm$Mean[to.label], 
				y = x$gene.summary.stats.norm$SD[to.label], 
				labels = gsub('hsa-','', rownames(x$gene.summary.stats.norm)[to.label]), 
				pos = c(1,2,3,4),
				col = 'grey30',
				cex = .7
				);
			}

		}

	#############################################################################################
	### gene. plot the density of the coefficient of variation before and after normalization ###
	#############################################################################################
	
	if ('cv' %in% plot.type) {

		# which genes have expression counts greater than one after normalization.  removal helps with some strange values.
		expressed.genes <- x$gene.summary.stats.norm$Mean > 1;

		# get the density of the gene distributions
		density.raw  <- density(na.omit(x$gene.summary.stats.raw[expressed.genes, 'CV']), from = 0);
		density.norm <- density(na.omit(x$gene.summary.stats.norm[expressed.genes,'CV']), from = 0);

		# setup the plotting dimensions
		plot(
			x = 0,
			ylab = 'Density',
			xlab = 'CV %',
			main = if (title == TRUE) 'Gene: Coefficient of Variation' else NA,
			type = 'n',
			ylim = c(0, max(density.raw$y,density.norm$y)),
			xlim = c(0, max(density.raw$x,density.norm$x))
			);

		# add density the lines
		lines(density.raw, lwd = 4, col = col[1]); 
		lines(density.norm, lwd = 4, col = col[2]); 
		box(col = 'grey60', lwd = 3);

		# add the legend
		cv.pos.before <- round(mean(x$gene.summary.stats.raw[x$normalized.data$Name %in% c("POS_A(128)", "POS_B(32)", "POS_C(8)", "POS_D(2)","POS_A", "POS_B", "POS_C", "POS_D"), 'CV'], na.rm = TRUE),1);
		cv.pos.after <- round(mean(x$gene.summary.stats.norm[x$normalized.data$Name %in% c("POS_A(128)", "POS_B(32)", "POS_C(8)", "POS_D(2)","POS_A", "POS_B", "POS_C", "POS_D"), 'CV'], na.rm = TRUE),1);

		cv.hk.before <- round(mean(x$gene.summary.stats.raw[x$normalized.data$Code.Class %in% c("housekeeping", "Housekeeping", "Control", "control"), 'CV'], na.rm = TRUE),1);
		cv.hk.after <- round(mean(x$gene.summary.stats.norm[x$normalized.data$Code.Class %in% c("housekeeping", "Housekeeping", "Control", "control"), 'CV'], na.rm = TRUE),1);

		legend(
			x = 'topright',
			legend = c('', '', paste('Pos' ,cv.pos.before, '% -->', cv.pos.after, '%'), paste('HKs', cv.hk.before, '% -->', cv.hk.after, '%')),
			text.col = 'grey30',
			bty = 'n',
			xjust = 1
			);

		legend(
			x = 'topright',
			legend = c('Before Normalization', 'After Normalization'),
			col = c(col[1], col[2]),
			text.col = 'grey30',
			lwd = 4,
			bty = 'n'
			);
		}

	#########################################################################################################
	### gene. plot a volcano plot for the differential expression of traits supplied in the normalization ###
	#########################################################################################################

	if ('volcano' %in% plot.type & ncol(x$gene.summary.stats.norm) > 4) {
		# get the trait names
		trait.names <- gsub(
			pattern = 'FC_',
			replacement = '',
			x = colnames(x$gene.summary.stats.norm)[grepl('FC_', colnames(x$gene.summary.stats.norm))]
			);

		# loop over each trait and create a separate plot
		for (i in 1:length(trait.names)) {

			trait.name <- trait.names[i];

			# extract the pvalue and fold-change
			trait.p <- -log10(x$gene.summary.stats.norm[,paste('P_', trait.name, sep = '')]);
			trait.fc <- x$gene.summary.stats.norm[,paste('FC_', trait.name, sep = '')];

			# define the point color depending on magnitude of foldchange
			trait.col <- rep(col[1], length(trait.p));
			trait.col[abs(trait.fc) > 2] <- col[2];

			# what is the adjusted alpha for plotting a threshold line
			bonferroni.alpha <- -log10(0.05/nrow(x$gene.summary.stats.norm));

			# define the xlimits
			trait.xlim <- max(abs(trait.fc), na.rm = TRUE);
			if (trait.xlim < 2) {
				trait.xlim <- 2.1;
				}

			# which genes have the most significant pvalues.
			is.gene.trait.top <- rank(-trait.p, ties.method = 'first') <= label.n; 

			# remove duplicate labels
			if ('genes' %in% names(label.ids)) {
				label.ids$genes <- unique(label.ids$genes[!label.ids$genes %in% rownames(x$gene.summary.stats.norm)[is.gene.trait.top]]);
				}

			# define the ylimits.  if very large truncate the top of the plot
			trait.ylim.max <- max(trait.p, na.rm = TRUE);
			trait.ylim <- trait.ylim.max;
			if (trait.ylim > 10) {
				trait.ylim <- 11;
				trait.p[trait.p > 10] <- 10.5;
				}
			else if (trait.ylim > bonferroni.alpha & trait.ylim < 10) {
				trait.ylim <- trait.ylim + .5;
				}
			else if (trait.ylim > 1.8 & trait.ylim < bonferroni.alpha) {
				trait.ylim <- bonferroni.alpha + .5;
				}
			else {
				trait.ylim <- 1.8;
				}

			# scale function to change data to 0-1 for plotting
			scale.data <- function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
	(x - from[1])/diff(from) * diff(to) + to[1];
	}

			# size the points according to mean and sig and foldchange
			trait.cex <- rep(.5,length(trait.p));
			trait.cex[(trait.p > 1.3 | abs(trait.fc) > 2) & !is.na(trait.p)] <- 1 + scale.data(x$gene.summary.stats.norm[(trait.p > 1.3 | abs(trait.fc) > 2) & !is.na(trait.p), 'Mean'], c(0,2));

			# group counts
			n1 <- sum(x$traits[,trait.name]==1,na.rm = TRUE);
			n2 <- sum(x$traits[,trait.name]==2,na.rm = TRUE);

			fold.change.label <- ifelse(as.logical(x$normalization.workflow[['is.log']]) | as.logical(x$normalization.workflow[['take.log']]), expression(paste('-lo',g[2],'Fold-Change', sep = '')), 'Fold-change');

			# setup the plotting environment
			plot(
				x = trait.fc,
				y = trait.p,
				xlab = fold.change.label,
				ylab = expression(paste('-lo',g[10],'P', sep = '')),
				main = if (title == TRUE) paste('Gene: Differential Expression\n',trait.name, '(n =',n2,'vs',n1,')') else NA,
				pch = 20,
				col = trait.col,
				cex = trait.cex,
#				xlim = c(-max(abs(trait.fc), na.rm = TRUE), max(abs(trait.fc), na.rm = TRUE)),
				xlim = c(-trait.xlim,trait.xlim),
				ylim = c(0, trait.ylim)
				);

			# add a more informative axis label if the pvalues are truncated
			if (trait.ylim == 11) {
				xlim.offset <- .05 * trait.xlim;
				lines(x = c(-trait.xlim,-trait.xlim + xlim.offset,-trait.xlim,-trait.xlim-xlim.offset,-trait.xlim) - (.04 * 2 * trait.xlim),y = seq(10,11,.25), lwd = 2, col = 'grey60', xpd = NA);
				axis(side = 2, at = 11, labels = ceiling(trait.ylim.max), col = 'grey30', col.ticks = 'black', cex.axis = 1.2);
				}

			is.gene.label.manual <- rownames(x$gene.summary.stats.norm) %in% label.ids$genes;

			# label points or as legend
			if (label.as.legend) {
				labels.best.guess <- 1:sum(is.gene.trait.top);

				if (length(label.ids$genes) > 0) {
					labels.manual <- (sum(is.gene.trait.top) + 1):(sum(is.gene.trait.top) + length(label.ids$genes));
					}
				else {
					labels.manual <- '';
					}

				# formatting the legend labels
				labels.for.legend <-rownames(x$normalized.data)[rownames(x$normalized.data) %in% c(rownames(x$gene.summary.stats.norm)[is.gene.trait.top], rownames(x$gene.summary.stats.norm)[is.gene.label.manual])];
				labels.for.legend <- gsub('hsa-', '', labels.for.legend);
				labels.for.legend <- labels.for.legend[order(trait.p[is.gene.trait.top | is.gene.label.manual], decreasing = TRUE)];
				labels.for.legend <- paste(1:length(labels.for.legend), ": ", labels.for.legend, sep = "");

				# don't add legend if missing pvalue
				labels.for.legend <- labels.for.legend[!is.na(trait.p[is.gene.trait.top | is.gene.label.manual][order(trait.p[is.gene.trait.top | is.gene.label.manual], decreasing = TRUE)])];

				legend.ylim.offset <- ifelse(trait.ylim == 11, 1, 0);

				# the label legend
				legend(
					x = -trait.xlim - trait.xlim * .14,
					y = trait.ylim - legend.ylim.offset,
					legend = labels.for.legend,
					bty = 'n',
					cex = .7,
					title = expression(bold("Top Genes:"))
					);
				}
			else {
				is.gene.label.manual <- rownames(x$gene.summary.stats.norm) %in% label.ids$genes;
				labels.best.guess <- gsub('hsa-','',rownames(x$gene.summary.stats.norm)[is.gene.trait.top]);
				labels.manual <- gsub('hsa-','',rownames(x$gene.summary.stats.norm)[is.gene.label.manual]);
			}

			# if best guess labeling then plot the top ranked genes
			if (label.best.guess == TRUE) {
				# change the relative position of the label depending on if using a legend
				if (label.as.legend) {
					pos.value <- NULL;
					}
				else {
					pos.value <- 1;
					}

				# label the best guess i.e. most sig
				text(
					x = trait.fc[is.gene.trait.top][order(trait.p[is.gene.trait.top], decreasing = TRUE)], 
					y = trait.p[is.gene.trait.top][order(trait.p[is.gene.trait.top], decreasing = TRUE)], 
					labels = labels.best.guess[1:sum(is.gene.trait.top)], 
					adj = ifelse(label.as.legend, .5, 1),
					pos = pos.value,
					col = 'grey30',
					cex = .6,
					srt = ifelse(label.as.legend,0,90)
					);
				}

			# label the genes explicitly
			if ('genes' %in% names(label.ids)) {
				text(
					x = trait.fc[is.gene.label.manual][order(trait.p[is.gene.label.manual], decreasing = TRUE)],
					y = trait.p[is.gene.label.manual][order(trait.p[is.gene.label.manual], decreasing = TRUE)], 
					labels = labels.manual,
					adj = ifelse(label.as.legend, .5, 1),
					pos = pos.value,
					col = 'grey30',
					cex = .6,
					srt = ifelse(label.as.legend,0,90)
					);
				}

			# add some lines
			abline(h = bonferroni.alpha, lty = 2, lwd = .5, col = 'grey30');
			abline(h = 1.3, lty = 2, lwd = .5, col = 'grey30');
			abline(v = -2, lty = 2, lwd = .5, col = 'grey30');
			abline(v = 2, lty = 2, lwd = .5, col = 'grey30');
			box(col = 'grey60', lwd = 3);

			}
		}

	################################################################################################
	### sample: mean raw vs norm missing. note that outliers could be a product of study design. ###
	################################################################################################

	if ('missing' %in% plot.type) {
		if (all(x$sample.summary.stats.norm$Sample.Missing == 0)) {
			no.missing <- TRUE;
			}
		else {
			no.missing <- FALSE;
			}

		# calculate the raw sample mean
		raw.sample.Mean <- apply(
			X = x$raw.data[grepl('Endogenous',x$raw.data$Code.Class), !colnames(x$raw.data) %in% c('Code.Class','Name', 'Accession')],
			MARGIN = 2,
			FUN = mean,
			na.rm = TRUE
			);
		
		# plot the data
		plot(
			x = raw.sample.Mean,
			y = x$sample.summary.stats.norm$Sample.Missing,
			xlab = 'Mean of Samples in Raw Data',
			ylab = 'Proportion of Missing in Normalized Data',
			main = if (title == TRUE) 'Sample: Missing' else NA,
			ylim = c(0, 1),
			xlim = c(0, max(raw.sample.Mean, na.rm = TRUE)),
			col = col[1]
			);
		
		if (no.missing == FALSE) {
			# what samples are outliers for missing
			outlier.missing <- (x$sample.summary.stats.norm$Sample.Missing - mean(x$sample.summary.stats.norm$Sample.Missing, na.rm = TRUE)) / sd(x$sample.summary.stats.norm$Sample.Missing, na.rm = TRUE);
			outlier.missing.threshold.pos <- 3 * sd(x$sample.summary.stats.norm$Sample.Missing, na.rm = TRUE) + mean(x$sample.summary.stats.norm$Sample.Missing, na.rm = TRUE);
			outlier.missing.threshold.neg <- -3 * sd(x$sample.summary.stats.norm$Sample.Missing, na.rm = TRUE) + mean(x$sample.summary.stats.norm$Sample.Missing, na.rm = TRUE);
			}			
		else {
			outlier.missing <- rep(0, nrow(x$sample.summary.stats.norm));
			outlier.missing.threshold.pos <- 0;
			outlier.missing.threshold.neg <- 0;
			}
		
		# what samples are outliers for mean
		outlier.mean <- (raw.sample.Mean - mean(raw.sample.Mean, na.rm = TRUE)) / sd(raw.sample.Mean, na.rm = TRUE);
		outlier.mean.threshold.pos <- 3 * sd(raw.sample.Mean, na.rm = TRUE) + mean(raw.sample.Mean, na.rm = TRUE);
		outlier.mean.threshold.neg <- -3 * sd(raw.sample.Mean, na.rm = TRUE) + mean(raw.sample.Mean, na.rm = TRUE);

		# if best guess label then label outliers
		if (label.best.guess == TRUE) {

			# label the samples with lots of missing.
			if (any(abs(outlier.missing) > 3)) {
				text(
					x = raw.sample.Mean[abs(outlier.missing) > 3], 
					y = x$sample.summary.stats.norm$Sample.Missing[abs(outlier.missing) > 3], 
					labels = rownames(x$sample.summary.stats.norm[abs(outlier.missing) > 3,]), 
					pos = c(1,2,3,4),
					col = 'grey30',
					cex = .7
					);
				}

			# label the samples with low mean
			if (any(abs(outlier.mean) > 3) ) {
# & !'samples' %in% names(label.ids)) {
				text(
					x = raw.sample.Mean[abs(outlier.mean) > 3], 
					y = x$sample.summary.stats.norm$Sample.Missing[abs(outlier.mean) > 3], 
					labels = rownames(x$sample.summary.stats.norm[abs(outlier.mean) > 3,]), 
					pos = c(1,2,3,4),
					col = 'grey30',
					cex = .7
					);
				}
			}
		# explicitly label the samples
		if ('samples' %in% names(label.ids)) {
			is.gene.manual.label <- rownames(x$sample.summary.stats.norm) %in% label.ids$samples;
			text(
				x = raw.sample.Mean[is.gene.manual.label], 
				y = x$sample.summary.stats.norm$Sample.Missing[is.gene.manual.label], 
				labels = rownames(x$sample.summary.stats.norm[is.gene.manual.label,]), 
				pos = c(1,2,3,4),
				col = 'grey30',
				cex = .7
				);
			}

		# add some lines
		abline(h = outlier.missing.threshold.pos, lty = 2, lwd = .5, col = 'grey30');
		abline(h = outlier.missing.threshold.neg, lty = 2, lwd = .5, col = 'grey30');
		abline(v = outlier.mean.threshold.pos, lty = 2, lwd = .5, col = 'grey30');
		abline(v = outlier.mean.threshold.neg, lty = 2, lwd = .5, col = 'grey30');
		abline(lowess(x = raw.sample.Mean, y = x$sample.summary.stats.norm$Sample.Missing), lwd = 4, col = 'grey60');
		box(col = 'grey60', lwd = 3);

		}

	#####################################################################################################################################################
	### sample. evaluate if sampleContent estimates are consistent.  HK genes do not use ligation for some assays.  outliers could indicate problems. ###
	#####################################################################################################################################################

	if ('RNA.estimates' %in% plot.type) {

	if (!any(grepl('[Hh]ousekeeping',x$raw.data$Code.Class))) {
		cat("Plot.NanoStringNorm: No housekeeping genes so RNA.estimates plot is skipped\n");
		}
	else {
		
		# get the geometric mean of the HK and endogenous genes
		RNA.hk <- apply(
			X = x$raw.data[grepl('[Hh]ousekeeping',x$raw.data$Code.Class), !colnames(x$raw.data) %in% c('Code.Class','Name', 'Accession')],
			MARGIN = 2,
			FUN = get.geo.mean
			);

		RNA.top <- apply(
			X = x$raw.data[grepl('[Ee]ndogenous',x$raw.data$Code.Class) & x$gene.summary.stats.norm$Mean > quantile(x$gene.summary.stats.norm$Mean,.8), !colnames(x$raw.data) %in% c('Code.Class','Name', 'Accession')],
			MARGIN = 2,
			FUN = get.geo.mean
			);

		# plot a scatterplot of the points
		plot(
			x = RNA.top,
			y = RNA.hk,
			xlab = 'Top Expressed Genes in Raw Data',
			ylab = '',
			#ylab = 'Housekeeping Genes in Raw Data',
			main = if (title == TRUE) 'Sample: RNA Content Estimates' else NA,
			col = col[1]
			)

		title(ylab = 'Housekeeping Genes in Raw Data', line = 3.5);

		# fit a linear model
		fit.RNA <- lm(RNA.hk ~ RNA.top);

		# what samples are outliers in terms of residuals
		fit.RNA.outliers <- resid(fit.RNA) / sd(resid(fit.RNA));

		# add the formula of the best fit line to the title
		mtext(paste('y =', round(fit.RNA$coef[2], 2), 'x +', round(fit.RNA$coef[1], 2)), side = 3, line = .2);

		# add the best fit line
		abline(fit.RNA, lwd = 4, col = 'grey60');
		box(col = 'grey60', lwd = 3);

		# if best buess label then label outliers in terms of residuals
		if (label.best.guess == TRUE) {
# & !'samples' %in% names(label.ids)) {
			if (any(abs(fit.RNA.outliers) > 3)) {
				text(
					x = RNA.top[abs(fit.RNA.outliers) > 3], 
					y = RNA.hk[abs(fit.RNA.outliers) > 3], 
					labels = rownames(x$sample.summary.stats.norm[abs(fit.RNA.outliers) > 3,]), 
					pos = c(1,2,3,4),
					col = 'grey30',
					cex = .7
					);
				}
			}
		# label samples explicitly
		if ('samples' %in% names(label.ids)) {
			to.label <- rownames(x$sample.summary.stats.norm) %in% label.ids$samples;
			text(
				x = RNA.top[to.label], 
				y = RNA.hk[to.label], 
				labels = rownames(x$sample.summary.stats.norm[to.label,]), 
				pos = c(1,2,3,4),
				col = 'grey30',
				cex = .7
				);
			}
		}
		}

	##############################################
	### sample. batch effects ####################
	##############################################

	#if ('batch.effects' %in% plot.type & dim(x$batch.effects)[1]>6) {
	if ('batch.effects' %in% plot.type) {
		if (all(is.na(x$batch.effects))) {
			cat("Plot.NanoStringNorm: no traits or plates for batch effects.  skipping plot.\n");
			}
		else {

			trait.names <- unique(x$batch.effects$trait.name);
			n.traits <- length(trait.names);
			
			sample.statistics <- c('Mean','SD','Missing');
			if (x$normalization.workflow['CodeCount'] != 'none' ) {
				sample.statistics <- c(sample.statistics,'PositiveControls');
				}
			if (x$normalization.workflow['Background'] != 'none' ) {
				sample.statistics <- c(sample.statistics,'NegativeControls');
				}
			if (x$normalization.workflow['SampleContent'] != 'none' ) {
				sample.statistics <- c(sample.statistics,'RNA Content');
				}

			if (length(trait.names) > 10) {
				op.batch.effects <- par(mfcol = c(3,1), mar=c(2,1,2,0), oma = c(6,5,5,1));
				}
			else {
				op.batch.effects <- par(mfcol = c(3,2), mar=c(2,1,2,2), oma = c(6,5,5,1));
				}

			for (sample.statistic in sample.statistics) {
				batch.data <- x$batch.effects[x$batch.effects$sample.statistics %in% sample.statistic,];
				batch.col <- rep(col[2], n.traits);
				batch.col[batch.data$p.ttest < 0.05] <- col[1];

				batch.cex <- -log10(batch.data$p.ttest);
				batch.cex[batch.data$p.ttest > 0.05] <- 1;
			#batch.cex[batch.data$p.ttest < 0.05] <- punif(batch.cex[batch.data$p.ttest > 0.05]);
				#batch.cex[batch.data$p.ttest < 0.05] <- punif(batch.cex[batch.data$p.ttest > 0.05]);

				batch.diff <- batch.data$mean.grp2 - batch.data$mean.grp1;

				plot(
					x = 1:n.traits,
					y = batch.diff,
					cex = batch.cex,
					col = batch.col,
					xlab = sample.statistic,
					ylim = c(-max(abs(batch.diff), na.rm=TRUE),max(abs(batch.diff), na.rm=TRUE)),
					#main = if (title == TRUE) sample.statistic else NA,
					main = sample.statistic,
					xaxt = 'n'
					);

				axis(1, at = 1:n.traits, labels = NA, col.axis = 'grey30');
				size.x.labels <- ifelse(max(nchar(as.character(trait.names))>10),.7,1);

				if (length(trait.names) > 10) {
					if (sample.statistic == 'Missing' | sample.statistic == 'RNA Content') {
						if (title == TRUE) mtext('Sample: Batch Effects', side = 3, cex = 2, col = 'grey30', outer = TRUE, line = .8);
						mtext('Trait', side = 1, cex = 1.5, col = 'grey30', outer = TRUE, line = 4.5);
						mtext('Mean of Group2 relative to Group1', side = 2, cex = 1.5, col = 'grey30', outer = TRUE, line = 2.3, las = 0);
						axis(1, at = 1:n.traits, labels = trait.names, col.axis = 'grey30', las = 3, cex.axis=size.x.labels);
						}
					}
				else{
					if (sample.statistic == 'Missing' | sample.statistic == 'RNA Content') {
						axis(1, at = 1:n.traits, labels = trait.names, col.axis = 'grey30', las = 3);
						}
					if (sample.statistic == 'Missing') {
						if (title == TRUE) mtext('Sample: Batch Effects', side = 3, cex = 2, col = 'grey30', outer = TRUE, line = .8);
						mtext('Trait', side = 1, cex = 1.5, col = 'grey30', outer = TRUE, line = 4.5);
						mtext('Mean of Group2 relative to Group1', side = 2, cex = 1.5, col = 'grey30', outer = TRUE, line = 2, las = 0);
						}
					}

				abline(h = 0, lwd = 2, lty =  1, col = 'grey30');
				box(col = 'grey60', lwd = 3);

				}
			par(op.batch.effects);
			}
		}
	
	######################################
	### sample. plot the raw controls  ###
	######################################

#	if ('raw.controls' %in% plot.type) {
#
#		# setup the plotting environment
#		op.multi.plot <- par(mfrow = c(2,1), mar=c(1,0,1,0), oma = c(5,5,5,2));
#
#		# how many plots are needed
#		n.plots <- ceiling(nrow(x$raw.data)/24);
#
#		codeClasses <- c("Positive", "Negative", "Housekeeping", "Endogenous");
#
#		# loop over codeclass
#		for (codeClass in codeClasses) {
#
#			# loop over each set of 24 samples
#			for (n.plot in 1:n.plots) {
#
#				# get the start and stop rows for plotting in bins of 24 samples.  special exception for the last bin.
#				first.sample.to.plot <- (n.plot - 1) * 24 + 1; 
#				last.sample.to.plot  <- first.sample.to.plot + 24;
#				if (n.plot == n.plots) last.sample.to.plot <- ncol(x$raw.data); 
#
#				# which samples are to be plotted
#				samples.to.plot <- rep(FALSE, ncol(x$raw.data));
#				samples.to.plot[first.sample.to.plot:last.sample.to.plot] <- TRUE;
#				genes.to.plot <- x$data.raw[grepl(codeClass, x$data.raw), samples.to.plot];
#
#				# setup plotting environment
#				boxplot(
#					x = genes.to.plot,
#					ylab = 'Raw Counts',
#					xlab = 'Samples',
#					xaxt = 'n',
#					);
#
#				# what samples are outliers i.e. greater than 100% from the mean ~ 3sd
#				outlier.samples <- abs(
#					normalization.factors.to.plot.scaled.n48[,'positiveControls']) > 100 | 
#					abs(normalization.factors.to.plot.scaled.n48[,'negativeControls']) > 100 |
#					abs(normalization.factors.to.plot.scaled.n48[,'sampleContent']) > 100;
#				
#				# if best guess label the outliers
#				if (label.best.guess == TRUE & !'samples' %in% names(label.ids)) {
#					outlier.samples.labels <- names(outlier.samples);
#					outlier.samples.labels[outlier.samples == FALSE] <- NA;
#					axis(side = 1, at = seq(2, 4*48, by = 4), labels = FALSE, col.axis = 'grey30', las = 3, cex.axis = .8);
#					axis(side = 1, at = seq(2, 4*48, by = 4), labels = outlier.samples.labels, col.axis = 'grey30', las = 3, tick = FALSE, cex.axis = .5, line = -1.5, hadj = 0);
#					}
#				# explicitly label
#				else if ('samples' %in% names(label.ids)) {
#					to.label <- rownames(x$sample.summary.stats.norm) %in% label.ids$samples;
#					sample.labels <- rownames(x$sample.summary.stats.norm)[to.label];
#					sample.labels[!to.label] <- NA;
#					sample.labels <- c(sample.labels, rep(NA,48-length(sample.labels)));
#					axis(side = 1, at = seq(2, 4*48, by = 4), labels = FALSE, col.axis = 'grey30', las = 3, cex.axis = .8);
#					axis(side = 1, at = seq(2, 4*48, by = 4), labels = sample.labels, col.axis = 'grey30', las = 3, tick = FALSE, cex.axis = .5, line = -1.5, hadj = 0);
#					}
#
#				# add legend and axis labels to the first plot on every page
#				if (n.plot %in% seq(1,100,2)) {
#					legend(
#						x = 96, y = 190,
#						legend = c('Positive Controls', 'Negative Controls', 'RNA Sample Content'),
#						col = c(col[1], col[2], 'grey60'),
#						text.col = 'grey30',
#						lwd = 3,
#						cex = .9,
#						horiz = TRUE,
#						bty = 'n',
#						xpd = NA,
#						y.intersp = 1,
#						xjust = 0.5
#						);
#
#					if (title == TRUE) mtext('Sample: Normalization Parameters', side = 3, cex = 2, col = 'grey30', outer = TRUE, line = .8);
#					mtext('Samples', side = 1, cex = 1.5, col = 'grey30', outer = TRUE, line = 2);
#					mtext('Percent', side = 2, cex = 1.5, col = 'grey30', outer = TRUE, line = 3, las = 0);
#					}
#
#				# add some lines
#				abline(h =  100, lty = 2, lwd = .5, col = 'grey30');
#				abline(h = -100, lty = 2, lwd = .5, col = 'grey30');
#				abline(h =    0, lwd = 2, lty =  1, col = 'grey30');
#				box(col = 'grey60', lwd = 3);
#				}
#			}
#		# reset the plotting parameters
#		par(op.multi.plot);
#		}

	##############################################
	### sample. plot the normalization factors ###
	##############################################

	if ('norm.factors' %in% plot.type) {

		# skip if no normalization factors
		if (x$normalization.workflow['CodeCount'] == 'none' && x$normalization.workflow['Background'] == 'none' && x$normalization.workflow['SampleContent'] == 'none') {
			cat("Plot.NanoStringNorm: No normalization factors to plot.\n");
			}
		else {
			# for simplicity add collate all the relevent data into a new object
			normalization.factors.to.plot <- data.frame(
				positiveControls = if (x$normalization.workflow['CodeCount'] != 'none') x$sample.summary.stats.norm$pos.norm.factor else rep(1,nrow(x$sample.summary.stats.norm)), 
				negativeControls = if (x$normalization.workflow['Background'] != 'none') (mean(x$sample.summary.stats.norm$background.level) / x$sample.summary.stats.norm$background.level) else rep(1,nrow(x$sample.summary.stats.norm)), 
				sampleContent = if (x$normalization.workflow['SampleContent'] != 'none') x$sample.summary.stats.norm$sampleContent.norm.factor else rep(1,nrow(x$sample.summary.stats.norm)),
				row.names = rownames(x$sample.summary.stats.norm)
				);

			normalization.factors.to.plot <- as.matrix(normalization.factors.to.plot);
			#colnames(normalization.factors.to.plot) <- c('positiveControls', 'negativeControls','sampleContent');
			#rownames(normalization.factors.to.plot) <- rownames(x$sample.summary.stats.norm);

			# first convert normalizaton factors to percentages above and below mean.  remember a high NF reflects a low value.
			normalization.factors.to.plot.scaled <- normalization.factors.to.plot;
			normalization.factors.to.plot.scaled[normalization.factors.to.plot >= 1] <- 100 * (normalization.factors.to.plot[normalization.factors.to.plot >= 1] - 1);
			normalization.factors.to.plot.scaled[normalization.factors.to.plot < 1 ] <- -100 * (1/(normalization.factors.to.plot[normalization.factors.to.plot < 1]) - 1) ;

			# how many plots are needed
			n.plots <- ceiling(nrow(normalization.factors.to.plot)/48);

			# setup the plotting environment
			if (n.plots > 1) {
				op.multi.plot <- par(mfrow = c(2,1), mar=c(1,0,1,0), oma = c(5,5,5,2));
				legend.y.pos <- 190;
				segment.offset.y.pos <- 2;
				}
			else {
				op.multi.plot <- par(mfrow = c(1,1), mar=c(1,0,1,0), oma = c(5,5,5,2));
				legend.y.pos <- 170;
				segment.offset.y.pos <- 1;
				}

			# loop over each set of 48 samples
			for (n.plot in 1:n.plots) {

				# get the start and stop rows for plotting in bins of 48 samples.  special exception for the last bin.
				first.sample.to.plot <- (n.plot - 1) * 48 + 1; 
				last.sample.to.plot  <- first.sample.to.plot + 47;
				if (n.plot == n.plots) last.sample.to.plot <- nrow(normalization.factors.to.plot.scaled); 

				# which samples are to be plotted
				samples.to.plot <- rep(FALSE, nrow(normalization.factors.to.plot.scaled));
				samples.to.plot[first.sample.to.plot:last.sample.to.plot] <- TRUE;
				normalization.factors.to.plot.scaled.n48 <- subset(normalization.factors.to.plot.scaled, samples.to.plot);

				# for simplicity just add some dummy data to the end in order to plot 48 samples
				if (nrow(normalization.factors.to.plot.scaled.n48) < 48){
					normalization.factors.to.plot.scaled.n48 <- rbind(
						normalization.factors.to.plot.scaled.n48,
						matrix(NA, nrow = 48 - nrow(normalization.factors.to.plot.scaled.n48), ncol = 3)
						);
					}

				# setup plotting environment
				plot(
					x = NA,
					type = 'n',
					lwd = 6,
					ylab = 'Percent',
					xlab = 'Samples',
					xaxt = 'n',
					xlim = c(1, 4*48),
					ylim = c(-135, 135)
					);

				# plot each normalization parameter as a thick segment
				segments(
					x0 = seq(1, 4*48, by = 4),
					y0 = sign(normalization.factors.to.plot.scaled.n48[,'positiveControls']) * segment.offset.y.pos,
					x1 =  seq(1, 4*48, by = 4),
					y1 = normalization.factors.to.plot.scaled.n48[,'positiveControls'],
					col = col[1],
					lwd = 2.8
					);

				segments(
					x0 =  seq(2, 4*48, by = 4),
					y0 = sign(normalization.factors.to.plot.scaled.n48[,'negativeControls']) * segment.offset.y.pos,
					x1 =  seq(2, 4*48, by = 4),
					y1 = normalization.factors.to.plot.scaled.n48[,'negativeControls'],
					col = col[2],
					lwd = 2.8
					);

				segments(
					x0 = seq(3, 4*48, by = 4),
					y0 = sign(normalization.factors.to.plot.scaled.n48[,'sampleContent']) * segment.offset.y.pos,
					x1 = seq(3, 4*48, by = 4),
					y1 = normalization.factors.to.plot.scaled.n48[,'sampleContent'],
					col = 'grey60',
					lwd = 2.8
					);

				# what samples are outliers i.e. greater than 100% from the mean ~ 3sd
				outlier.samples <- abs(
					normalization.factors.to.plot.scaled.n48[,'positiveControls']) > 100 | 
					abs(normalization.factors.to.plot.scaled.n48[,'negativeControls']) > 100 |
					abs(normalization.factors.to.plot.scaled.n48[,'sampleContent']) > 100;
				
				# if best guess label the outliers
				if (label.best.guess == TRUE) {
# & !'samples' %in% names(label.ids)) {
					outlier.samples.labels <- names(outlier.samples);
					outlier.samples.labels[outlier.samples == FALSE] <- NA;
					axis(side = 1, at = seq(2, 4*48, by = 4), labels = FALSE, col.axis = 'grey30', las = 3, cex.axis = .8);
					axis(side = 1, at = seq(2, 4*48, by = 4), labels = outlier.samples.labels, col.axis = 'grey30', las = 3, tick = FALSE, cex.axis = .5, line = -1.5, hadj = 0);
					}
				# explicitly label
				if ('samples' %in% names(label.ids)) {
					to.label <- rownames(x$sample.summary.stats.norm) %in% label.ids$samples;
					sample.labels <- rownames(x$sample.summary.stats.norm)[to.label];
					sample.labels[!to.label] <- NA;
					sample.labels <- c(sample.labels, rep(NA,48-length(sample.labels)));
					axis(side = 1, at = seq(2, 4*48, by = 4), labels = FALSE, col.axis = 'grey30', las = 3, cex.axis = .8);
					axis(side = 1, at = seq(2, 4*48, by = 4), labels = sample.labels, col.axis = 'grey30', las = 3, tick = FALSE, cex.axis = .5, line = -1.5, hadj = 0);
					}

				# add legend and axis labels to the first plot on every page
				if (n.plot %in% seq(1,100,2)) {

					legend(
						x = 96, y = legend.y.pos,
						legend = c('Positive Controls', 'Negative Controls', 'RNA Sample Content'),
						col = c(col[1], col[2], 'grey60'),
						text.col = 'grey30',
						lwd = 3,
						cex = .9,
						horiz = TRUE,
						bty = 'n',
						xpd = NA,
						y.intersp = 1,
						xjust = 0.5
						);
					if (title == TRUE) mtext('Sample: Normalization Parameters', side = 3, cex = 2, col = 'grey30', outer = TRUE, line = .8);
					mtext('Samples', side = 1, cex = 1.5, col = 'grey30', outer = TRUE, line = 2);
					mtext('Percent', side = 2, cex = 1.5, col = 'grey30', outer = TRUE, line = 3, las = 0);
					}

				# add some lines
				abline(h =  100, lty = 2, lwd = .5, col = 'grey30');
				abline(h = -100, lty = 2, lwd = .5, col = 'grey30');
				abline(h =    0, lwd = 2, lty =  1, col = 'grey30');
				box(col = 'grey60', lwd = 3);
				}

			# reset the plotting parameters
			par(op.multi.plot);
			}
		}

	###################################################################################
	### sample. the relationship between the observed vs expected positive controls ###
	###################################################################################

	if ('positive.controls' %in% plot.type) {

		n.plots <- ncol(x$normalized.data) - 3;

		# setup the plotting environment
		if (n.plots > 6) {
			op.samples <- par(mfrow = c(4,3), mar=c(1,0,1,0), oma = c(5,5,5,2));
			y.legend.pos <- 6.9;
			x.legend.pos <- 14.5;
			}
		else {
			op.samples <- par(mfrow = c(3,2), mar=c(1,0,1,0), oma = c(5,5,5,2));
			y.legend.pos <- 5;
			x.legend.pos <- 6;
			}
		# setup the plotting envirnoment

		# the expected concentration
		pos.control.conc <- log2(c(.125, .5, 2, 8, 32, 128));

		# the y axis limit 
		max.obs.pos <- max(log2(x$raw.data[grepl('[Pp]ositive',x$raw.data$Code.Class),-c(1:3)]), na.rm = TRUE);

		# which plots should have x axis and labels (bottom of page)
		if (ncol(x$normalized.data) >= 15) {
			xlab.plots <- c(seq(13, ncol(x$normalized.data), by = 12), seq(14, ncol(x$normalized.data), by = 12), seq(15, ncol(x$normalized.data), by = 12));
			}
		else {
			xlab.plots <- c(10,11,12);
			}

		# loop over each sample
		for (i in 4:ncol(x$normalized.data)) {

			sample.name <- colnames(x$normalized.data)[i];

			# the observed counts
			pos.control.count <- rev(log2(x$raw.data[grepl('[Pp]ositive', x$raw.data$Code.Class),i]));

			# setup the plotting environment
			plot(
				x = pos.control.conc,
				y = pos.control.count,
				ylim = c(0, max.obs.pos),
				xlim = c(-4, 8),
				type = 'n',
				yaxt = if (i %in% seq(4, ncol(x$normalized.data), by = 3)) 's' else 'n',
				xaxt = 'n'
				);

			# add the positive control points
			points(
				x = pos.control.conc,
				y = pos.control.count,
				col = col[1],
				cex = 1.5
				);

			# add the negative control points at 0 to get some perspective. the background should be around 50-60 positive control counts (5.5-6 on log2).
			points(
				x = jitter(rep(0,length(x$raw.data[grepl('[Nn]egative', x$raw.data$Code.Class),i])), amount = .5),
				y = log2(x$raw.data[grepl('[Nn]egative', x$raw.data$Code.Class),i]),
				col = col[2],
				cex = 1.5
				);

			# add the x axis
			if (i %in% xlab.plots) {
				axis(side = 1, at = seq(-3,7,2))
				}

			# fit a line to the plot 
			fit.pos <- lm(pos.control.count ~ pos.control.conc);

			# add the best fit formula to the title
			mtext(text = paste(sample.name, '\n','y =', round(fit.pos$coef[2], 2), 'x +', round(fit.pos$coef[1], 2)), side = 3, col = 'grey30', cex = .7, line = -1);

			# add some lines
			abline(fit.pos, lwd = 2, col = 'grey60');
			box(col = 'grey60', lwd = 3);

			# only add the legend, main and axis labels to the first plot on a page
			if (i %in% seq(4,ncol(x$normalized.data),12)) {
				legend(
					#x = 14.5, y = 22.3,
					x = x.legend.pos, y = max.obs.pos + y.legend.pos,
					legend = c('Positive Controls', 'Negative Controls'),
					col = c(col[1], col[2]),
					text.col = 'grey30',
					lwd = 3,
					cex = 1.5,
					horiz = TRUE,
					bty = 'n',
					xjust = 0.5,
					xpd = NA
					);

				if (title == TRUE) mtext('Sample: Positive Controls', side = 3, outer = TRUE, col = 'grey30', cex = 2, line = 2);
				mtext(expression(paste('Expected Concentration lo',g[2], ' fM', sep = '')), side = 1, outer = TRUE,col = 'grey30', line = 2);
				mtext(expression(paste('Observed lo',g[2], ' Counts'), sep = ''), side = 2, outer = TRUE, las = 0, col = 'grey30',line = 3);
				}
		
			}

		# reset the plotting parameters
		par(op.samples);

		}

	# reset the plotting parameters
	par(op.default);
	invisible(1)
	}
