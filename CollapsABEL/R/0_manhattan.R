#' Prepare data for Manhattan plot. 
#' 
#' @param chr integer. Chromosome vector.
#' @param bp integer. Position vector.
#' @param p numeric. P-value vector.
#' @param snp character. SNP name vector.
#' @param color_vec character/factor. Color vector. Doesn't have to be color names, any categorical variable will be fine.
#' @param sort_chr_bp logical. Whether to sort the whole data frame by CHR and BP before return.
#' @return A list with the following members (1) A data frame with columns including CHR, SNP, BP, P, etc. (2) Total number of SNPs. (3) A vector of unique chromosomes.
#' 
#' @author Kaiyin Zhong
#' @export
manhattanData = function(chr, bp, p, snp, color_vec = NULL, sort_chr_bp = TRUE) {
	# number of snps
	nsnps = length(chr)
	# vectors should have the same length
	stopifnot(
			length(bp) != nsnps || length(p) != nsnps || length(snp) != nsnps || length(color_vec) != nsnps
	)
	if(!is.numeric(chr)) {
		chr = as.integer(chr)
	}
	if(!is.numeric(bp)) {
		bp = as.integer(bp)
	}
	if(!is.numeric(p)) {
		p = as.integer(p)
	}
	if(!is.character(snp)) {
		snp = as.character(snp)
	}
	mlogp = mLogP(p)
	# add color vec if it's provided
	if(!is.null(color_vec)) {
		dat = data.frame(
				CHR = chr, 
				BP = bp, 
				P = p, 
				SNP = snp,
				MLOGP = mlogp,
				COLOR = color_vec
		)
	} else {
		dat = data.frame(
				CHR = chr, 
				BP = bp, 
				P = p, 
				SNP = snp,
				MLOGP = mlogp
		)
	}
	# order by chr and bp
	if(sort_chr_bp) dat = dat[order(dat$CHR, dat$BP), ]
	# update nsnps
	nsnps = nrow(dat)
	# unique chromosomes
	chr_unique = sort(unique(dat$CHR))
	# apparent chromosomes
	dat$ACHR = appChr(dat$CHR, nsnps, chr_unique)
	dat$XPOS = scaleByChr(dat$ACHR, dat$CHR, dat$BP, nsnps, chr_unique)
	list(
			data = dat,
			nsnps = nsnps,
			chr_unique = chr_unique
			)
}

mLogP = function(p) {
	# QC on p values. Should between 1e-300 and 1
	p[which(p < 1e-300)] = 1e-300
	p[which(p > 1)] = 1
	# -log of p
	-1 * log10(p)
}

# calculate apparent chr
appChr = function(chr, nsnps = NULL, chr_unique = NULL) {
	if(is.null(nsnps)) nsnps = length(chr)
	if(is.null(chr_unique)) chr_unique = sort(unique(chr))
	chr_diff = diff(chr)
	# if CHR is a sequence of natural numbers, then no change is needed
	if(all(chr_diff == 1)) {
		first_chr = chr_unique[1]
		if(first_chr == 1) return(chr)
		else {
			return(chr - first_chr + 1)
		}
	}
	achr = rep(0, nsnps)
	for(i in chr_unique) {
		achr = achr + (chr >= i)
	}
	achr
}

scaleByChr = function(achr, chr, bp, nsnps = NULL, chr_unique = NULL, zero_div_adjust = TRUE) {
	if(is.null(nsnps)) nsnps = length(chr)
	if(is.null(chr_unique)) chr_unique = sort(unique(chr))
	all_chr_scaled_pos = rep(NA, nsnps)
	for(chr_iter in chr_unique) {
		chr_check = chr_iter == chr
		chr_idx = which(chr_check)
		chr_nsnps = length(chr_idx)
		first_pos = bp[chr_idx[1]]
		pos_diff = bp[chr_idx] - first_pos
		scaled_pos = pos_diff / (pos_diff[chr_nsnps] + ifelse(zero_div_adjust, 0.05, 0))
		all_chr_scaled_pos[chr_idx] = scaled_pos
	}
	achr + all_chr_scaled_pos
}

#' Produce Manhattan plot
#' 
#' @param mh_dat_res list. Result from \code{manhattanData}
#' @param hlines numeric. Horizontal lines to draw.
#' @return ggplot object.
#' @importFrom ggplot2 ggplot geom_point scale_x_continuous aes xlab scale_color_discrete scale_y_continuous geom_hline ylab geom_segment annotate
#' 
#' @author Kaiyin Zhong
#' @export
manhattanPlot = function(mh_dat_res, hlines = NULL) {
	mh_data = mh_dat_res$data
	nsnps = mh_dat_res$nsnps
	chr_unique = mh_dat_res$chr_unique
	# limits of the y-axis
	mlogp_range = range(mh_data$MLOGP)
	minlogp = mlogp_range[1]
	maxlogp = ceiling(mlogp_range[2])
	# if nchr > 1, x axis labeled by scaled BP (sbp)
	# if nchr == 1, x axis labeled by BP
	if(length(chr_unique) > 1) {
		chr_color = TRUE
		myplot = ggplot(mh_data, aes(x=XPOS, y=MLOGP)) +
				xlab("Chromosomes") +
				scale_x_continuous(breaks=sort(unique(mh_data$ACHR)), minor_breaks=NULL, labels=chr_unique)
	} else {
		chr_color = FALSE
		myplot = ggplot(mh_data, aes(x=BP / 1e6, y=MLOGP)) +
				xlab(sprintf("Position on CHR %i (Mb)", mh_data$CHR[1]))
	}
	# do we need to use different colors for neighboring chromosomes? if there is only one chromosome, then no.
	if(chr_color) {
		if(!is.null(mh_data$COLOR)) {
			myplot = myplot + suppressWarnings(geom_point(aes(color = factor(paste(COLOR, ACHR %% 2))), alpha=.5)) +
					scale_color_discrete(name="")
		} else {
			myplot = myplot + suppressWarnings(geom_point(aes(color=factor(ACHR %% 2)), alpha=.5)) +
					scale_color_discrete(guide=FALSE)
		}
	} else {
		if(!is.null(mh_data$COLOR)) {
			myplot = myplot + suppressWarnings(geom_point(aes(color=factor(COLOR)), alpha=.5)) +
					scale_color_discrete(name = "")
		} else {
			myplot = myplot + suppressWarnings(geom_point(alpha=.5))
		}
	}
	# set y-axis limits
	myplot = myplot + scale_y_continuous(limits = c(minlogp, maxlogp), minor_breaks = NULL)
	# add horizontal lines if requested
	if(!is.null(hlines)) {
		hlines = -log10(hlines)
		for(hline in hlines) {
			myplot = myplot + geom_hline(yintercept = hline, alpha = .3, color = "blue")
		}
	}
	# set y-label
	myplot = myplot + ylab("-log P")
	myplot
}

#' Prepare data for \code{contrastPlot}
#' @param chr integer. Chromosome vector.
#' @param bp integer. Position vector.
#' @param p numeric. P-value vector.
#' @param gcdh_p numeric. GCDH p-value vector.
#' @param snp character. SNP name vector.
#' 
#' @author Kaiyin Zhong
#' @export
contrastData = function(chr, bp, p, gcdh_p, snp) {
	len = length(chr)
	stopifnot(
			length(bp) == len && length(p) == len && length(gcdh_p) == len
			)
	chr = rep(chr, 2)
	bp = rep(bp, 2)
	p = c(p, gcdh_p)
	color_vec = rep(c("Single SNP", "GCDH"), each = len)
	manhattanData(chr, bp, p, snp, color_vec, TRUE)
}

#' Produce contrast Manhattan plot
#' 
#' Overlays p-values from single-SNP method and GCDH. 
#' 
#' @param chr integer. Chromosome vector.
#' @param bp integer. Position vector.
#' @param p numeric. P-value vector.
#' @param gcdh_p numeric. GCDH p-value vector.
#' @param snp character. SNP name vector.
#' @param ... passed to \code{manhattanPlot}
#' @return ggplot object.
#' 
#' @author Kaiyin Zhong
#' @export
contrastPlot = function(chr, bp, p, gcdh_p, snp, ...) {
	cdata = contrastData(chr, bp, p, gcdh_p, snp)
	manhattanPlot(cdata, ...)
}

#annoContrastRegion = function(ggp, snp1, snp2) {
#	anno_dat = p$data[p$data$SNP == snp1, ]
#	anno_dat = anno_dat[anno_dat$colorvec == "GCDH", ]
#	sbp2 = bp2 * anno_dat$sbp / anno_dat$BP
#	anno_dat2 = within(anno_dat, {BP = bp2; SNP = snp2; sbp=sbp2})
#	p1 = p + geom_point(data = anno_dat2, aes(BP/1e6, mlogp), shape=3) +
#			geom_segment(aes(x = anno_dat$BP / 1e6,
#							y = anno_dat$mlogp,
#							xend = anno_dat2$BP / 1e6,
#							yend = anno_dat2$mlogp),
#					linetype=2)
#	p1
#}

#' Annotate a pair of SNPs in the contrast Manhattan plot
#' 
#' @param cplot ggplot object. The contrast Manhattan plot to be annotated.
#' @param snp1 character. First SNP. 
#' @param snp2 character. Second SNP.
#' @param linetype See \code{ggplot2::geom_segment}. Default to "dotted".
#' @param hjust  See \code{ggplot2::annotate}. Default to 0. 
#' @param text_size  See \code{ggplot2::annotate}. Default to 3. 
#' @return ggplot object.
#' 
#' @author Kaiyin Zhong
#' @export
connectSnpPair = function(cplot, snp1, snp2, linetype = "dotted", hjust = 0, text_size =  3) {
	ann_dat = cplot$data[cplot$data$SNP %in% c(snp1, snp2), ]
	top_x = left_x = ann_dat$BP[ann_dat$SNP == snp1][1] / 1e6
	right_x        = ann_dat$BP[ann_dat$SNP == snp2][1] / 1e6
	top_y   = ann_dat$MLOGP[ann_dat$SNP == snp1 & ann_dat$COLOR == "GCDH"]
	left_y  = ann_dat$MLOGP[ann_dat$SNP == snp1 & ann_dat$COLOR == "Single SNP"]
	right_y = ann_dat$MLOGP[ann_dat$SNP == snp2 & ann_dat$COLOR == "Single SNP"]
	data = data.frame(
			x = c(top_x, top_x), 
			y = c(top_y, top_y), 
			xend = c(left_x, right_x), 
			yend = c(left_y, right_y))
	cplot = cplot + geom_segment(
					data = data,
					mapping = aes(x = x,
							y = y,
							xend = xend,
							yend = yend
					), linetype = linetype, size = 0.1) +
			annotate("text", top_x, top_y, 
					label = paste(" ", snp1, "/", snp2, sep = ""), hjust = hjust, size = text_size)
	cplot
}




#' Manhattan plot the simple way
#' @param gcdh_report data.frame, from a GCDH analysis
#' @param outfile output image filepath. Any type (.png, .pdf, etc) supported by ggplot2::ggsave. Default to NULL. When it's not NULL, this function will try to save the plot to the specified path.
#' @return A ggplot object
#' 
#' @author kaiyin
#' @export
mh = function(gcdh_report, outfile=NULL) {
	cdata = contrastData(
			gcdh_report$CHR1, 
			gcdh_report$BP1, 
			gcdh_report$P1, 
			gcdh_report$P, 
			gcdh_report$SNP1
	)
	plot_res = manhattanPlot(cdata)
	if(!is.null(outfile)) {
		ggplot2::ggsave(outfile, plot_res, width = 11, height = 5)
	}
	plot_res
}
