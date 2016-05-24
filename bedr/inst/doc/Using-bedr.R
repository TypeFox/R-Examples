## ---- results = "hold", message = TRUE, eval = TRUE, errors = TRUE-------
# load bedr library
library("bedr");

## ---- results = "hold", message = TRUE, eval = TRUE, errors = TRUE-------

if (check.binary("bedtools")) {

	# get example regions
	index <- get.example.regions();
	a <- index[[1]];
	b <- index[[2]];

	# region validation
	is.a.valid  <- is.valid.region(a);
	is.b.valid  <- is.valid.region(b);
	a <- a[is.a.valid];
	b <- b[is.b.valid];
		
	# print
	cat(" REGION a: ", a, "\n");
	cat(" REGION b: ", b, "\n");
	}

## ---- results = "hold", message = TRUE, eval = TRUE, errors = TRUE-------
if (check.binary("bedtools")) {

	# check if already sorted
	is.sorted <- is.sorted.region(a);

	# sort lexographically
	a.sort <- bedr.sort.region(a);

	# sort naturally
	a.sort.natural <- bedr.sort.region(a, method = "natural");

	# sort - explicit call using primary API function bedr()
	b.sort <- bedr(
		engine = "bedtools", 
		input = list(i = b), 
		method = "sort", 
		params = ""
		);

	# print
	cat(" REGION a: ", a.sort, "\n");
	cat(" REGION b: ", a.sort.natural, "\n");
	cat(" REGION c: ", b.sort, "\n");
	}

## ---- results = "hold", message = TRUE, eval = TRUE, errors = TRUE-------
if (check.binary("bedtools")) {

	# check if already merged (non-overlapping regions)
	is.merged <- is.merged.region(a.sort);
	is.merged <- is.merged.region(b.sort);

	# merge
	a.merge <- bedr.merge.region(a.sort);

	# merge - explicit call using primary API function bedr()
	b.merge <- bedr(
		engine = "bedtools", 
		input = list(i = b.sort), 
		method = "merge", 
		params = ""
		);

	# print
	cat(" REGION a: ", a.merge, "\n");
	cat(" REGION b: ", b.merge, "\n");
	}

## ---- results = "hold", message = TRUE, eval = TRUE, errors = TRUE-------
if (check.binary("bedtools")) {

	# subtract
	a.sub1 <- bedr.subtract.region(a.merge, b.merge);

	# subtract - explicit call using primary API function bedr()
	a.sub2 <- bedr(
		input = list(a = a.merge, b = b.merge), 
		method = "subtract", 
		params = "-A"
		);

	# print
	cat(" REGION a - sub1: ", a.sub1, "\n");
	cat(" REGION a - sub2: ", a.sub2, "\n");
	}

## ---- results = "hold", message = TRUE, eval = TRUE, errors = TRUE-------
if (check.binary("bedtools")) {

	# check if present in a region
	is.region <- in.region(a.merge, b.merge);

	# or alternatively R-like in command
	is.region <- a.merge %in.region% b.merge

	# print
	cat(" is.region: ", is.region, "\n");
	}

## ---- results = "hold", message = TRUE, eval = TRUE, errors = TRUE-------
if (check.binary("bedtools")) {

	# intersect / join
	a.int1 <- bedr.join.region(a.merge, b.merge);
	a.int2 <- bedr(
		input = list(a = a.sort, b = b.sort), 
		method = "intersect", 
		params = "-loj -sorted"
		);

	# multiple join
	d <- get.random.regions(15, chr="chr1", sort = TRUE);
	a.mult <- bedr.join.multiple.region(
		x = list(a.merge, b.merge, bedr.sort.region(d))
		);

	# print
	cat(" REGION a intersect: \n"); print(a.int1); cat("\n");
	cat(" REGION multi (a,b,c) intersect: \n"); print(a.mult); cat("\n");
	}

## ---- results = "hold", message = TRUE, eval = TRUE, errors = TRUE-------
if (check.binary("bedtools")) {

	# compare a and b set of sequences
	jaccard.stats <- jaccard(a.sort, b.sort);
	reldist.stats <- reldist(a.sort, b.sort);

	# print
	cat(" JACCARD a & b: \n"); print(jaccard.stats); cat("\n");
	cat(" RELDIST a & b: \n"); print(reldist.stats); cat("\n");

	# even better way to run both jaccard and reldist, as well as estimate P value through random permutations
	jaccard.reldist.stats <- test.region.similarity(a.sort, b.sort, n = 40);

	cat(" JACCARD/RELDIST a & b: \n"); print(jaccard.reldist.stats); cat("\n");
	}

## ---- results = "hold", message = TRUE, eval = TRUE, errors = TRUE-------
if (check.binary("bedtools")) {

	# read example regions file
	regions.file <- system.file("extdata/example-a-region.bed", package = "bedr");
	a <- read.table(regions.file, header = FALSE, stringsAsFactors = FALSE);
	colnames(a) <- c("a.CHROM", "a.START", "a.END", "Score");

	# sort
	a <- bedr.sort.region(a);

	# group by (on first three columns) the score in column 4. Concatenate scores
	a.collapsed <- bedr(
		input = list(i = a), 
		method = "groupby", 
		params = "-g 1,2,3 -c 4 -o collapse"
		);

	# group by (on first three columns) the score in column 4. Compute mean
	a.mean <- bedr(
		input = list(i = a), 
		method = "groupby", 
		params = "-g 1,2,3 -c 4 -o mean"
		);

	# print
	cat(" REGION a groupby (collapsed): \n"); print(a.collapsed); cat("\n");
	cat(" REGION a groupby (mean): \n"); print(a.mean); cat("\n");
	}

## ---- results = "hold", message = TRUE, eval = TRUE, errors = TRUE-------
if (check.binary("bedtools")) {
	VALID.SV.TYPES <- c('BND', 'CNV', 'DEL', 'DUP', 'INS', 'INV');
	POSITION.COLUMNS <- c('CHROM', 'POS', 'END');

	callerA.filename <- system.file("extdata/callerA.vcf.gz", package = "bedr");
	callerB.filename <- system.file("extdata/callerB.vcf.gz", package = "bedr");

	# read the VCF file
	callerA <- read.vcf(callerA.filename, split.info = TRUE)$vcf;
	callerB <- read.vcf(callerB.filename, split.info = TRUE)$vcf;

	# focus on SVs
	callerA <- callerA[which(callerA$SVTYPE %in% VALID.SV.TYPES), ];
	callerB <- callerB[which(callerB$SVTYPE %in% VALID.SV.TYPES), ];

	# convert to zero-based coordinates
	callerA$POS <- callerA$POS - 1;
	callerB$POS <- callerB$POS - 1;

	# find all overlapping pairs, retrieve size of overlap (bp)
	overlapping.pairs <- bedr.join.region(
		callerA[, POSITION.COLUMNS],
		callerB[, POSITION.COLUMNS],
		report.n.overlap = TRUE,
		check.chr = FALSE
		);
	colnames(overlapping.pairs) <- c(
		'a.CHROM', 'a.POS', 'a.END',
		'b.CHROM', 'b.POS', 'b.END',
		'Overlap'
	    );
	overlapping.pairs$b.POS <- as.numeric(overlapping.pairs$b.POS);
	overlapping.pairs$b.END <- as.numeric(overlapping.pairs$b.END);

	# compute a distance between overlapping pairs
	min.breakpoint.distances <- cbind(
	    overlapping.pairs$a.POS - overlapping.pairs$b.POS,
	    overlapping.pairs$a.END - overlapping.pairs$b.END
	    );

	min.breakpoint.distances <- apply(
	    abs(min.breakpoint.distances),
	    1,
	    min
	    );
	a.length <- overlapping.pairs$a.END - overlapping.pairs$a.POS;
	b.length <- overlapping.pairs$b.END - overlapping.pairs$b.POS;

	overlapping.pairs$distance  <- (min.breakpoint.distances + abs(a.length - b.length)) / 2;

	# print
	cat(" OVERLAPPING PAIRS: \n"); print(head(overlapping.pairs)); cat("\n");
	}

## ---- results = "hold", message = TRUE, eval = TRUE, errors = TRUE-------
if (check.binary("bedtools")) {

	# get Human RefSeq genes (Hg19) in BED format
	refseq.file <- system.file("extdata/ucsc.hg19.RefSeq.chr1-2.txt.gz", package = "bedr");
	refseq <- read.table(refseq.file, header = FALSE, stringsAsFactors = FALSE);

	# sort Refseq and remove chr prefix
	refseq.sorted <- bedr.sort.region(refseq[, 1:4]);
	colnames(refseq.sorted) <- c(POSITION.COLUMNS, "Gene");
	refseq.sorted$CHROM <- gsub("^chr", "", refseq.sorted$CHROM);

	# reuse the SV calls from workflow 1 and add gene identifiers to callerA
	callerA.annotated <- bedr.join.region(
		callerA[, POSITION.COLUMNS],
		refseq.sorted,
		report.n.overlap = TRUE,
		check.chr = FALSE
		);
	colnames(callerA.annotated) <- c(
        'a.CHROM', 'a.POS', 'a.END',
        'Gene.CHROM', 'Gene.POS', 'Gene.END', 'Gene', 'Overlap'
        );

	# reuse the SV calls from workflow 1 and add gene identifiers to callerB
	callerB.annotated <- bedr.join.region(
		callerB[, POSITION.COLUMNS],
		refseq.sorted,
		report.n.overlap = TRUE,
		check.chr = FALSE
		);
	colnames(callerB.annotated) <- c(
        'b.CHROM', 'b.POS', 'b.END',
        'Gene.CHROM', 'Gene.POS', 'Gene.END', 'Gene', 'Overlap'
        );

	# reinstate chr prefix to chromosome names
	callerA.annotated$a.CHROM <- paste('chr', callerA.annotated$a.CHROM, sep = "");
	callerB.annotated$b.CHROM <- paste('chr', callerB.annotated$b.CHROM, sep = "");

	# print
	cat(" CALLER A GENES (chr 1,2): \n"); print(head(callerA.annotated)); cat("\n");
	cat(" CALLER B GENES (chr 1,2): \n"); print(head(callerB.annotated)); cat("\n");
	}

## ---- results = "hold", message = FALSE, warnings = FALSE, eval = TRUE, errors = TRUE----
options("warn" = -1);
if (check.binary("bedtools")) {

	# collapse column 7 (Gene) against unique composite key (column 1, 2 and 3)
	callerA.annotated.grouped <- bedr(
		input = list(i = callerA.annotated), 
		method = "groupby", 
		params = "-g 1,2,3 -c 7 -o collapse"
		);

	# collapse column 7 (Gene) against unique composite key (column 1, 2 and 3)
	callerB.annotated.grouped <- bedr(
		input = list(i = callerB.annotated), 
		method = "groupby", 
		params = "-g 1,2,3 -c 7 -o collapse"
		);

	# print
	cat(" CALLER A GENES (chr 1,2) GROUPED: \n"); print(head(callerA.annotated.grouped)); cat("\n");
	cat(" CALLER B GENES (chr 1,2) GROUPED: \n"); print(head(callerB.annotated.grouped)); cat("\n");
	}

## ---- results = "hide", message = FALSE, eval = TRUE, errors = TRUE, fig.width = 7, fig.height = 7----
if (check.binary("bedtools")) {

	# merge and plot overlapping regions between the two callers with genes 
	# on chromosome 1 and 2
	callerA.merged <- callerA.annotated[, c('a.CHROM', 'a.POS', 'a.END')];
	callerB.merged <- callerB.annotated[, c('b.CHROM', 'b.POS', 'b.END')];
	callerA.merged <- bedr.merge.region(callerA.merged);
	callerB.merged <- bedr.merge.region(callerB.merged);

	# plot (sub-regions exclusive to a and b, and sub-regions in common)
	bedr.plot.region(
		input = list(
			a = callerA.merged, 
			b = callerB.merged
			),
		params = list(lty = 2, label.col = "black", main = "Genes Overlap"),
		feature = 'interval',
		verbose = FALSE
		);
	}

## ---- results = "hide", message = TRUE, eval = TRUE, errors = TRUE-------
if (check.binary("bedtools")) {

	# read copy number segmented (regions) data
	cna.file <- system.file("extdata/CNA.segmented.txt.gz", package = "bedr");
	cna.data <- read.table(cna.file, header = TRUE, stringsAsFactors = FALSE);
	cna.data$Chromosome <- paste("chr", cna.data$Chromosome, sep = "");

	# check if regions are valid
	valid.segments <- is.valid.region(cna.data[, c("Chromosome", "Start", "End")]);
	cat(" VALID REGIONS: ", length(which(valid.segments) ==  TRUE), "/", length(valid.segments), "\n");

	# restrict to copy number regions with |log-ratio| > 0.10
	cna.data <- cna.data[which(abs(cna.data$Segment_Mean) > 0.10), ];

	# create a list data-structure for every patient's copy number data
	# and sort them
	cna.data.gain <- list();
	cna.data.loss <- list();
	pga.gain <- list();
	pga.loss <- list();
	sample.ids <- unique(cna.data$Sample);
	hg19.size <- 3137161264;
	for (sample.id in sample.ids) {

		# extract sample specific gains
		gain.index <- which(cna.data$Sample == sample.id & cna.data$Segment_Mean > 0);
		if (length(gain.index) > 0) {

			cna.data.gain[[sample.id]] <- cna.data[
				gain.index,
				2:4
				];

			# sort
			cna.data.gain[[sample.id]] <- bedr.sort.region(
				x = cna.data.gain[[sample.id]],
				method = "natural",
				verbose = FALSE
				);

			# estimate percent genome gained
			pga.gain[[sample.id]] <- sum(
				apply(
					cna.data.gain[[sample.id]],
					1,
					FUN = function(x) { return( (as.numeric(x[3]) - as.numeric(x[2])) ); }
					)
				);
			pga.gain[[sample.id]] <- pga.gain[[sample.id]]/hg19.size*100;
			}

		# extract sample specific losses
		loss.index <- which(cna.data$Sample == sample.id & cna.data$Segment_Mean < 0);
		if (length(loss.index) > 0) {

			cna.data.loss[[sample.id]] <- cna.data[
				loss.index,
				2:4
				];

			# sort
			cna.data.loss[[sample.id]] <- bedr.sort.region(
				x = cna.data.loss[[sample.id]],
				method = "natural",
				verbose = FALSE
				);

			# estimate percent genome loss
			pga.loss[[sample.id]] <- sum(
				apply(
					cna.data.loss[[sample.id]],
					1,
					FUN = function(x) { return( (as.numeric(x[3]) - as.numeric(x[2])) ); }
					)
				);
			pga.loss[[sample.id]] <- pga.loss[[sample.id]]/hg19.size*100;
			}
		}
	}

## ---- results = "hide", message = FALSE, eval = TRUE, errors = TRUE, fig.width = 7, fig.height = 5----
if (check.binary("bedtools")) {

	# set graphics params
	par(mfrow = c(1, 2), las = 1);

	# plot histograms of Gain and Loss frequencies
	hist(unlist(pga.gain), xlab = "Percent Gain", main = "Percent Genome Altered (Gain)");
	hist(unlist(pga.loss), xlab = "Percent Loss", main = "Percent Genome Altered (Loss)");
	}

## ---- results = "hold", message = TRUE, eval = TRUE, errors = TRUE-------
if (check.binary("bedtools")) {

	# find minimal common regions (gain)
	mcr.gain <- bedr.join.multiple.region(
		x = cna.data.gain,
		species = "human",
		build = "hg19",
		check.valid = FALSE,
		check.sort = FALSE,
		check.merge = FALSE,
		verbose = FALSE
		);

	# find minimal common regions (loss)
	mcr.loss <- bedr.join.multiple.region(
		x = cna.data.loss,
		species = "human",
		build = "hg19",
		check.valid = FALSE,
		check.sort = FALSE,
		check.merge = FALSE,
		verbose = FALSE
		);

	# reorder by frequency of recurrence
	mcr.gain <- mcr.gain[order(as.numeric(mcr.gain$n.overlaps), decreasing = TRUE), ];
	mcr.loss <- mcr.loss[order(as.numeric(mcr.loss$n.overlaps), decreasing = TRUE), ];

	# print
	cat(" RECURRENT CNA GAINS \n"); print(head(mcr.gain[, 1:5])); cat("\n");
	cat(" RECURRENT CNA LOSSES \n"); print(head(mcr.loss[, 1:5])); cat("\n");
	}

