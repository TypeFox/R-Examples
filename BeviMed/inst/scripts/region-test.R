get_sample_names <- function(file, connection_type=base::file) {
	cnx <- connection_type(file, open="r")
	on.exit(close(cnx))
	headings <- NULL
	while (is.null(headings)) {
		l <- readLines(cnx, n=1)
		if (grepl(x=l, pattern="^#CHROM")) headings <- strsplit(l, split="\t")[[1]]
	}
	headings[10:length(headings)]
}

compressed_vcf_sample_names <- function(vcf_file_name) {
	get_sample_names(paste0("zcat ", vcf_file_name), connection_type=pipe)
}

just_counts <- function(parts, description_columns=9) {
	y <- sapply(parts, "[", -(1:description_columns))
	structure(grepl(x=y, pattern="^[^0.][/|].") + grepl(x=y, pattern="^.[/|][^0.]"), dim=c(if (length(parts) > 0) length(parts[[1]])-description_columns else 0, length(parts)))
}

var_info <- function(parts, description_columns=9) structure(dimnames=list(NULL, c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","FORMAT","INFO")), structure(dim=c(length(parts),description_columns), t(sapply(parts, "[", 1:description_columns))))

test_region <- function(vcf_file_name, chr, from, to, case_IDs, samples=compressed_vcf_sample_names(vcf_file_name), description_columns=9, min_ac=1, ...) {
	file_samples <- compressed_vcf_sample_names(vcf_file_name) 
	sample_inds <- match(samples, file_samples)
	cmd <- paste("tabix ", vcf_file_name, " ", chr, ":", from, "-", to, sep="")
	z <- pipe(cmd)
	lines <- grep(value=TRUE, pattern="^#", invert=TRUE, x=readLines(z))
	close(z)
	parts <- strsplit(lines, split="\t")
	counts <- t(just_counts(parts, description_columns))
	info <- var_info(parts, description_columns)
	x <- summary(bevimed(
		y=samples %in% case_IDs, 
		G=counts[,sample_inds,drop=FALSE], 
		min_ac=min_ac,
		...
	))
	x[["Z"]] <- setNames(x[["Z"]], paste0(info[,1],":",info[,2]))
	x
}
