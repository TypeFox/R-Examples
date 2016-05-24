#' Calculate the 'enrichment score' for gene-set amongst ranking/correlations
#'
#' @param N Integer value
#' @param S Ranks of gene set
#' @param p Weighting of ranking/correlations
#' @param r Rank/correlation scores
#' @return Numeric value - enrichment score
es <- function(N, S, p=1, r=N:1/N) {
	stopifnot(N == length(r))
	stopifnot(all(S) %in% names(r) | is.integer(S))
	.Call(
		"R_es",
		sort((if (is.character(S)) match(S, names(r)) else S)-1),
		p,
		r,
		PACKAGE="gsEasy"
	)
}

#' Gene set enrichment test
#'
#' @param S Ranks of gene set
#' @param N Integer value. Only required if \code{r} is not specified.
#' @param r Rank/correlation scores. If \code{S} is \code{character}, then \code{r} must be named by gene or be a character vector itself of the gene names (necessarily containing \code{S}) in rank order.
#' @param p Weighting of ranking/correlations, see Subramanian et. al 2005.
#' @param min_its Minimum number of null permutations to compare.
#' @param max_its Maximum number of null permutations to compare.
#' @param significance_threshold Maximum p-value of significant result.
#' @param log_dismiss Threshold log probability of returning a significant result, below which function returns current p-value.
#' @return Numeric value - p-value of enrichment.
#' @examples 
#' gset(S=1:5 * 2, N=1000)
#' gset(S=letters[1:3], r=letters)
#' @export
#' @importFrom Rcpp evalCpp
#' @importFrom stats setNames
#' @useDynLib gsEasy
gset <- function(
	S,
	N=NULL,
	r=NULL,
	p=1,
	min_its=2e2,
	max_its=1e5,
	significance_threshold=0.05,
	log_dismiss=-10
) {
	if (is.null(N) & is.null(r))
		stop("Must specify either N or r!")

	if (is.null(N))
		N <- length(r)

	if (is.null(r))
		r <- (N:1)/N

	r_sorted <- if (is.character(r)) setNames(nm=r, (length(r):1/length(r))) else sort(decreasing=TRUE, r)
	stopifnot(is.vector(S) & (all(S %in% names(r_sorted)) | is.numeric(S)))
	.Call(
		"R_gset",
		N,
		sort((if (is.character(S)) match(S, names(r_sorted)) else as.integer(S))-1),
		p,
		r_sorted,
		min_its,
		max_its,
		significance_threshold,
		log_dismiss
	)
}

#' Create list of gene sets defined by ontological annotation 
#'
#' @param ontology \code{ontology_index} object.
#' @param gene Character vector of genes.
#' @param term Character vector of term IDs annotated to corresponding genes.
#' @param min_genes Minimum number of genes in gene sets.
#' @param max_genes Maximum number of genes in gene sets.
#' @return List of character vectors of term IDs.
#' @export
#' @importFrom ontologyIndex get_ancestors
get_ontological_gene_sets <- function(
	ontology,
	gene,
	term,
	min_genes=1,
	max_genes=500
) {
	gene.anno <- lapply(split(term, gene), get_ancestors, ontology=ontology)
	genes.by.term <- lapply(FUN=as.character, X=split(unlist(mapply(SIMPLIFY=FALSE, FUN=rep, names(gene.anno), sapply(gene.anno, length))), unlist(gene.anno)))
	Filter(x=genes.by.term, f=function(x) length(x) <= max_genes & length(x) >= min_genes)
}

#' Create list of gene sets defined by GO term annotation
#'
#' Note, this function takes several minutes to execute.
#'
#' @param GO_annotation_file File path of annotation file, which should contain a column of genes and a column of terms. Can be downloaded from at http://geneontology.org/gene-associations/gene_association.goa_human.gz.
#' @param GO_file File path of gene ontology.
#' @param min_genes Minimum number of genes in gene sets.
#' @param max_genes Maximum number of genes in gene sets.
#' @param verbose Print progress.
#' @return List of character vectors of term IDs.
#' @export
#' @importFrom ontologyIndex get_ontology
#' @importFrom utils read.table
get_GO_gene_sets <- function(
	GO_annotation_file,
	GO_file="http://purl.obolibrary.org/obo/go.obo",
	min_genes=15,
	max_genes=500,
	verbose=TRUE
) {
	if (verbose) cat("reading ontology file...\n")
	go <- get_ontology(GO_file, qualifier="GO")

	if (verbose) cat("reading annotation file...\n")
	anno.df <- read.table(GO_annotation_file, sep="\t", quote="", stringsAsFactors=FALSE, comment.char="!", skipNul=TRUE)

	if (verbose) cat("creating gene set list...\n")
	get_ontological_gene_sets(ontology=get_ontology(GO_file, "GO"), term=anno.df[,5], gene=anno.df[,3], min_genes=min_genes, max_genes=max_genes)
}
