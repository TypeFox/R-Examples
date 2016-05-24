#'@title Variation in Cell Abundance Loci
#'@name vocal
#'@aliases vocal
#'
#'@description Probing immune system genetics via gene expression.
#' VoCAL is a deconvolution-based method that utilizes transcriptome data
#' to infer the quantities of immune-cell types, and then uses these quantitative
#' traits to uncover the underlying DNA loci (iQTLs) assuming homozygosity
#' (such as in the case of recombinent inbred strains).
#'
#'@param reference_data a data frame representing immune cell expression profiles.
#'  Each row represents an expression of a gene, and each column represents a
#'  different immune cell type. \code{colnames} contains the name of each immune cell
#'  type and the \code{rownames} includes the genes' symbol. The names of each immune
#'  cell type and the symbol of each gene should be unique. Any gene with
#'  missing expression values must be excluded.
#'
#'@param expression_data a data frame representing RNA-seq or microarray
#'  gene-expression profiles of a given complex tissue across a population of
#'  genetically distinct (genotyped) individuals. Each row represents an
#'  expression of a gene, and each column represents a genetically distinct
#'  individual. \code{colnames} contain the name of each individual, as written in the
#'  \code{genotyping_data}, and \code{rownames} includes the genes' symbol.
#'  The name of each individual sample and the symbol of each gene should be unique.
#'  Any gene with missing expression values should be excluded.
#'
#'@param genotyping_data a data frame where each row represents a different
#'  locus, and each column represents a genetically distinct individual.
#'  The genotype should be taken from homozygous individuals only.
#'  Where the genotype is unknown \code{NA} should be used.
#'  The first six columns contain the following information: (1) The sequential
#'  identifier of the locus; (2) The name of each locus Chr; (3) Chromosome
#'  position; (4) Start genome position; (5) End genome position;
#'  (6) position in cM.
#'
#'@param normalize_data normalization type. The data will be normalized by either:
#'  (1) "All" - subtraction of the mean expression of all strains;
#'  (2) "None" - data is already normalized, do nothing;
#'  (3) name of individual included in \code{colnames} of \code{expression_data};
#'
#'@param T.i numerical. significant iQTL association score \code{(-log10(Pvalue))}
#'  cutoff for the refinement step of the VoCAL algorithm.
#'
#'@param T.e numerical. significant eQTL association score \code{(-log10(Pvalue))}
#'  cutoff for the refinement step of the VoCAL algorithm.
#'
#'@param eqtl_association_scores (optional) a data frame where each entry
#'  represents an association score for a gene given the genotype of all the
#'  individuals that appear in the expression_data data frame, in a specific locus.
#'  This eQTL analysis should be peformed over the normalized expression_data.
#'  \code{colnames} contain the UID (as written in the genotyping_data) and
#'  \code{rownames} includes the genes' symbol (as written in the expression_data).
#'  The symbol of each gene should be unique. These scores should be in -log10(P value).
#'  Default is NULL, meaning that eQTL analysis will be performed.
#'
#'@param ... one or more data frames of one column, each one represents a
#'  preselected marker set that likely discriminate well between the immune-cell
#'  types given in the reference data. The number of data frames defines the
#'  number of association scores that would be combined to generate the final
#'  iQTL association score.
#'
#'@return a list of two martices
#'  \item{final_association_score}{a matrix that contains the output iQTL association
#'  score after applying the iterative filteration procedure. Each row represents the genome
#'  wide-association result for a specific immune trait over a range of DNA loci.
#'  \code{rownames} provides the identifier of the locus and \code{colnames} contains the
#'  immune-cell type names. Each entry provides the \code{-log10(P value)} of an iQTL
#'  association score.}
#'  \item{marker_info}{the names of all the markers removed from the different
#'  marker sets provided}
#'
#'
#'@usage vocal(...,reference_data,expression_data,genotyping_data,normalize_data,
#'  T.i=5,T.e=10,eqtl_association_scores=NULL)
#'
#'@references Steuerman Y and Gat-Viks I.
#'  Exploiting Gene-Expression Deconvolution to Probe the Genetics of the Immune System (2015), Submitted.
#'
#'@examples
#'data(commons)
#'data(vocalEx)
#'\dontrun{
#'  results <- vocal(DCQ_mar, reference_data=immgen_dat, expression_data=lung_dat,
#'  genotyping_data=gBXD, normalize_data="B6", eqtl_association_scores=eQTL_res)
#'}
#'
#' @export
vocal <- function(..., reference_data,
                  expression_data,
                  genotyping_data,
                  normalize_data,
                  T.i=5, T.e=10,
                  eqtl_association_scores=NULL) {

  models <- list(...)
  ### print number of marker sets
  num_of_marker_sets <- length(models)
  print(paste("marker sets:", num_of_marker_sets, sep=" "))

  marker_sets <- c()
  for (i in seq(1:num_of_marker_sets)) {
    marker_sets <- c(marker_sets, models[[i]][,1])
  }

  normed_data <- .normalizeData(orig_data = expression_data, norm_method = normalize_data)

  if (is.null(eqtl_association_scores)) {
    eqtl_association_scores <- .qtlAnalysis(exp_data = normed_data, genes_to_use = marker_sets,
                                            genotype_data = genotyping_data, transpose.file=T)$logPval
  }

  list_of_results_mat <- list()

  for (i in seq(1:num_of_marker_sets)) {

    marker_set_name <- colnames(models[[i]][1])
    print(marker_set_name)

    list_of_results_mat[[marker_set_name]] <- .create_iQTL_association_scores(reference_data = reference_data,
                                                                           mix_data = normed_data,
                                                                           marker_set = models[[i]][1],
                                                                           genotype_data = genotyping_data)
  }

  fisher_mat_log <- .combine_iQTL_association_scores(combined_results = list_of_results_mat)

  finalRes <- .iterative_filteration(models,
                                     fisher_combined_res=fisher_mat_log,
                                     res_by_marker_set = list_of_results_mat,
                                     eQTL_marker_set = eqtl_association_scores,
                                     T.i = T.i, T.e = T.e,
                                     reference_data = reference_data,
                                     expression_data = normed_data,
                                     genotyping_data = genotyping_data)

  return(finalRes)
}
