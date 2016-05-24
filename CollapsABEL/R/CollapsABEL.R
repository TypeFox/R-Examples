#' CollapsABEL: an R library for detecting compound heterozygote alleles in genome-wide association or sequencing studies
#' 
#' Compound Heterozygosity (CH) in classical genetics is the presen-ce of two
#' different recessive mutations at a particular gene locus, one on each chromosome
#' . The presence of CH has been found for nearly all
#' autosomal recessive disorders as well as other phenotypes such as red hair
#' color. A relaxed form of CH, i.e., in which the genetic variants are not
#' necessarily coding, rare, and deleterious, is likely involved in a wide range of
#' human polygenic traits and referred to as generalized CH (GCH). Howev-er,
#' individually analyzing a large number of DNA sequence vari-ants, as being the
#' routine in genome-wide association studies (GWAS), has limited power to detect
#' genetic associations caused by GCH, which may be partially responsible for the
#' currently still "missing heritability".
#' Existing tools specifically designed for detecting GCH alleles are scarce, in
#' particular for the analysis of densely imputed Single Nucleotide Polymorphism
#' (SNP) array data or whole genome se-quencing data. Previously, we developed a
#' collapsed double heter-ozygosity (CDH) test for detecting the association
#' between CH genotypes and binary traits by applying a chi-squared statistic to
#' pseudo-genotypes collapsed from a pair of SNPs, which was
#' implemented as a function in the GenABEL R package .
#' Here, we implement a generalized CDH (GCDH) method to overcome previous
#' limitations and allow (1) fast analysis of densely imputed SNP data or whole
#' genome se-quencing data; (2) flexible analysis of binary and quantitative traits
#' with covariates; (3) empirical power estimation and type-I error control; and
#' (4) easy interface with graphical utilities

#' @docType package
#' @name CollapsABEL
#' @aliases collapsabel collapsabel-package CollapsABEL-package
#' @param phe_file character. Phenotype file.
#' @return FALSE when the file is invalid, or a data.frame when it is. 
#' @author Kaiyin Zhong, Fan Liu
NULL