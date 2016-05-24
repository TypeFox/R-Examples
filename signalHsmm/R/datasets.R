#' @name benchmark_dat
#' @title Benchmark data set for signalHsmm
#' @description Lists eukaryotic proteins added to UniProt database release 
#' 2015_06 between 1.01.2010 and 1.06.2015 (140 proteins with signal peptide and 280 randomly 
#' sampled proteins without signal peptide). 
#' @docType data
#' @usage benchmark_dat
#' @format a list of \code{\link[seqinr]{SeqFastaAA}} objects. 
#' Slot \code{sig} contains the range of signal peptide (if any).
#' @source \href{http://www.uniprot.org/}{UniProt}
#' @examples summary(benchmark_dat)
#' @keywords datasets
NULL

#' @name aaaggregation
#' @title Reduced amino acid alphabet
#' @description Amino acids are grouped together in larger sets based on their 
#' physicochemical properties important in  the recognition of signal peptide.
#' @docType data
#' @usage aaaggregation
#' @format a list of length four. Each element contains a \code{character} vector 
#' of amino acid names (one-letter abbreviations).
#' @keywords datasets
NULL