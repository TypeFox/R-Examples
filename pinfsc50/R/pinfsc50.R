#' pinfsc50: A package containing the sequence, annotations and variants for \emph{P. infestans}.
#'
#' The pinfsc50 package contains data from \emph{Phytophthora infestans} intended to be used as example data.
#' 
#' @section pinfsc50 functions:
#' This package contains no functions.
#' 
#' @section Files in inst/extdata:
#' 
#' \strong{pinf_sc50.fasta} - fasta format file containing the nucleotide sequence for \emph{Phytophthora infestans} T30-4 supercontig_1.50.
#' This sequence is a subset from a file downloaded from the \href{http://www.broadinstitute.org/annotation/genome/phytophthora_infestans/MultiHome.html}{Broad's \emph{P. infestans}} page.
#' This data was published in Haas et al. (2009).
#' 
#' \strong{pinf_sc50.gff} - gff format file containing annotations for supercontig_1.50
#' This file is a subset from a file downloaded from the \href{http://www.broadinstitute.org/annotation/genome/phytophthora_infestans/MultiHome.html}{Broad's \emph{P. infestans}} page called .
#' 
#' \strong{pinf_sc50.vcf.gz} - gzipped vcf format file containing variant information for supercontig_1.50.
#' This file was created with the GATK's haplotype caller.
#' The data were then phased with beagle4.
#' beagle4 returns a vcf file with which lacks much of the diagnostic information contained in the input file.
#' I therefore stripped the unphased genotypes from the original file and pasted on the phased genotypes from beagle4 to create a vcf file with the genotypes from beagle4, but all of the other information contained in the GATK's haplotype caller's file.
#'
#'
#' Short read data for sample t30-4 was downloaded from the \href{http://www.broadinstitute.org/annotation/genome/phytophthora_infestans/MultiHome.html}{Broad's \emph{P. infestans}} page.
#' Short read data from sample blue13 was published in Cooke et al. (2012).
#' Short read data from samples DDR7602, LBUS5, NL07434, P10127, P10650, P11633, P12204, P13527, P1362, P13626, P17777us22, P6096 and P7722 were published in Yoshida et al. (2013).
#' Short read data from samples BL2009P4_us23, IN2009T1_us22, RS2009P1_us8 were published in Martin et al (2013).
#'
#'
#'
#'
#'
#' @references
#' 
#' Cooke, D. E., Cano, L. M., Raffaele, S., Bain, R. A., Cooke, L. R., Etherington, G. J., ... & Kamoun, S. (2012). Genome analyses of an aggressive and invasive lineage of the Irish potato famine pathogen. PLoS Pathog, 8(10), e1002940.
#' 
#' Haas, B. J., Kamoun, S., Zody, M. C., Jiang, R. H., Handsaker, R. E., Cano, L. M., ... & Liu, Z. (2009). Genome sequence and analysis of the Irish potato famine pathogen \emph{Phytophthora infestans}. Nature, 461(7262), 393-398.
#' 
#' Martin, M. D., Cappellini, E., Samaniego, J. A., Zepeda, M. L., Campos, P. F., Seguin-Orlando, A., ... & Gilbert, M. T. P. (2013). Reconstructing genome evolution in historic samples of the Irish potato famine pathogen. Nature communications, 4.
#' 
#' \emph{Phytophthora infestans} Sequencing Project, Broad Institute of Harvard and MIT (\url{http://www.broadinstitute.org/}). 
#'
#' Yoshida, K., Schuenemann, V. J., Cano, L. M., Pais, M., Mishra, B., Sharma, R., ... & Burbano, H. A. (2013). The rise and fall of the \emph{Phytophthora infestans} lineage that triggered the Irish potato famine. Elife, 2, e00731.
#'
#'
#'
#' @examples 
#' 
#' \dontrun{
#' dna <- system.file("extdata", "pinf_sc50.fasta", package = "pinfsc50")
#' dna <- ape::read.dna(dna, format="fasta")
#' gff <- system.file("extdata", "pinf_sc50.gff", package = "pinfsc50")
#' gff <- read.table(gff, header=FALSE, sep="\t", quote = "")
#' vcf <- system.file("extdata", "pinf_sc50.vcf.gz", package = "pinfsc50")
#' vcf <- vcfR::read.vcf(vcf)
#' }
#'
#' @docType package
#' @name pinfsc50
NULL
#> NULL


