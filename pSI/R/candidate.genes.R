#'@title
#'Candidate Gene Lists
#' 
#'@description
#'This list contains 5 candidate gene lists. They are provided as sample gene lists to be used with \code{fisher.iteration} & \code{candidate.overlap} functions
#'found in this package.
#'\itemize{
#'  \item{\code{candidate.genes$AutDB}}\cr{Autism spectrum disorder candidate gene list from AutDB (N=328)}
#'  \item{\code{candidate.genes$protein.disrupting.rdnv}}\cr{Autism spectrum disorder candidate gene list collected from 4 studies published in 2012 (N = 122)}
#'  \item{\code{candidate.genes$silent.rdnv}}\cr{Autism spectrum disorder negative control gene list collected from 4 studies published in 2012 (N = 122)}
#'  \item{\code{candidate.genes$hcrt.genes}}\cr{Narcolepsy Candidate Gene List (N=63)}
#'  \item{\code{candidate.genes$retinopathy.genes}}\cr{Human Congenital Retinopathies Disease Gene List (N=120)}
#'}
#'NOTE:Supplementary data (human & mouse expression sets, calculated pSI datasets, etc.) can be found in \code{pSI.data} package located at the following URL:
#'\url{http://genetics.wustl.edu/jdlab/psi_package/}

#'
#'@details
#'\itemize{
#'  \item{\code{candidate.genes$AutDB}}\cr{Hand-curated list of Autism Spectrum Disorder (ASD) candidate genes derived from human genetics studies downloaded from AutDB (N=328)}
#'  \item{\code{candidate.genes$protein.disrupting.rdnv}}\cr{List of Protein-Disrupting rare de novo variant affected genes in ASD Probands (N = 122)}
#'  \item{\code{candidate.genes$silent.rdnv}}\cr{List of Silent rare de novo variant affected genes in ASD unaffected siblings (N = 122)}
#'  \item{\code{candidate.genes$hcrt.genes}}\cr{List of differentially dysregulated genes from narcoleptic mice with Hcrt neuron ablation versus control (N=63)}
#'  \item{\code{candidate.genes$retinopathy.genes}}\cr{List of genes identified in human congenital retinopathies downloaded from the curated RetNet database (N=120)}
#'}
#'
#'@examples
#'data(candidate.genes)
#'
#'names(candidate.genes)
#'
#'candidate.genes[[5]]
#'
#'@docType data
#'@keywords datasets
#'@format 5 candidate gene lists, each in the form of a character vector, which are contained within one R list.
#'@name candidate.genes
#'
#'@source
#'\code{AutDB}\cr
#'Basu SN, Kollu R, Banerjee-Basu S (2009): AutDB: a gene 
#'reference resource for autism research. Nucleic Acids Research. 37:D832-D836.
#'\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13379}
#'
#'@source
#'\code{protein.disrupting.rdnv} & \code{silent.rdnv}\cr
#'Iossifov I, Ronemus M, Levy D, Wang Z, Hakker I, Rosenbaum J, et al. (2012): De novo gene disruptions in children on the autistic spectrum. Neuron. 74:285-299.
#'
#'Neale BM, Kou Y, Liu L, Ma'ayan A, Samocha KE, Sabo A, et al. (2012): Patterns and rates of exonic de novo mutations in autism spectrum disorders. Nature. 485:242-245.
#'
#'Sanders SJ, Murtha MT, Gupta AR, Murdoch JD, Raubeson MJ, Willsey AJ, et al. (2012): De novo mutations revealed by whole-exome sequencing are strongly associated with autism. Nature. 485:237-241.
#'
#'O'Roak BJ, Vives L, Girirajan S, Karakoc E, Krumm N, Coe BP, et al. (2012): Sporadic autism exomes reveal a highly interconnected protein network of de novo mutations. Nature. 485:246-250.
#'
#'
#'@source
#'\code{hcrt.genes}\cr
#'Honda M, Eriksson KS, Zhang S, Tanaka S, Lin L, Salehi A, et al. (2009): IGFBP3 colocalizes with 
#'and regulates hypocretin (orexin). PLoS One. 4:e4254.
#'\url{http://www.plosone.org/article/info:doi/10.1371/journal.pone.0004254}
#'
#'@source
#'\code{retinopathy.genes}\cr
#'Daiger, SP. RetNet, the Retinal Information Network.
#'\url{https://sph.uth.edu/RetNet/}
NULL