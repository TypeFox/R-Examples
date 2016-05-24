#' Estimate allelic frequencies in reference database
#' 
#' @param x profiles object
#' @param labels (optional) list of per-locus labels of alleles (repeat numbers for integers, like factor levels)
#' @details The allele frequencies are estimated using the counting method. That is, the empirical fraction of each allele is taken as an estimate of the frequency.
#'            
#'          Since alleles are stored as integer, labels can be supplied that map the integer to a repeat number (similar to a factor level). See below for an example.
#' @examples
#' data(freqsNLsgmplus)
#' 
#' set.seed(123)
#' 
#' # sample a small reference db
#' x <- sample.profiles(N = 1e3,freqs=freqsNLsgmplus)
#' 
#' # estimate frequencies
#' fr0 <- est.freqs(x,labels = lapply(get.freqs(x),names))
#'
#' # mean absolut difference between fr0 and freqsNLsgmplus is small
#' mean(abs(c(fr0,recursive = TRUE)-c(freqsNLsgmplus,recursive=TRUE)))
#' @export
est.freqs <- function(x,labels){
  x <- Zassure.matrix(x)
  x.loci <- Znames.to.loci(colnames(x))
  Amax <- sapply(x.loci, function(L) max(x[,paste(L,1:2,sep = ".")],na.rm = TRUE))  

  if (missing(labels)){
    labels <- lapply(Amax, function(A) as.character(seq_len(A)))
  }else{
    # check the labels
    if (any(is.na(match(x.loci,names(labels))))){
      stop("Please supply labels for all loci.")
    }
    if (any(sapply(labels[x.loci],length)<Amax)){
      stop("Higher allele numbers in x than labels for at least one locus.")
    }
  }
  
  Lmax <- sapply(labels[x.loci],length)
  counts <- Zcountallelescpp(x = x, Amax = max(Lmax),w = rep(1,nrow(x)))
  
  ret <- list()
  for( l in seq(x.loci)){
    c0 <- counts[1:Lmax[l],l]
    ret[[x.loci[l]]] <- setNames(c0/sum(c0),labels[[x.loci[l]]])
  }
  
  ret
}
