#' Compute allele specific homozygosity.
#'
#' A function to compute allele specific homozygosity values for haplotype frequency data. 
#' The allele specific homozygosity is the homozygosity statistic computed for alleles at 
#' one locus that are found on haplotypes with a specific allele (the focal allele) at the 
#' other locus (the focal locus). 
#' 
#' @param dat	A data.frame with 5 required variables (having the names listed below):
#'      \tabular{ll}{
#'           \code{haplo.freq} \tab A numeric vector of haplotype frequencies.\cr
#'            \code{locus1} \tab A character vector indentifying the first locus.\cr
#'            \code{locus2} \tab A character vector indentifying the second locus.\cr
#'            \code{allele1} \tab A character vector indentifying the allele at locus 1.\cr
#'            \code{allele2} \tab A character vector indentifying the allele at locus 1.\cr
#'            }
#' @param tolerance	A threshold for the sum of the haplotype frequencies.
#'   If the sum of the haplotype frequencies is greater than 1+tolerance or less
#'   than 1-tolerance an error is returned. The default is 0.01.
#' @param sort.var	a vector of variable names specifying the "sort by" variables. 
#'                      The default is c("focal","allele").
#' @param sort.asc	a vector of TRUE/FALSE values, with the same length as "sort.var", 
#'                      indicating whether sorting of each variable is in ascending order. 
#'                      The default order is ascending.
#'
#' @return The return value is a dataframe with the following components:
#'  \tabular{ll}{
#'  \code{loci}	\tab The locus names separated by "-".\cr
#'  \code{focal}	\tab The name of the focal locus (locus conditioned on).\cr
#'  \code{allele}	\tab The name of the focal allele (allele conditioned on).\cr
#'  \code{allele.freq}	\tab The frequency of the focal allele.\cr
#'  \code{as.homz}	\tab The allele specific homozygosity (on haplotypes with the focal allele).\cr
#'  }
#'
#' @section Details:
#' A warning message is given if the sum of the haplotype frequencies is greater than 1.01 or less
#' than 0.99 (regardless of the \code{tolerance} setting). The haplotype frequencies that are 
#' passed to the function are normalized within the function to sum to 1.0 by dividing each 
#' frequency by the sum of the passed frequencies.
#'
#' @examples
#' library(asymLD)
#' 
#' # An example using haplotype frequencies from Wilson(2010)
#' data(hla.freqs)
#' hla.dr_dq <- hla.freqs[hla.freqs$locus1=="DRB1" & hla.freqs$locus2=="DQB1",]
#' compute.ALD(hla.dr_dq)
#' compute.AShomz(hla.dr_dq, sort.var=c("focal","allele"), sort.asc=c(TRUE,TRUE))
#' compute.AShomz(hla.dr_dq, sort.var=c("focal","allele.freq"), sort.asc=c(FALSE,FALSE))
#' # Note that there is substantially less variablity (higher ALD) for HLA*DQB1 
#' # conditional on HLA*DRB1 than for HLA*DRB1 conditional on HLA*DQB1, indicating 
#' # that the overall variation for DQB1 is relatively low given specific DRB1 alleles.
#' # The largest contributors to ALD{DQB1|DRB1} are the DRB1*0301 and DRB1*1501 focal 
#' # alleles, which have high allele frequencies and also have high allele specific 
#' # homozygosity values. 
#' @export

compute.AShomz <- function(dat, tolerance=.01, sort.var=c("focal","allele"), 
  sort.asc=rep(TRUE, length(sort.var))) {
  names.req <- c("locus1", "locus2", "allele1", "allele2", "haplo.freq")
  check.names <- names.req %in% names(dat)
  if (sum(!check.names) > 0)
    stop("The following required variables are missing from the dataset: ",
      paste(names.req[!check.names], collapse = " "))
  if (sum(dat$haplo.freq) > 1 + tolerance)
    stop("Sum of haplo freqs > 1; sum=", sum(dat$haplo.freq))
  if (sum(dat$haplo.freq) < 1 - tolerance)
    stop("Sum of haplo freqs < 1; sum=", sum(dat$haplo.freq))
  if (abs(1 - sum(dat$haplo.freq)) > 0.01)
    warning("Sum of haplo freqs is not 1; sum=", sum(dat$haplo.freq))
  locus1 <- unique(dat$locus1)
  locus2 <- unique(dat$locus2)
  
  # normalize haplotype freqs to sum to 1.0
  dat$haplo.freq <- dat$haplo.freq/sum(dat$haplo.freq)

 #generate allele freqs based on haplo freqs
  by.vars1 <- list(dat$allele1) 
  by.vars2 <- list(dat$allele2)
    names(by.vars1) <- c("allele1")
    names(by.vars2) <- c("allele2")
  af1 <- aggregate(dat$haplo.freq, by=by.vars1, FUN=sum)
  af2 <- aggregate(dat$haplo.freq, by=by.vars2, FUN=sum)
    names(af1)[length(names(af1))] <- "allele.freq1"
    names(af2)[length(names(af2))] <- "allele.freq2"
  mrg1 <- merge(dat,  af1, by.x=c("allele1"), by.y=c("allele1"), all.x=TRUE, all.y=FALSE) 
  mrg2 <- merge(mrg1, af2, by.x=c("allele2"), by.y=c("allele2"), all.x=TRUE, all.y=FALSE) 
  dat <- mrg2  
  
  locus.name <- paste(locus1, locus2, sep="-") 

  rbind.as.homz <- NULL
  for (a1 in unique(dat$allele1)) {
    af.1 <- unique(dat$allele.freq1[dat$allele1 == a1])
    den <- sum(dat$haplo.freq[dat$allele1 == a1])
    rel.freq <- (dat$haplo.freq[dat$allele1 == a1]) / den
    as1.homz <- sum(rel.freq^2)
    as1.homz.vec <- c(locus.name, locus1, a1, af.1, as1.homz)
    names(as1.homz.vec) <- c("loci", "focal", "allele", "allele.freq", "as.homz")
    rbind.as.homz <- rbind(rbind.as.homz, as1.homz.vec)
  }
  for (a2 in unique(dat$allele2)) {
    af.2 <- unique(dat$allele.freq2[dat$allele2 == a2])
    den <- sum(dat$haplo.freq[dat$allele2 == a2])
    rel.freq <- (dat$haplo.freq[dat$allele2 == a2]) / den
    as2.homz <- sum(rel.freq^2)
    as2.homz.vec <- c(locus.name, locus2, a2, af.2, as2.homz)
    names(as1.homz.vec) <- c("loci", "focal", "allele", "allele.freq", "as.homz")
    rbind.as.homz <- rbind(rbind.as.homz, as2.homz.vec)
  }   
  rownames(rbind.as.homz) <- NULL
  as.homz.dat <- as.data.frame(rbind.as.homz)
  for(col.no in 4:5) as.homz.dat[,col.no] <- round(as.numeric(as.character(as.homz.dat[,col.no])),10)
  ASF <- lsort(as.homz.dat, by=sort.var, asc=sort.asc)
  return(ASF)  
}
 