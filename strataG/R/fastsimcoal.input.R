#' @name fastsimcoal.input
#' @title Input functions for fastsimcoal parameters
#' @description Functions to create \code{pop.info}, \code{locus.params}, and 
#'   \code{hist.ev} input matrices for fastsimcoal function.
#'  
#' @note SNPs are simulated as a diploid DNA sequence with a single nucleotide,
#'   with no recombination (each on its own "chromosome"). If \code{mut.rate} 
#'   has more than one value then this many SNPs are simulated. Otherwise, 
#'   \code{num.loci} independent SNPs are simulated.
#'   
#' @references Excoffier, L. and Foll, M (2011) fastsimcoal: a continuous-time 
#'   coalescent simulator of genomic diversity under arbitrarily complex 
#'   evolutionary scenarios Bioinformatics 27: 1332-1334.\cr
#'   \url{http://cmpg.unibe.ch/software/fastsimcoal2/}
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{fastsimcoal}}
#' 
NULL


#' @rdname fastsimcoal.input
#' 
#' @param pop.size a vector giving size of each populaiton.
#' @param sample.size a vector giving the number of samples to take from each
#'   population.
#' @param sample.times a vector giving the number of generations in the past
#'   at which samples are taken.
#' @param growth.rate a vector giving the growth rate of each population.
#' 
#' @export
#' 
fscPopInfo <- function(pop.size, sample.size, sample.times = 0, growth.rate = 0) {
  cbind(
    pop.size = pop.size, sample.size = sample.size,
    sample.times = sample.times, growth.rate = growth.rate
  )
}


#' @rdname fastsimcoal.input

#' @param locus.type a character representation of what type of marker to simulate.
#'   Can be "dna", "msat", or "snp".
#' @param sequence.length \code{dna}: number of DNA base pairs to use.
#' @param num.loci \code{msat, snp}: number of loci to simulate.
#' @param mut.rate \code{dna, msat}: per base pair or locus mutation rate.
#' @param transition.rate dna: fraction of substitutions that are transitions. 
#'   Set to 1 (all transitions) for SNPs.
#' @param gsm.param \code{msat}: Value of the geometric parameter for a
#'   Generalized Stepwise Mutation (GSM) model. This value represents the
#'   proportion of mutations that will change the allele size by more than
#'   one step. Values between 0 and 1 are required. A value of 0 is for a
#'   strict Stepwise Mutation Model (SMM).
#' @param range.constraint \code{msat}: Range constraint (number of different
#'   alleles allowed). A value of 0 means no range constraint.
#' @param recomb.rate recombination rate between adjacent markers. No effect for 
#'   SNPs.
#' @param chromosome number or character identifying which chromosome the marker
#'   is on.
#' @param num.chrom a value giving the number of chromosomes that the
#'   \code{locus.params} marker specifications should be copied for. If
#'   \code{NULL}, then chromosome assignment is taken from the \code{chromosome} 
#'   column. Any non-\code{NULL} integer will cause the value in \code{chromosome}.
#' @param ploidy positive integer giving the ploidy of the marker type to 
#'   be simulated.
#'   
#' @export
#' 
fscLocusParams <- function(locus.type = c("dna", "msat", "snp"),
                           sequence.length = NULL, num.loci = NULL, 
                           mut.rate = NULL, transition.rate = 1 / 3, 
                           gsm.param = 0, range.constraint = 0, recomb.rate = 0, 
                           chromosome = NULL, num.chrom = NULL, ploidy = 2) {
  
  createLocusParams <- function(chr, type, num.markers, recomb.rate, param.4,
                                param.5, param.6, ploidy, num.chrom) {
    if(is.null(chr)) chr <- 1
    df <- data.frame(
      chromosome = chr, type = type, num.markers = num.markers,
      recomb.rate = recomb.rate, param.4 = param.4, param.5 = param.5,
      param.6 = param.6, stringsAsFactors = FALSE
    )
    df <- df[order(df$chromosome), ]
    attr(df, "num.chrom") <- num.chrom[1]
    attr(df, "ploidy") <- ploidy
    return(df)
  }
  
  df <- switch(
    match.arg(locus.type),
    dna = createLocusParams(
      chromosome, "DNA", sequence.length, recomb.rate, mut.rate,
      transition.rate, NA, 1, num.chrom
    ),
    msat = createLocusParams(
      chromosome, "MICROSAT", num.loci, recomb.rate, mut.rate, gsm.param,
      range.constraint, 2, num.chrom
    ),
    snp = if(length(mut.rate) == 1) {
      createLocusParams(
        1, "DNA", 1, 0, mut.rate, 1, NA, 2, num.loci
      )
    } else {
      createLocusParams(
        1:length(mut.rate), "DNA", 1, 0, mut.rate, 1, NA, 2, NULL
      )
    }
  )
  attr(df, "opts") <- if(locus.type == "snp") "-s" else ""
  return(df)
}


#' @rdname fastsimcoal.input
#' 
#' @param num.gen Number of generations, t, before present at which the
#'   historical event happened.
#' @param source.deme Source deme (the first listed deme has index 0)
#' @param sink.deme Sink deme
#' @param prop.migrants Expected proportion of migrants to move from source to sink.
#' @param new.sink.size New size for the sink deme, relative to its size at
#'   generation t.
#' @param new.sink.growth New growth rate for the sink deme.
#' @param new.mig.mat New migration matrix to be used further back in time.
#' 
#' @export
#' 
fscHistEv <- function(num.gen = 0, source.deme = 0, sink.deme = 0,
                      prop.migrants = 1, new.sink.size = 1,
                      new.sink.growth = 0, new.mig.mat = 0) {
  cbind(
    num.gen = num.gen, source.deme = source.deme, sink.deme = sink.deme,
    prop.migrants = prop.migrants, new.sink.size = new.sink.size,
    new.sink.growth = new.sink.growth, new.mig.mat = new.mig.mat
  )
}