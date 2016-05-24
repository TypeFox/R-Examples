gquad_overlap_base <- function(a){
  a4 <- "(?=(G{2,7}?[A|C|G|T|U|N]{1,36}?G{2,7}?[A|C|G|T|U|N]{1,36}?G{2,7}?[A|C|G|T|U|N]{1,36}?G{2,7}?))"
  a5 <- gregexpr(a4, a, ignore.case = TRUE, perl = TRUE)
  sequence_position <- a5[[1]][1:length(a5[[1]])]
  sequence_length <- as.vector(attr(a5[[1]], "capture.length"))
  a_end = sequence_position + sequence_length - 1
  sequence <- substring(a, sequence_position, a_end)
  a10 <- cbind(sequence_position, sequence, sequence_length)
  a11 <- a10[as.integer(as.character(a10[,3]))<46,]

  if(length(a11) > 3){
    a12 <- as.data.frame(a11, stringsAsFactors=FALSE)
    a12 [,1] <- as.numeric(as.character(a12 [,1]))
    a12 [,3] <- as.integer(as.character(a12 [,3]))
    a13 <- a12[order(a12[,1], a12[,3]),]
    return(a13)
  }

  if((length(a11) == 3)){
    if(as.integer(as.character(a11[[3]])) < 3 ){
      a14 <- data.frame("sequence_position" = "-", "sequence" = "-", "sequence_length" = "-")
      return(a14)
    }else{
      a14 <- data.frame("sequence_position" = a11[[1]], "sequence" = a11[[2]], "sequence_length" = a11[[3]])
      return(a14)
    }
  }

  if(length(a11) < 3){
    a14 <- data.frame("sequence_position" = "-", "sequence" = "-", "sequence_length" = "-")
    return(a14)
  }
}


gquadOverlap_main <- function(b){
  if(length(b) == 1){
    b1 <- gquad_overlap_base(b)
    return(b1)
  }else{
    input_pos = 0
    q <- data.frame("input_ID" = integer(0), "sequence_position" = character(0), "sequence" = character(0), "sequence_length" = character(0))
    for(i in b){
      b1 <- gquad_overlap_base(i)
      input_pos = input_pos + 1
      b2 <- cbind(input_ID = input_pos, b1)
      b2[,c(2,4)] <- sapply(b2[,c(2,4)],as.character)
      q <- rbind(q, b2)
    }
    return(q)
  }
}


#' Predicting G quadruplex motif(s) including overlaps
#'
#' This function predicts G quadruplex motif(s)
#' in 'x' (nucleotide sequence(s)) like the gquad function, but includes overlaps.
#' Nucleotide sequence can be provided in raw or fasta format or as GenBank accession number(s).
#' Internet is needed to connect to GenBank database, if accession number(s) is given as argument.
#'
#' @param x nucleotide sequence(s) in raw format or a fasta file or a GenBank accession number(s); from which G quadruplex motif(s) (including overlaps) will be predicted.
#'  If the fasta file name does not contain an absolute path, the file name is relative to the current working directory.
#' @param xformat a character string specifying the format of x : default (raw), fasta, GenBank accession number(s).
#' @return A dataframe of G quadruplex motif(s) position, sequence and length. If more than one nucleotide sequence is provided as argument, an input ID is returned for motif(s) predicted from each input sequence.
#' @author Hannah O. Ajoge
#' @details
#' This function predicts G quadruplex motif(s) in nucleic (both DNA and RNA) sequences, including overlaps and provide the position, sequence and length of the predicted motif(s), if any.
#' @export
#' @importFrom ape read.GenBank
#' @importFrom seqinr read.fasta
#' @importFrom seqinr getSequence
#' @references paper under review
#' @seealso gquad
#' @examples
#' ## Predicting G quadruplex motif(s) (including overlaps) from raw nucleotide sequences
#' E1 <- c("TCTTGGGCATCTGGAGGCCGGAAT", "taggtgctgggaggtagagacaggatatcct")
#' gquadO(E1)
#'
#' ## Predicting G quadruplex motif(s) (including overlaps) from nucleotide sequences in fasta file
#' ## Not run: gquadO(x="Example.fasta", xformat = "fasta")
#'
#' ## Predicting G quadruplex motif(s) (including overlaps) from nucleotide sequences,
#' ## using GenBank accession numbers.
#' ## Internet connectivity is needed for this to work.
#' ## Not run: gquadO(c("BH114913", "AY611035"), xformat = "GenBank accession number(s)")

gquadO <- function(x, xformat = "default"){
  if(xformat == "default"){
    x1 <- gquadOverlap_main(x)
    return(x1)
  }

  if(xformat == "GenBank accession number(s)"){
    x2 <- read.GenBank(x)
    x3 <- sapply(x2, paste, collapse="")
    x4 <- gquadOverlap_main(x3)
    return(x4)
  }

  if(xformat == "fasta"){
    x5 <-read.fasta(x)
    x6 <- getSequence(x5, as.string = TRUE)
    x7 <- unlist(x6)
    x8 <- gquadOverlap_main(x7)
    return(x8)
  }

}
