#' @name mega
#' @title Read and Write MEGA
#' @description Read and write MEGA formatted files.
#' 
#' @param file a MEGA-formatted file of sequences.
#' @param g a \linkS4class{gtypes} object.
#' @param title title for data in file.
#' @param line.width width of sequence lines.
#' @param locus number or name of locus to write.
#' 
#' @return for \code{read.mega}, a list of:
#' \tabular{ll}{
#'   \code{title} \tab title of MEGA file.\cr
#'   \code{dna.seq} \tab a list of DNA sequences.\cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @aliases mega, MEGA
#' @export
#' 
read.mega <- function(file) {  
  mega.file <- scan(file, what = "character", sep = "\n", 
                    strip.white = TRUE, quiet = TRUE)
  markers <- c(grep("#", mega.file), length(mega.file) + 1)
  title <- paste(mega.file[2:(markers[2] - 1)], collapse = " ")
  seq.df <- as.data.frame(t(sapply(2:(length(markers) - 1), function(i) {
    id <- sub("#", "", mega.file[markers[i]])
    seq.start <- markers[i] + 1
    seq.stop <- markers[i + 1] - 1
    dna.seq <- paste(mega.file[seq.start:seq.stop], collapse = "")
    c(id = id, sequence = dna.seq)
  })), stringsAsFactors = FALSE)
  dna <- strsplit(seq.df$dna.seq, "")
  names(dna) <- dna$id
  list(title = title, dna = as.DNAbin(dna))
}

#' @rdname mega
#' @export
#' 
write.mega <- function(g, file = NULL, title = NULL, line.width = 60, locus = 1) {
  if(is.null(file)) {
    file <- paste(description(g), ".meg", sep = "")
    file <- gsub("[[:punct:]]", ".", file)
  }
  if(is.null(title)) title <- description(g)
  dna <- sequences(g, locNames(g)[locus])
  dna <- as.character(as.matrix(dna))
  
  write("#MEGA", file)
  write(paste("title:", title, sep = ""), file, append = TRUE)
  write("", file, append = TRUE)
  for(x in names(dna)) {
    write(paste("#", x, sep = ""), file, append = TRUE)
    mt.seq <- dna[x, ]
    for(j in seq(1, length(mt.seq), by = line.width)) {
      seq.line <- paste(mt.seq[j:(j + line.width - 1)], collapse = "")
      write(seq.line, file, append = TRUE)
    }
    write("", file, append = TRUE)
  }
}