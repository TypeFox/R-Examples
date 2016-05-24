### For reading.

read.seq <- function(file.name, forceDNAtolower = FALSE,
    convertDNAtoupper = TRUE){
  ret <- seqinr::read.fasta(file.name, forceDNAtolower = forceDNAtolower)

  ### Make sure everything is in upper case.
  if(convertDNAtoupper){
    ret <- lapply(ret, function(x){ dna.low2up(x) })
  }

  ret
} # End of read.seq().

read.phi.df <- function(file.name, header = TRUE, sep = "\t", quote = ""){
  ret <- read.table(file.name, header = header, sep = sep, quote = quote,
                    stringsAsFactors = FALSE)
  ret$phi <- as.double(ret$phi)
  ret
} # End of read.phi.df().

write.seq <- function(seq.data, file.name){
  seqinr::write.fasta(seq.data, names(seq.data), file.name)
  invisible()
} # End of write.seq().

write.phi.df <- function(phi.df, file.name){
  phi.df$phi <- as.double(phi.df$phi)
  write.table(phi.df, file.name, quote = FALSE, sep = "\t", row.names = FALSE)
  invisible()
} # End of write.phi.df().

