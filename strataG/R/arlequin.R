#' @name arlequin
#' @title Read and Write Arlequin Files
#' @description Read and write an Arlequin-formatted files.
#' 
#' @param file filename for output file.
#' @param g a \linkS4class{gtypes} object.
#' @param title title for data in file.
#' @param data.type type of data. Can be "DNA", "RFLP", or "MICROSAT".
#' @param locus numeric or character designation of which locus to write.
#' 
#' @references Excoffier, L.G. Laval, and S. Schneider (2005) 
#'   Arlequin ver. 3.0: An integrated software package for population genetics 
#'   data analysis. Evolutionary Bioinformatics Online 1:47-50.\cr
#'   Available at \url{http://cmpg.unibe.ch/software/arlequin3/}
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom reshape2 dcast
#' @export
#' 
read.arlequin <- function(file) {
  arp <- scan(file, what = "character", quiet = TRUE)
  hd <- grep("HaplotypeDefinition", arp, ignore.case = TRUE)[1]
  left.brace <- grep("[{]", arp)
  right.brace <- grep("[}]", arp)
  seq.start <- min(left.brace[left.brace >= hd]) + 1
  seq.end <- min(right.brace) - 1
  seq.lt <- strsplit(arp[seq(seq.start + 1, seq.end, 2)], "")
  names(seq.lt) <- toupper(arp[seq(seq.start, seq.end - 1, 2)])
  seq.lt <- lapply(seq.lt, tolower)
  seq.lt <- lapply(seq.lt, function(x) {
    x[x == "?"] <- "n"
    x
  })
  sample.name <- grep("SampleName", arp, ignore.case = TRUE)
  left.brace <- left.brace[left.brace > seq.end + 1]
  right.brace <- right.brace[right.brace > seq.end + 1]
  eq.symbol <- grep("=", arp)
  freq.df <- do.call(rbind, lapply(sample.name, function(i) {
    strata.start <- min(eq.symbol[eq.symbol >= i])
    strata.end <- min(eq.symbol[eq.symbol > strata.start]) - 1
    strata <- paste(arp[strata.start:strata.end], collapse = " ")
    strata <- gsub("^[[:alnum:]=[:blank:]]*[\"]", "", strata)
    strata <- gsub("[\"]", "", strata)
    freq.start <- min(left.brace[left.brace > i]) + 1
    freq.end <- min(right.brace[right.brace > i]) - 1
    haps <- arp[seq(freq.start, freq.end - 1, by = 2)]
    freqs <- as.numeric(arp[seq(freq.start + 1, freq.end, by = 2)])
    data.frame(strata = strata, haps = haps, freqs = freqs, 
               stringsAsFactors = FALSE)
  }))
  freq.df <- dcast(freq.df, haps ~ strata, sum, value.var = "freqs")
  freq.df$haps <- toupper(freq.df$haps)
  list(freq.df = freq.df, dna.seq = seq.lt)
}

#' @rdname arlequin
#' @export
#' 
write.arlequin <- function(g, file = "gtypes.prj", title = "gtypes from R", 
                           data.type = c("DNA", "MICROSAT", "RFLP"),
                           locus = 1) {
  
  data.type <- match.arg(data.type)
  if(ploidy(g) > 1 & data.type == "DNA") {
    stop("'data.type' cannot be DNA if 'g' is not haploid.")
  }
  if(ploidy(g) == 1 & data.type %in% c("RFLP", "MICROSAT")) {
    stop("'data.type' cannot be RFLP or MICROSAT if 'g' is haploid.")
  }
  
  write("[Profile]", file = file)
  write(paste("Title = \"", title, "\"", sep = ""), file = file, append = TRUE)
  write(paste("NbSamples =", nInd(g)), file = file, append = TRUE)
  write(paste("DataType =", data.type), file = file, append = TRUE)
  write(paste("GenotypicData =", ifelse(ploidy(g) == 1, 0, 1)), 
        file = file, append = TRUE)
  
  write("[Data]", file = file, append = TRUE)
  if(!is.null(sequences(g))) {
    write("[[HaplotypeDefinition]]", file = file, append = TRUE)
    write("HaplListName=\"Haplotypes\"", file = file, append = TRUE)
    write("HaplList={", file = file, append = TRUE)
    dna <- as.character(as.matrix(sequences(g, locus)))
    for(x in rownames(dna)) {
      x.seq <- paste(dna[x, ], collapse = "")
      write(paste(x, x.seq), file = file, append = TRUE)
    }
    write("}", file = file, append = TRUE)
  }
  
  invisible(NULL)
}