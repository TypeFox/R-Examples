random.msa <- function (nb.seq = 100, id = "SEQ", nb.pos = 100, gap = FALSE, aa.strict = FALSE,align = NULL, align.replace = TRUE) {

  #one letter codes for amino acids 
  aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
    "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "J", "X")
  replace <- TRUE
  #remove ambiguous amino acids
  if (aa.strict)
    aa <- aa[1:20]

  if (gap)
    aa <- c(aa, "-")

  if(!is.null(align)){
    if (!inherits(align, "align")) 
        stop("mmds is not a 'align' object")
    aa <-as.vector(unlist(align))
    replace<-align.replace
  }

  msa <- lapply(seq_len(nb.seq), function (i) {sample(aa, nb.pos, replace = replace)})

  msa.names <- paste(id, seq_len(nb.seq), sep = "")

  names(msa) <- msa.names

  return(msa)
}