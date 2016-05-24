"write.fasta" <-
function(alignment=NULL,
         ids=NULL,
         seqs=alignment$ali,
         file,
         append = FALSE) {
  
  if (is.null(seqs))
    stop("write.fasta: please provide a 'seqs' or 'alignment' input object")
  
  if (!is.null(alignment)) {
    if (is.null(alignment$id) | is.null(alignment$ali)) {
      stop("write.fasta: 'alignment' should be a list with '$id' and '$ali'components")
    }
    if (is.null(ids)) {
      ids=alignment$id
    }
  } else {
    if (is.null(ids)) {
      n.ids <- nrow(seqs)
      if(is.null(n.ids)) { n.ids=1 }
      ids=seq( 1, length=n.ids )
    }
  } 

  if (!append) {
    ##file.remove(file, showWarnings = FALSE)
    suppressWarnings( file.remove(file) )
  }
  nseqs <- length(ids)
  if (nseqs == 1) {
     # change for shortening lines (<=60) - Xinqiu
    cat(">", ids, "\n", file = file, append = TRUE, sep = "")
    cat(seqs, file = file, append = TRUE, sep = "", fill = 60)
#    cat(">", ids, "\n", seqs, "\n", file = file,
#        append = TRUE, sep = "")
  }
  else {
    for (i in 1:nseqs) {
       cat(">", ids[i], "\n", file = file, append = TRUE, sep = "")
       cat(seqs[i,], file = file, append = TRUE, sep = "", fill = 60)
#      cat(">", ids[i], "\n", seqs[i,],
#          "\n", file = file, append = TRUE, sep = "")
    }
  }
}
