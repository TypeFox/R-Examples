"write.pir" <-
function(alignment = NULL,
         ids = NULL,
         seqs = alignment$ali,
         pdb.file = NULL,
         chain.first = NULL,
         resno.first = NULL,
         chain.last = NULL, 
         resno.last = NULL, 
         file,
         append = FALSE) {
  
  if (is.null(seqs))
    stop("write.pir: please provide a 'seqs' or 'alignment' input object")
  
  if (!is.null(alignment)) {
    if (is.null(alignment$id) | is.null(alignment$ali)) {
      stop("write.pir: 'alignment' should be a list with '$id' and '$ali'components")
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

  if(is.null(pdb.file)) pdb.file = rep("", nrow(seqs))
  if(is.null(chain.first)) chain.first = rep("", nrow(seqs))
  if(is.null(resno.first)) resno.first = rep("", nrow(seqs))
  if(is.null(chain.last)) chain.last = rep("", nrow(seqs))
  if(is.null(resno.last)) resno.last = rep("", nrow(seqs))

  if (!append) {
    ##file.remove(file, showWarnings = FALSE)
    suppressWarnings( file.remove(file) )
  }
  nseqs <- length(ids)
  if (nseqs == 1) {
    head = "sequence"
    if(pdb.file!="") {
       head = "structureX"
       if(chain.first == "") chain.first = "@"
       if(resno.first == "") resno.first = "FIRST"
       if(resno.last == "") resno.last = paste("+", sum(!is.gap(seqs)), sep="") 
    }
    head = paste(head, pdb.file, resno.first, chain.first,
            resno.last, chain.last, "", "", "", "", sep=":")
 
    # change for shortening lines (<=60) - Xinqiu
    cat(">P1;", ids, "\n", file = file, append = TRUE, sep = "")
    cat(head, "\n", file = file, append = TRUE)
    cat(seqs, "*", file = file, append = TRUE, sep = "", fill = 60)
  }
  else {
    for (i in 1:nseqs) {
       head = "sequence"
       if(pdb.file[i]!="") {
          head = "structureX"
          if(chain.first[i] == "") chain.first[i] = "@"
          if(resno.first[i] == "") resno.first[i] = "FIRST"
          if(resno.last[i] == "") resno.last[i] = paste("+", sum(!is.gap(seqs[i,])), sep="") 
       }
       head = paste(head, pdb.file[i], resno.first[i], chain.first[i],
               resno.last[i], chain.last[i], "", "", "", "", sep=":")
 
       cat(">P1;", ids[i], "\n", file = file, append = TRUE, sep = "")
       cat(head, "\n", file = file, append = TRUE)
       cat(seqs[i,], "*", file = file, append = TRUE, sep = "", fill = 60)
    }
  }
}
