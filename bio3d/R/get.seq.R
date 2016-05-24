`get.seq` <-
function(ids, outfile="seqs.fasta", db="nr") {
  ## Download FASTA format sequences from the NR or
  ## SWISSPROT/UNIPROT databases via their gi or
  ## SWISSPROT identifer number

  if( !(db %in% c("nr", "swissprot", "uniprot")) )
    stop("Option database should be one of nr, swissprot or uniprot")


  ids <- unique(ids)
  
  if(db=="nr") {
    get.files <- paste("http://www.ncbi.nlm.nih.gov/",
                       "sviewer/viewer.fcgi?db=protein&val=",
                       ids,"&report=fasta&retmode=text", sep="")

##    ## Old pre Oct-18th-2010 format URL
##    get.files <- paste("http://www.ncbi.nlm.nih.gov/entrez/",
##                       "viewer.fcgi?db=protein&val=",
##                       ids, "&dopt=fasta&sendto=t", sep="")

  } else {
    if(any(nchar(ids) != 6)) {
      warning("ids should be standard 6 character SWISSPROT/UNIPROT formart: trying first 6 char...")
      ids <- substr(basename(ids),1,6)
    }
    ids <- unique(ids)
    get.files <- file.path("http://www.uniprot.org/uniprot",
                           paste(ids, ".fasta", sep="") )
  }

  ## Remove existing file
  if(file.exists(outfile)) {
    warning(paste("Removing existing file:",outfile))
    unlink(outfile)
  }

  rtn <- rep(NA, length(ids))
  for(k in 1:length(ids)) {
    rtn[k] <- download.file( get.files[k], outfile, mode="a" )
  }

  names(rtn) <- ids
  if(all(!rtn)) {
    return(read.fasta(outfile))
  } else {
    warning("Not all downloads were sucesfull, see returned values")
    return(rtn)
  }
}


