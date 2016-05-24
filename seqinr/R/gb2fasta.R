########################################################################
#
#                                   gb2fasta
#
#                   single entry genbank to fasta conversion
#
########################################################################
gb2fasta <- function(source.file = "ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Acinetobacter_ADP1_uid61597/NC_005966.gbk",
 destination.file = "Acinetobacter_ADP1_uid61597.fasta")
{
  input <- readLines(source.file)
  head <- input[1]
  head <- unlist(strsplit(head, split=" "))
  head <- head[nchar(head) > 0]
  seqname <- head[2]
  seqsize <- as.integer(head[3])
  outheader <- sprintf(">%s %d bp", seqname, seqsize)

  confile <- file(destination.file, open="w")
  writeLines(outheader, confile)
  #
  # Look for sequence position:
  #
  debut <- which(substring(input,1,6)=="ORIGIN") + 1
  if( length( debut ) > 1 )
    stop("Multiple entries not yet implemented !")
  fin <- which(substring(input,1,2)=="//") - 1
  if( length( fin ) > 1 )
    stop("Multiple entries not yet implemented !")
    
  input <- input[debut:fin]
  input <- sapply(input, function(x) {
    return(paste(substr(x,11,20),substr(x,22,31),substr(x,33,42),substr(x,44,53),
                 substr(x,55,64),substr(x,66,75),sep="",collapse="")) } )
  names(input)<-NULL
  writeLines(input, confile )
  close(confile)
}


