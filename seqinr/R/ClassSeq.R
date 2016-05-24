
##########################################################################
#
# SeqFastadna:
# 

as.SeqFastadna <- function(object, name = NULL, Annot = NULL){
  attributes(object) <- list(name = name, Annot = Annot)
  class(object) <- "SeqFastadna"	
  return(object)
}

is.SeqFastadna <- function(object) inherits(object, "SeqFastadna")

summary.SeqFastadna <- function(object, alphabet = s2c("acgt"), ...){
  length <- getLength(object)
  compo <- count(object, 1, alphabet = alphabet)
  return(list(length = length , composition = compo, GC = GC(object)))
}

#
# SeqFastaAA:
#

as.SeqFastaAA <- function(object, name = NULL, Annot = NULL){
  attributes(object) <- list(name = name, Annot= Annot)
  class(object) <- "SeqFastaAA"	
  return(object)
}

is.SeqFastaAA <- function(object) inherits(object, "SeqFastaAA")

summary.SeqFastaAA <- function(object,...){
  length <- getLength(object)
  compo <- table(factor(object, levels = levels(SEQINR.UTIL$CODON.AA$L)))
  return(list(length = length, composition=compo/length, AA.Property=AAstat(object,plot=FALSE)[[2]]))
}


#
# SeqAcnucWeb:
#

as.SeqAcnucWeb <- function(object, length, frame, ncbigc){
  attributes(object) <- list(length = as.numeric(length),
    frame = as.numeric(frame), ncbigc = as.numeric(ncbigc))
  class(object) <- "SeqAcnucWeb"
  return(object)
}

is.SeqAcnucWeb <- function(object) inherits(object, "SeqAcnucWeb")


print.SeqAcnucWeb <- function(x, ...)
{
  res <- c(x, attr(x, "length"), attr(x, "frame"), attr(x, "ncbigc"))
  names(res) <- c("name", "length", "frame", "ncbicg")
  print(res, ...)
}

#
# SeqFrag:
#


as.SeqFrag <- function(object, begin, end, name){
  attr(object, "seqMother") <- name
  attr(object, "begin") <- begin
  attr(object,"end") <- end
  class(object) <- "SeqFrag"
  return(object)
}

is.SeqFrag <- function(object) inherits(object, "SeqFrag")

#
# Query Acnuw Web (qaw class)
#

print.qaw <- function(x, ...)
{
  cat(x$nelem, x$type, "for", list1$call$query)
}
	
