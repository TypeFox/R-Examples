seq2aln <-
function(seq2add, aln, id="seq", file = "aln.fa", ...) {
  ##- Add a sequence 'seq2add' to an existing alignment 'aln'
  ##  Adds at the bottom of alignment
  cl <- match.call()
   
  if(!inherits(aln, "fasta")) 
     stop("Input 'aln' should be a 'fasta' object")

  tmp.seq = as.fasta(seq2add)
  if(nrow(tmp.seq$ali) > 1)
     warning("Multiple sequences in 'seq2add' should be pre-aligned")

  dots = list(...)
  # fixed arguments
  dots$profile = aln
  dots$outfile = file
    
  args = c(list(aln=seq2add, id=id), dots)
  naln = do.call(seqaln, args)

  naln$call = cl
  return(naln)
}

