SaveSequence = function(startYear, sequence, outFile) {
   nSeq = length(sequence)
   stopifnot(nSeq>=1)
   
   indexedSeq = matrix(0, nrow=nSeq, ncol=2)
    
   appendFlag = FALSE
   for (ndx in 1:nSeq) {
      if (is.na(sequence[ndx])) sequence[ndx] = 0
      cat(c(startYear + ndx - 1, ", ", sequence[ndx], "\n"), file=outFile, append=appendFlag)
      appendFlag = TRUE
      
      indexedSeq[ndx, 1] = startYear+ndx-1
      indexedSeq[ndx, 2] = sequence[ndx]
   }
   
   return(indexedSeq)
}