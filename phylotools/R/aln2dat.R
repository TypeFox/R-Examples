#### Function aln2dat as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

aln2dat <-
function(aln){
   aln <- aln[-1]
   aln2 <- aln[(regexpr(" ", aln) > 0)&(!grepl("[*]", aln))]
   seqs <- gsub(" ","",substring(aln2, regexpr(" ", aln2), nchar(aln2)))
   nam <- substring(aln2, 1, regexpr(" ", aln2)-1)
   seqnam <- unique(nam)
   nline <- length(aln2)
   nsp <- length(seqnam)
   blocks <- nline/nsp
   species <- c()
   for (i in 1:nsp){
      linessp <- i + nsp*(1:blocks)
      species[i] <- paste(seqs[linessp], collapse = "")
   }
   result <- data.frame(names = seqnam, sequences = species)
   result
}

