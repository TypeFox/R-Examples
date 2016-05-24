## Adds variables to the data element of a seqelist object
##   ntrans : the length of the event sequence (number of transitions)
##   nevent : the total number of events
##   avg.occ: average occurrences of the subsequence by sequence

seqentrans <- function(fsubseq, avg.occ=FALSE){
     fsubseq$data$ntrans <- sapply(strsplit(as.character(fsubseq$subseq), "-"),
     length)
     fsubseq$data$nevent <- fsubseq$data$ntrans - 1 +
     sapply(strsplit(as.character(fsubseq$subseq), ","), length)
     if (avg.occ){
         wtot <- sum(seqeweight(fsubseq$seqe))
         fsubseq$data$Avg.occ <- fsubseq$data$Count/wtot
     }
     return(fsubseq)
}
