diffsred <-function(x,nobj){
################################################################
# calculate differences from non-paired comparison responses
# (i.e. ratings/Likert or rankings)
# Sinclair ordering: (1-2)(1-3)(2-3)(1-4)
# reduces diffs to -1/0/1
################################################################

   ncomp<-nobj*(nobj-1)/2
   diffs<-matrix(c(0:0),nrow(x),ncomp)
   col<-1
   for (j in 2:nobj) {
       for (i in 1:(j-1) ){
           diffs[, col]  <- ifelse(x[,i]>x[,j],1,ifelse(x[,i]<x[,j],-1,0))
           col <- col + 1
       }
   }
   diffs
}
