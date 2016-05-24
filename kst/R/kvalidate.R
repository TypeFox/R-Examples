### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### kvalidate.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### dependencies: library(sets), library(proxy), library(relations)
###
### 2008-04-14: created
###

kvalidate <- function(x, rpatterns=NULL, method=c("gamma","percent","VC","DA")) {

   ### check x
   if (!inherits(x, "kstructure")) {
      stop(sprintf("%s must be of class %s.", dQuote("x"), dQuote("kstructure")))
   }

   ### check rpatterns
   dom <- kdomain(x)
   rp <- rpatterns
   rp <- rp[,order(colnames(rp))]
   if(is.null(rp) || as.set(colnames(rp))!=dom) {
      stop("Invalid response patterns.")
   }
   relmat <- relation_incidence(as.relation(x))
   relmat <- relmat[,order(colnames(relmat))]

   ### check method
   if (method!="gamma" & method!="percent" & method!="VC" & method!="DA") {
      stop("Invalid validation method.")
   } 

   if (method=="gamma") {
   ### gamma index
      nc <- 0
      nd <- 0
      for (i in seq_len(nrow(relmat))) {
         for (j in seq_len(ncol(relmat))) {
            if (relmat[i,j]==1) {
               nc <- nc+sum(rp[,i]==1 & rp[,j]==0)
               nd <- nd+sum(rp[,i]==0 & rp[,j]==1)
            }
         }
      }
      validate <- (nc-nd)/(nc+nd)

   } else if (method=="percent") {
   ### percent
      validate <- as.data.frame(colSums(rp)/nrow(rp)*100)
      colnames(validate) <- "%"

   } else if (method=="VC") {
   ### violational coefficient
      nd <- 0
      for (i in seq_len(nrow(relmat))) {
         for (j in seq_len(ncol(relmat))) {
            if (relmat[i,j]==1) {
               nd <- nd+sum(rp[,i]==0 & rp[,j]==1)
            }
         }
      }
      validate <- (1/(nrow(rp)*(sum(relmat)-ncol(rp))))*nd

   } else if (method=="DA") {
   ### distance agreement coefficient
      kmatrix <- 0+(t(sapply(x, function(z) dom %in% z)))
      colnames(kmatrix) <- dom
      Distances <- apply(dist(rp, kmatrix, method="Manhattan"), 1, min)
      ddat <- mean(Distances)
      ddat_dist <- table(Distances)
      validate <- NULL
      validate$ddat <- ddat
      validate$ddat_dist <- ddat_dist
      ps <- as.list(set(0,1)^length(kmatrix[1,]))
      psm <- mat.or.vec(length(ps),length(dom))
      colnames(psm) <- dom
      for (i in 1:length(ps)) {
         psm[i,] <- unlist(ps[[i]])
      }
      Distances <- apply(dist(psm, kmatrix, method="Manhattan"), 1, min)
      dpot <- mean(Distances)
      dpot_dist <- table(Distances)
      validate$dpot <- dpot
      validate$dpot_dist <- dpot_dist
      da <- ddat/dpot
      validate$DA <- da
   }

   ### return results
   validate

}
