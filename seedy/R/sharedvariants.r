sharedvariants <- function(deepseqmat) {
      sharers <- apply(deepseqmat,1,function(x){sum(!x%in%c(0,1))})
      if (sum(sharers==0)>0) {
         sharers <- sharers[-which(sharers==0)]
      }
      sharedpolylocs <- which(sharers>1)
      polymorphic.loci <- as.numeric(sapply(strsplit(rownames(deepseqmat), " "), function(x)x[2]))
      sharedpolys <- polymorphic.loci[sharedpolylocs]
      cand <- list()
      n <- ncol(deepseqmat)
      scoremat <- matrix(0,n,n) # number of shared variants between hosts
      for (i in 1:n) {
        if (sum(polymorphic.loci[which(!deepseqmat[,i]%in%c(0,1))]%in%sharedpolys)>1) { 
          # shared polys for this host
          B <- polymorphic.loci[which(!deepseqmat[,i]%in%c(0,1))]
          B <- which(polymorphic.loci%in%B[which(B%in%sharedpolys)])
          candidates <- NULL
          for (j in 1:length(B)) {
            candidates <- c(candidates, which(!deepseqmat[B[j],]%in%c(0,1)))
          }
          candidates <- unique(candidates)
          candidates <- candidates[-which(candidates==i)]
          if (length(candidates)>0) {
            score <- numeric(length(candidates))
            for (j in 1:length(score)) {
              score[j] <- length(intersect(which(!deepseqmat[,candidates[j]]%in%c(0,1)),
                                    which(!deepseqmat[,i]%in%c(0,1))))
            }
            cand[[i]] <- candidates
            scoremat[i,candidates] <- score
          }
        }
      }
      rownames(scoremat) <- paste("P",as.numeric(sapply(strsplit(colnames(deepseqmat), " "), function(x)x[2])), sep="")
      colnames(scoremat) <- rownames(scoremat)
      return(scoremat)
}