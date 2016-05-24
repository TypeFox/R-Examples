########################################################################
##get enrichment score for each miR set
EnrichmentScore <- function(miR.list, miR.set, weighted.score.type = 1, 
correl.vector = NULL) {  


   tag.indicator <- sign(match(miR.list, miR.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
   no.tag.indicator <- 1 - tag.indicator 
   N <- length(miR.list) 
   Nh <- length(miR.set) 
   Nm <-  N - Nh 
   if (weighted.score.type == 0) {
      correl.vector <- rep(1, N)
   }
   alpha <- weighted.score.type
   correl.vector <- abs(correl.vector**alpha)
   sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
   norm.tag    <- 1.0/sum.correl.tag
   norm.no.tag <- 1.0/Nm
   RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
   max.ES <- max(RES)
   min.ES <- min(RES)
   if (max.ES > - min.ES) {
      ES <- signif(max.ES, digits = 5)
      arg.ES <- which.max(RES)
   } else {
      ES <- signif(min.ES, digits=5)
      arg.ES <- which.min(RES)
   }
   return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
}



#############################################################################
EnrichmentScore2 <- function(miR.list, miR.set, weighted.score.type = 1, correl.vector = NULL) {  


   N <- length(miR.list) 
   Nh <- length(miR.set) 
   Nm <-  N - Nh 

   loc.vector <- vector(length=N, mode="numeric")
   peak.res.vector <- vector(length=Nh, mode="numeric")
   valley.res.vector <- vector(length=Nh, mode="numeric")
   tag.correl.vector <- vector(length=Nh, mode="numeric")
   tag.diff.vector <- vector(length=Nh, mode="numeric")
   tag.loc.vector <- vector(length=Nh, mode="numeric")

   loc.vector[miR.list] <- seq(1, N)
   tag.loc.vector <- loc.vector[miR.set]

   tag.loc.vector <- sort(tag.loc.vector, decreasing = FALSE)

   if (weighted.score.type == 0) {
      tag.correl.vector <- rep(1, Nh)
   } else if (weighted.score.type == 1) {
      tag.correl.vector <- correl.vector[tag.loc.vector]
      tag.correl.vector <- abs(tag.correl.vector)
   } else if (weighted.score.type == 2) {
      tag.correl.vector <- correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
      tag.correl.vector <- abs(tag.correl.vector)
   } else {
      tag.correl.vector <- correl.vector[tag.loc.vector]**weighted.score.type
      tag.correl.vector <- abs(tag.correl.vector)
   }

   norm.tag <- 1.0/sum(tag.correl.vector)
   tag.correl.vector <- tag.correl.vector * norm.tag
   norm.no.tag <- 1.0/Nm
   tag.diff.vector[1] <- (tag.loc.vector[1] - 1) 
   tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
   tag.diff.vector <- tag.diff.vector * norm.no.tag
   peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
   valley.res.vector <- peak.res.vector - tag.correl.vector
   max.ES <- max(peak.res.vector)
   min.ES <- min(valley.res.vector)
   ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)

   return(list(ES = ES))

}
