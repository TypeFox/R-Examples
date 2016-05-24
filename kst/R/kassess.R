### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### kassess.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dependencies: library(sets)
###
### 2008-04-17: created
###

kassess <- function(x, rpatterns=NULL, method="deterministic") {

   ### check x
   if (!inherits(x, "kstructure")) {
      stop(sprintf("%s must be of class %s.", dQuote("x"), dQuote("kstructure")))
   }

   ### check rpatterns
   dom <- kdomain(x)
   rp <- rpatterns
   if(is.null(rp) || as.set(colnames(rp))!=dom) {
      stop("Invalid response patterns.")
   }

   method <- match.arg(method)

   ### deterministic assessment
   if (method=="deterministic") {
      assess <- NULL
      for (i in seq_len(nrow(rp))) {
         pstates <- as.list(x)
         while (length(pstates)>1) {
            problem <- which.min((table(unlist(pstates))-length(pstates)/2)^2)
            states <- grep(names(problem), pstates)
            if (rp[i,names(problem)]==1) {
               pstates <- pstates[states]
            } else {
               pstates <- pstates[-states]
            }
         }
         assess <- c(assess, pstates)
      }
      names(assess) <- paste("Respondent", seq_along(assess), sep="")
   }
   
   ### return results
   assess

}
