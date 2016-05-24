#' Prepare the data
#' Simple wrapper to make sure that the matrices are sorted accordingly
#' @param H Host distance matrix 
#' @param P Parasite distance matrix 
#' @param HP Host-parasite association matrix, hosts in rows
#' @return A list with objects H, P, HP
#' @export
#' @examples 
#' data(gopherlice)
#' library(ape)
#' gdist <- cophenetic(gophertree)
#' ldist <- cophenetic(licetree)
#' D <- prepare_paco_data(gdist, ldist, gl_links)
prepare_paco_data <- function(H, P, HP)
{
   if(NROW(H) != NCOL(H))
      stop("H should be a square matrix")
   if(NROW(P) != NCOL(P))
      stop("P should be a square matrix")
   if(NROW(H) != NROW(HP)){
      warning("The HP matrix should have hosts in rows. It has been translated.")
      HP <- t(HP)
   }
   if (!(NROW (H) %in% dim(HP)))
     stop ("The number of species in H and HP don't match")
   if (!(NROW (P) %in% dim(HP)))
     stop ("The number of species in P and HP don't match")
   if (!all(rownames(HP) %in% rownames (H)))
     stop ("The species names H and HP don't match")
   if (!all(colnames(HP) %in% rownames (P)))
     stop ("The species names P and HP don't match")
   H <- H[rownames(HP),rownames(HP)]
   P <- P[colnames(HP),colnames(HP)]
   HP[HP>0] <- 1
   D <- list(H=H, P=P, HP=HP)
   class(D) <- "paco"
   return(D)
}
