#' dist.entropic returns an index of quantitative asymmetry and complexity
#' @title Return an index of quantitative asymmetry and complexity
#' @author Marc Girondot
#' @return A numeric value
#' @param l1 Set of measures at one side of an organism
#' @param l2 Set of measures at the other side of an organism
#' @param details If TRUE, will show the details of computing
#' @description Return an index of quantitative asymmetry and complexity. The higher is the value, 
#' the higher is the complexity (number of objects) and diversity (difference between them).\cr
#' The indice is based on the product of the average angular distance of Edwards (1971) with the 
#' geometric mean of the Shannon entropy H for both sides.
#' @references 
#' Edwards, A.W.F., 1971. Distances between populations on the basis of gene frequencies. Biometrics 27, 873–881.\cr
#' Shannon C.E. 1948 A mathematical theory of communication. Bell System Technical Journal 27(3), 379-423.\cr
#' @examples
#' l1 <- c(0.1, 0.1, 0.05, 0.2, 0.3, 0.25)
#' l2 <- c(0.2, 0.3, 0.5)
#' dist.entropic(l1, l2)
#' 
#' l1 <- c(0.1, 0.1, 0.05, 0.2, 0.3, 0.25)
#' l2 <- c(0.1, 0.1, 0.05, 0.2, 0.3, 0.25)
#' dist.entropic(l1, l2)
#' 
#' l1 <- c(0.2, 0.3, 0.5)
#' l2 <- c(0.2, 0.3, 0.5)
#' dist.entropic(l1, l2)
#' 
#' l1 <- c(0.2, 0.2, 0.2, 0.2, 0.2)
#' l2 <- c(0.2, 0.3, 0.5)
#' dist.entropic(l1, l2)
#' 
#' l1 <- c(0.2, 0.2, 0.2, 0.2, 0.2)
#' l2 <- c(0.3333, 0.3333, 0.3333)
#' dist.entropic(l1, l2)
#' 
#' l1 <- c(0.2, 0.2, 0.2, 0.2, 0.2)
#' l2 <- c(0.2, 0.2, 0.2, 0.2, 0.2)
#' dist.entropic(l1, l2)
#' 
#' l1 <- c(0.3333, 0.3333, 0.3333)
#' l2 <- c(0.3333, 0.3333, 0.3333)
#' dist.entropic(l1, l2)
#' @export

dist.entropic <- function(l1, l2, details=FALSE) {
  l1 <- c(l1, rep(0, max(0, length(l2)-length(l1))))
  l2 <- c(l2, rep(0, max(0, length(l1)-length(l2))))
  
  l1 <- l1/sum(l1)
  l2 <- l2/sum(l2)
  
  pp <- getFromNamespace(".permutations", ns="HelpersMG")(v=l1, r=length(l1), n=length(l1), 
                                                          set=FALSE)
  p <- (mean(apply(pp, 1, function(x) {sqrt(1-(sum(x*l2)))})))
  p1 <- -sum(l1[l1!=0]*log(l1[l1!=0]))
  p2 <- -sum(l2[l2!=0]*log(l2[l2!=0]))
  if (details) {
    return(c(Average.Edwards=p, Geometric.mean.Shannon=sqrt(p1*p2), dist.entropic=p*sqrt(p1*p2)))
  } else {
    return(p*sqrt(p1*p2))
  }
}



