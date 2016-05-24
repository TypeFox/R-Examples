#' Contribution of individual links
#' @param D A list returned by proc_analysis
#' @param .parallel if \code{TRUE}, calculate the jacknife contribution in parallel using the backend provided by foreach
#' @return A list with added object jacknife, containing the mean and upper CI values for each link
#' @export
paco_links <- function(D, .parallel = FALSE)
{
   HP.ones <- which(D$HP > 0, arr.ind=TRUE)
   SQres.jackn <- matrix(rep(NA, sum(D$HP)^2), sum(D$HP))# empty matrix of jackknifed squared residuals
   colnames(SQres.jackn) <- paste(rownames(D$proc$X),rownames(D$proc$Yrot), sep="-") #colnames identify the H-P link
   t.critical = stats::qt(0.975,sum(D$HP)-1) #Needed to compute 95% confidence intervals.
   nlinks <- sum(D$HP)

   # In parallel

   plyr::adply(1:nlinks, 1, function(x) single_paco_link(D, HP.ones, x), .parallel=.parallel)

   SQres.jackn <- SQres.jackn^2 #Jackknifed residuals are squared
   SQres <- (stats::residuals(D$proc))^2 # Vector of original square residuals
   #jackknife calculations:
   SQres.jackn <- SQres.jackn*(-(sum(D$HP)-1))
   SQres <- SQres*sum(D$HP)
   SQres.jackn <- t(apply(SQres.jackn, 1, "+", SQres)) #apply jackknife function to matrix
   phi.mean <- apply(SQres.jackn, 2, mean, na.rm = TRUE) #mean jackknife estimate per link
   phi.UCI <- apply(SQres.jackn, 2, stats::sd, na.rm = TRUE) #standard deviation of estimates
   phi.UCI <- phi.mean + t.critical * phi.UCI/sqrt(sum(D$HP))
   D$jackknife <- list(mean = phi.mean, upper = phi.UCI)
   return(D)
}

#PACo setting the ith link = 0
single_paco_link <- function (D, HP.ones, i) {
  HP_ind <- D$HP
  HP_ind[HP.ones[i,1],HP.ones[i,2]]=0
  PACo.ind <- add_pcoord(list(H=D$H, P=D$P, HP=HP_ind))
  Proc.ind <- vegan::procrustes(X=PACo.ind$H_PCo, Y=PACo.ind$P_PCo)
  res.Proc.ind <- c(residuals.paco(Proc.ind))
  res.Proc.ind <- append(res.Proc.ind, NA, after= i-1)
}
