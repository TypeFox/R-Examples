#' Summary of a mrs object
#'
#' This function summarizes the output of the \code{\link{mrs}} function.
#' It provides the marginal prior and posterior null and 
#' the top regions of the representative tree.
#'
#' @param object A \code{mrs} object
#' @param rho Threshold for the posterior alternative probability. 
#'  All regions with posterior alternative probability larger 
#'  than \code{rho} are reported. Default is \code{rho = 0.5}. 
#' @param abs_eff Threshold for the effect size. All regions with 
#'  effect size larger than \code{abs_eff} in absolute value are reported. 
#'  Default is \code{abs_eff = 0}.
#' @param sort_by Define in which order the regions are reported. 
#'  The options are \code{sort_by = c("eff", "prob")} and 
#'  the default is \code{sort_by = "eff"}. 
#' @param ... Additional summary parameters. 
#' @return A \code{list} with information about the top regions.
#' @references Soriano J. and Ma L. (2014). Multi-resolution two-sample comparison 
#' through the divide-merge Markov tree. \emph{Preprint}. 
#'  \url{http://arxiv.org/abs/1404.3753}
#' @export
#' @S3method summary mrs
#' @examples
#' set.seed(1) 
#' n = 100
#' p = 2
#' X = matrix(c(runif(p*n/2),rbeta(p*n/2, 1, 4)), nrow=n, byrow=TRUE)
#' G = c(rep(1,n/2), rep(2,n/2))
#' object = mrs(X=X, G=G)
#' fit = summary(object, rho = 0.5, abs_eff = 0.1)
summary.mrs <- function(object, rho = 0.5, abs_eff = 0, sort_by = "eff", ...)
{
  if(class(object)!="mrs")
  {
    print("ERROR: object should be a mrs object.")
    return(0)
  }
  
  abs_max = apply( abs(object$RepresentativeTree$EffectSizes),1,max)
  top_regions = which( (object$RepresentativeTree$AltProbs > rho) & ( abs_max > abs_eff ) )
  
  if(length(top_regions)>0)
  {
    
    if(sort_by == "eff" ){
      plot_order = rev(top_regions[sort.int(abs_max[top_regions], index.return=TRUE)$ix])
    }else if(sort_by == "prob"){
      plot_order = rev(top_regions[sort.int(object$RepresentativeTree$AltProbs[top_regions], index.return=TRUE)$ix])
    }
    
    
    effect_size = matrix(object$RepresentativeTree$EffectSizes[plot_order,], nrow = length(plot_order))
    effect_size.names = rep(NA,object$Data$Groups)
    for( i in 1:object$Data$Groups )
      effect_size.names[i] = paste("Eff_Size",i, sep="_")
    colnames(effect_size) = effect_size.names  
    
    
    regions = matrix(object$RepresentativeTree$Regions[plot_order,], nrow = length(plot_order))
    regions.names = rep(NA,2*ncol(object$Data$X))
    for( i in 1:ncol(object$Data$X) )
      regions.names[c(2*i-1,2*i)] = c( paste("Min",i, sep="_"), paste("Max",i, sep="_"))    
    colnames(regions) = regions.names
    
    output = list(  Prior_Null = object$PriorGlobNull,
                    Posterior_Null = object$PostGlobNull,
                    Alt_Prob = object$RepresentativeTree$AltProbs[plot_order],
                    Effect_Size = effect_size,
                    Regions = regions,
                    Directions = object$RepresentativeTree$Directions[plot_order],
                    groups = object$Data$Groups,
                    p = ncol(object$Data$X),
                    Num_Regions = length(top_regions)                  
                  )
  }else{
    output = list(  Prior_Null = object$PriorGlobNull,
                    Posterior_Null = object$PostGlobNull,
                    Alt_Prob = NULL,
                    Effect_Size = NULL,
                    Regions = NULL,
                    Directions = NULL,
                    groups = object$Data$Groups,
                    p = ncol(object$Data$X),
                    Num_Regions = 0                  
    )
  }
  
  class(output) <- "summary.mrs"
  output
  
  
}
