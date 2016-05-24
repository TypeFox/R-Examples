#' Apply a differentiation statistic to a bootstrap sample
#'
#' This function applies a differentiation statistic (eg, D_Jost, Gst_Hedrick or 
#' Gst_Nei) to a list of genind objects, possibly produced with
#' chao_bootsrap or jacknife_populations. The resulting list contains a matrix
#' of values with the statistic for each locus as well as a global estimate 
#' for every object in the sample. Additionally, mean and 95% confidence 
#' intervals are calculated for each set of statisics A custom print method 
#' that displays these summaries is provided.
#' 
#'
#' @param bs list of genind objects
#' @param statistic differentiation statistic to apply (the function itself, 
#' as with apply family functions)
#' @family resample
#' @return per.locus:  matrix of statistics calculated for each locus and each 
#' bootstrap replication
#' @return global.het: vector of global estimates calculated from overall 
#' heterozygosity 
#' @return global.het: vector of global estimates calculated from harmonic
#' mean of statistic (only applied to D_Jost)
#' @return summary.loci: matrix containing mean, .025 and 0.975 percentile and
#' varaince of statisic for each locus
#' @return summary.global_het: mean, .025 and 0.975 percentile and variance for
#' global estimate variance of statistic for each locus based on heterozygosity
#' @return summary.global_harm: mean, .025 and 0.975 percentile and variance for
#' global estimate variance of statistic for each locus based on harmonic mean
# (only applies to D_Jost)
#' @importFrom stats var
#' @importFrom stats quantile
#' @export
#' @examples
#'\dontrun{  
#' data(nancycats)
#' bs <- chao_bootstrap(nancycats)
#' summarise_bootstrap(bs, D_Jost)
#'}


summarise_bootstrap <- function(bs, statistic){
  nreps <- length(bs)
  stats <- sapply(bs, statistic)
  loc_stats <- do.call(rbind, stats["per.locus",])

  res <-list("per.locus"= loc_stats,
             "global.het"=unlist(stats[2,])
            )
  #Only D_Jost has another estimate of global value
  if(identical(statistic, D_Jost)){
    res$global.harm <- unlist(stats[3,])
    }
  summarise <- function(x){
    return(c(mean=mean(x), 
             quantile(x, c(0.025, 0.975), na.rm=TRUE),
             variance=var(x) 
            ))
  }
  res$summary.loci <- apply(loc_stats, 2, summarise)
  res$summary.global.het <- summarise(res$global.het)
  if(identical(statistic, D_Jost)){
    res$summary.global.harm <- summarise(res$global.harm)
    }
  class(res) <- "summarised_bs"
  return(res)
}

#' @export 

print.summarised_bs <- function(x, ...){
  
  print_line <- function(x){
  x <- round(x, 4)
  return(paste(x[1], "\t(", x[2],"-", x[3], ")\n", sep="")) 
  }
  
  loc.names <- colnames(x$per.locus)
  loc.results <- apply(x$summary.loci, 2, print_line)
  cat("\nEstimates for each locus\n")
  cat("Locus\tMean\t 95% CI\n")
  for(i in 1:length(loc.names)){
    cat(paste(loc.names[i], loc.results[i], sep="\t"))
  }
  cat("\nGlobal Estimate based on average heterozygosity\n")
  cat(print_line(x$summary.global.het))
  if(!is.null(x$summary.global.harm)){
    cat("\nGlobal Estimate based on harmonic mean of statistic\n")
    cat(print_line(x$summary.global.harm))
  }
}

