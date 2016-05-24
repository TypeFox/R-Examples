#' Generate Goodness Of Fit plot from a GERGM object.
#'
#' @param GERGM_Object The object returned by the estimation procedure using the
#' GERGM function.
#' @param ... Additional Arguments can be passed in. Included for eventual compatibility with XERGM package.
#' @return A set of box plots where of simulated network statistics centered at the observed value for those statistics and normalized by their standard deviation. This aids in interpretation as the y-axis can be interpreted as the number of simulated-sample standard deviations above or below the observed statistic.
#' @export
GOF <- function(GERGM_Object,
                ...){
  #define colors
  UMASS_BLUE <- rgb(51,51,153,155,maxColorValue = 255)
  UMASS_RED <- rgb(153,0,51,255,maxColorValue = 255)

  if (GERGM_Object@simulation_only) {
    temp <- GERGM_Object@simulated_statistics_for_GOF

    boxplot(temp, medcol = UMASS_BLUE,
            xlab = "Network Statistic",
            ylab = "Statistic Values",
            main = "Simulated Network Statistics")
  } else {
    if (!GERGM_Object@directed_network) {
      #input is the GERGM_Object, we are going to normalize all statistics
      temp <- GERGM_Object@simulated_statistics_for_GOF

      temp2 <- apply(GERGM_Object@simulated_statistics_for_GOF,2,sd)
      for (i in 1:ncol(temp)) {
        temp[,i] <- temp[,i] - GERGM_Object@stats[2,i]
        temp[,i] <- temp[,i]/temp2[i]
      }

      #now we are only dealing with ttriads and twostars since this is an undirected network
      temp <- temp[,c(2,5,6)]
      colnames(temp)[1] <- "twostars"

      boxplot(temp, medcol = UMASS_RED,
              xlab = "Network Statistic",
              ylab = "Normalized Statistic Values",
              main = "Blue = Observed Statistic, Red = Simulated Mean")
      zero_line <- rep(0,length(GERGM_Object@stats[2, ]))
      zero_line <- zero_line[c(2,5,6)]
      zero_plot <- rbind(zero_line,zero_line)
      boxplot(zero_plot, add = T, medcol = UMASS_BLUE, names = F)
    } else {
      #input is the GERGM_Object, we are going to normalize all statistics
      temp <- GERGM_Object@simulated_statistics_for_GOF

      temp2 <- apply(GERGM_Object@simulated_statistics_for_GOF,2,sd)
      for(i in 1:ncol(temp)){
        temp[,i] <- temp[,i] - GERGM_Object@stats[2,i]
        temp[,i] <- temp[,i]/temp2[i]
      }

      boxplot(temp, medcol = UMASS_RED,
              xlab = "Network Statistic",
              ylab = "Normalized Statistic Values",
              main = "Blue = Observed Statistic, Red = Simulated Mean")
      zero_line <- rep(0,length(GERGM_Object@stats[2, ]))
      zero_plot <- rbind(zero_line,zero_line)
      boxplot(zero_plot, add = T, medcol = UMASS_BLUE, names = F)
    }
  }
}
