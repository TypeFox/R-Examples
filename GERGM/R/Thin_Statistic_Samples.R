Thin_Statistic_Samples <- function(statistics){

  if(nrow(statistics) > 1 & length(unique(statistics$ttriads)) > 1){
    # we are going to test using ttriads statistics
    ttriads <- statistics$ttriads
    ar1 <- stats::cor(ttriads[2:length(ttriads)],ttriads[1:(length(ttriads)-1)])
    if(is.na(ar1)){
      cat("There was perfect autocorrelation, unable to thin statistics...\n")
    }else{
      #print(ar1)

      thin <- 1
      while(ar1 > .01){
        thin = thin +1
        thinSeq <- round(seq(1,length(ttriads),by=thin))
        thinDens <- ttriads[thinSeq]
        ar1 <- stats::cor(thinDens[2:length(thinDens)],thinDens[1:(length(thinDens)-1)])
        #print(ar1)
      }

      thinSeq <- round(seq(1,nrow(statistics),by=thin))

      if(length(thinSeq) > 99){
        cat("Thinning statistics to correct for autocorrelation in calculating fit diagnostics...\n")
        cat("Statistics were thinned by a factor of ",thin,", resulting in ",length(thinSeq)," samples.\n", sep = "")
        statistics <- statistics[thinSeq,]
      }else{
        cat("Thinning statistics to correct for autocorrelation in calculating fit diagnostics resulted in two few samples (less than 100), consider increasing the number of iterations of MCMC sampling...\n")
        cat("Statistics were thinned by a factor of ",thin,", resulting in ",length(thinSeq)," samples. Sample size was manually increased to 100.\n", sep = "")
        thinSeq <- round(seq(1,nrow(statistics),length.out = 100))
        statistics <- statistics[thinSeq,]
      }
    }
  }else{
    cat("Could not thin statistics, no variation...\n")
  }

  return(statistics)
}
