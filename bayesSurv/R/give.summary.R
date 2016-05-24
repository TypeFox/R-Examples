#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2005)                              ####
####                                                 ####
#### FILE:       give.summary.R                      ####
####                                                 ####
#### FUNCTIONS:  give.summary                        ####
#########################################################

### ======================================
### give.summary
### ======================================
give.summary <- function(sample, probs=c(0.5, 0.025, 0.975)){
  if (is.null(dim(sample))){
    means <- mean(sample)
    quant <- quantile(sample, probs=probs)
    negative <- sum(sample < 0)/length(sample)
    positive <- sum(sample > 0)/length(sample)
    p.value <- 2*min(negative, positive)
    summ <- matrix(c(means, quant, p.value), ncol=1)
    rownames(summ) <- c("means", names(quant), "p.value")
    colnames(summ) <- "sample"
  }
  else{
    means <- apply(sample, 2, mean)
    quant <- apply(sample, 2, quantile, probs=probs)
    negative <- apply(sample < 0, 2, sum)/dim(sample)[1]
    positive <- apply(sample > 0, 2, sum)/dim(sample)[1]
    p.value <- 2*apply(rbind(negative, positive), 2, min)
    summ <- rbind(means, quant, p.value)    
  }
  
  return(summ)
}
