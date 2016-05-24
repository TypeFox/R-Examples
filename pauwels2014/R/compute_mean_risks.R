compute_mean_risks <- function(mean_risks, legend){
  Reduce(
    lapply(
      1:length(mean_risks), 
      FUN = function(k){
        data.frame(
          cost = mean_risks[[k]][2,], 
          risk = mean_risks[[k]][1,], 
          type = legend, 
          chain = as.factor(k)
        )
      }
    ), 
    f = rbind
  )
}
