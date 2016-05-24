## Summarize function used in generic methods
## Return mean, median, quantiles

summarize <-
function(x, .prob, .digits){
  out <-
    c(round(mean(x), .digits), 
      round(median(x), .digits),
      round(quantile(x, probs = c(0, .prob) + (1 - .prob) / 2), .digits))
  return(out)
}
