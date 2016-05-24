Theta.cutoff <- function(trt.effect, f.d){
  #calculate theta(cutoff) for each entry in trt.effect. 
  # trt.effect should already be sorted by 
 
  out <- cumsum(-trt.effect[order(f.d)])/length(trt.effect)
  out[rank(f.d)]
}