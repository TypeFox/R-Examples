lrSib = function(sib1, sib2, Freqs  = NULL, nLoci = length(sib1)/2,
                 f = NULL, n = NULL){
  
  if(is.null(Freqs) & (is.null(f) | is.null(n))){
    stop("Must specify either Freqs or both f and n")
  }
  
  if(is.null(Freqs)){
    loc = rep(1:length(n), n)
    Freqs = split(f, n)
  }
  
  
  lr = .lrSib(sib1, sib2, Freqs$freqs)
  
  return (lr)
}
