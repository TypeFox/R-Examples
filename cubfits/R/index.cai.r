### This file is to compute cai (codon adaptive index).
###
### CAI = exp( 1/L \sum_{l = 1}^L log( w_i(l) ) )
###
### w_i(l) = f_i / max( f_j )
###
### i, j \in synonymous codons for amino acid

calc_cai_values <- function(y, y.list, w = NULL){
  ### Get w
  if(is.null(w)){
    w <- do.call("c", lapply(1:length(y),
                        function(i.aa){
                          tmp <- colSums(y[[i.aa]])
                          tmp <- tmp / max(tmp)
                          names(tmp) <- colnames(y[[i.aa]])
                          tmp
                        }))
  }

  ### Update w's names.
  w.names <- names(w)
  if(is.null(w.names)){
    stop("w should be named.")
  } else{
    w.names <- codon.low2up(w.names)
    names(w) <- w.names
  }

  ### Sort w by y's codon names. 
  codon.names <- do.call("c", lapply(y.list[[1]], function(aa){ names(aa) }))
  new.w <- NULL 
  for(i in codon.names){
    id <- which(w.names == i)
    if(length(id) != 1){
      stop("match incorrect.")
    } else{
      new.w <- c(new.w, w[id])
    }
  }
  names(new.w) <- codon.names

  ### Comput CAI.
  log.w <- log(new.w)
  CAI.values <- do.call("c", lapply(1:length(y.list),
                               function(i.seq){
                                 tmp <- do.call("c", y.list[[i.seq]])
                                 exp(log.w %*% tmp / sum(tmp))
                               }))
  names(CAI.values) <- names(y.list)

  ### Return.
  ret <- list(CAI = CAI.values, w = new.w)
  ret
} # End of calc_cai_values().

