assignClusterNames <-
function(out, lab_no, beta_weights, vocab){
  
  members <- c()
  
  for(i in seq(length(out$children))){
    if (!('size' %in% names(out$children[[i]])))
      out$children[[i]] <- assignClusterNames(out$children[[i]], lab_no, 
                                              beta_weights, vocab)
  }
  
  for(i in seq(length(out$children))){
    members <- c(members, out$children[[i]]$topic_no)
    #labels <- c(labels, strsplit(out$children[[i]]$name, split=", ")[[1]])
  }
  
  out$topic_no <- members
  margins <- marginalize(members, beta_weights$beta, beta_weights$weights)
  labels <- vocab[margins$indices[1:lab_no]]
  
  out$name <- paste(sample(labels, lab_no), collapse=", ")
  
  return(out)
}
