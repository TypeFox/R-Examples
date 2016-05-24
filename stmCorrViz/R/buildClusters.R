buildClusters <-
function(out, current, merge_list, labels, full_labels, 
                          thoughts, topic_proportions){
  
  # Recursive definition
  for(i in seq(length(current))){
    out[[i]] <- list()
    out[[i]]$name <- current[i]
    if(current[i] > 0){
      out[[i]]$children <- buildClusters(list(), 
                                         merge_list[paste(current[i])][[1]],
                                         merge_list, labels=labels, full_labels,
                                         thoughts, topic_proportions)
      out[[i]]$name <- current[i]
    } else {
      out[[i]]$size <- 1800
      out[[i]]$name <- labels[-current[i]]
      out[[i]]$topic_no <- -current[i]
      out[[i]]$thought_1 <- iconv(thoughts$docs[-current[i]]$Topic[1], to='utf-8', sub="")
      out[[i]]$thought_2 <- iconv(thoughts$docs[-current[i]]$Topic[2], to='utf-8', sub="")
      out[[i]]$prob <- paste(full_labels$prob[-current[i],], collapse = ", ")
      out[[i]]$frex <- paste(full_labels$frex[-current[i],], collapse = ", ")
      out[[i]]$lift <- paste(full_labels$lift[-current[i],], collapse = ", ")
      out[[i]]$score <- paste(full_labels$score[-current[i],], collapse = ", ")
      out[[i]]$proportion <- format(round(topic_proportions[-current[i]], 2))
    }
  }
  return(out)
}
