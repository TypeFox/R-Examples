collapseMergeList <-
function(merge_list, merge_matrix, aggregate, K){
  calls <- which(merge_matrix %in% aggregate)
  
  # Generate deletion sequence
  old_refs <- c()
  row_indices <- c()
  col_indices <- c()
  for(i in seq(length(calls))){
    call <- calls[i]
    if(call <= K-1){
      row_indices[i] <- call
      col_indices[i] <- 1
    } else {
      row_indices[i] <- call - K + 1
      col_indices[i] <- 2
    }
    old_refs[i] <- merge_list[row_indices[i]][[1]][col_indices[i]]
  }
  deletion_seq <- data.frame(old_refs=old_refs,
                             row_indices=row_indices,
                             col_indices=col_indices)
  deletion_seq <- deletion_seq[order(old_refs),]
  
  # Perform collapsing: Insert new sequences
  for(i in seq(nrow(deletion_seq))){
    row_index <- deletion_seq$row_indices[i]
    col_index <- deletion_seq$col_indices[i]
    old_ref <- deletion_seq$old_refs[i]
    #merge_list[row_index][[1]] <- merge_list[row_index][[1]][-col_index]
    merge_list[row_index][[1]] <- c(merge_list[row_index][[1]],
                                    merge_list[old_ref][[1]])
  }
  
  # Perform collapsing: Delete old references
  for(row_index in names(merge_list)){
    delete <- which(merge_list[paste(row_index)][[1]] %in% aggregate)
    if(any(delete))
      merge_list[row_index][[1]] <- merge_list[row_index][[1]][-delete]
  }
  merge_list <- merge_list[-aggregate]
  return(merge_list)
}
