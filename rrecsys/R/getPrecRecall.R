# Recall, precision, TP, FP, atc. on a given user.

getPrecRecall <- function(data, rec_items, fold_items_x_user, goodRating) {
  res <- list(TP = 0, FP = 0, TN = 0, FN = 0, precision = 0, recall = 0, F1 = 0)

  if (length(rec_items) != 0) {
    # match test set items with recommended items
    match_TS <- which(rec_items %in% fold_items_x_user)
    difference_TS <- which(!fold_items_x_user %in% rec_items)

    res$TP <- sum(data[rec_items[match_TS]] >= goodRating)

    res$FN <- sum(data[fold_items_x_user[difference_TS]] >= goodRating)
    
    res$FP <- length(rec_items) - res$TP
    
    res$TN <- length(data) - res$TP - res$FN - res$FP
    
    
    if (sum(data[fold_items_x_user] >= goodRating) == 0) {#all the items in the fold disliked by the user
      res$precision <- 1
    }
    
    if (length(rec_items) != 0) {
      res$precision <- res$TP/length(rec_items)
    }else{
      res$precision <- 1
    }
    
    if ((res$TP + res$FN) != 0){ 
      res$recall <- res$TP/(res$TP + res$FN)
    }else{
      res$recall <- 1
    }
    
    if ((res$precision + res$recall)!= 0 )
    {
      res$F1 = 2 * (res$precision * res$recall)/ (res$precision + res$recall)
    }
    
  }
  
  res
} 