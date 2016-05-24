# SVD####
setMethod("predict", signature = c(model = "SVDclass"), function(model, Round = FALSE) {
  
  
  item_not_rated <- which(model@data@data == 0)
  
  # generate predictions
  p <- model@factors$U %*% t(model@factors$V) 
  
  model@data@data[item_not_rated] <- p[item_not_rated]
  
  roundData(model@data, Round)
})

# IBKNN####
setMethod("predict", signature = c(model = "IBclass"), function(model, Round = FALSE) {
  
  data <- model@data
  
  model_row <- nrow(data@data)
  item_names <- colnames(data@data)
  user_names <- rownames(data@data)
  
  rated_items <- which(data@data != 0)
  
  item_not_rated <- vector("list", model_row)
  
  # determining the set of items to be predicted
  item_not_rated <- lapply(1:model_row, function(i) which(data@data[i, ] == 0))
  item_not_rated <- lapply(1:model_row, function(i) unname(item_not_rated[[i]]))
  
  
  
  # generating predictions for the requested items. Weighted sum computation.
  
  globaMean <- sum(data@data)/numRatings(model@data) 
  
  for (i in 1:model_row) {
    for (item in item_not_rated[[i]]) {
      sim_item_ratings <- data@data[i, model@sim_index_kNN[item, ]]
      item_sim <- model@sim[item, model@sim_index_kNN[item, ]] 
      
      sim_sum <- sum(item_sim * (sim_item_ratings != 0), na.rm = TRUE )
      
      denom_w_sum <- sum(abs(item_sim) * (sim_item_ratings != 0), na.rm = TRUE)
      # sim_sum <- sum(abs(item_sim))
      
      if (sim_sum == 0) data@data[i, item] <- globaMean
      
      if (denom_w_sum == 0) {
        data@data[i, item] <- globaMean
      } else {
        data@data[i, item] <- sum(sim_item_ratings * item_sim , na.rm = TRUE)/denom_w_sum
      }
    }
  }
  
  roundData(data, Round)
  
})


# ALS####
setMethod("predict", signature = c(model = "wALSclass"), function(model, Round = FALSE) {
  
  item_not_rated <- which(model@data@data == 0)
  
  # generate predictions
  p <- model@factors$U %*% t(model@factors$V)
  p <- p * model@weightScheme
  
  # replacing not rated items
  model@data@data[item_not_rated] <- p[item_not_rated]
  
  roundData(model@data, Round)
  
})
# bpr####
setMethod("predict", signature = c(model = "BPRclass"), function(model, Round = FALSE) {
  
  item_not_rated <- which(model@data@data == 0)
  
  # generate predictions
  p <- model@factors$U %*% t(model@factors$V)
  
  model@data@data[item_not_rated] <- p[item_not_rated]
  
  roundData(model@data, Round)
})

# average####

setMethod("predict", signature = c(model = "algAverageClass"), function(model, Round = FALSE) {
  
  item_not_rated <- which(model@data@data == 0)
  
  model@data@data[item_not_rated] <- model@average[item_not_rated]
  
  
  roundData(model@data, Round)
}) 




#roundData############
roundData <- function(data, Round){
  if(Round){
    if (!data@binary & data@halfStar)
    {
      data@data <- round(data@data * 2)/2
    }else{
      data@data <- round(data@data)
    }
    data@data[data@data < data@minimum] <- 0
    data@data[data@data > data@maximum] <- data@maximum
  }
  data@data
}




