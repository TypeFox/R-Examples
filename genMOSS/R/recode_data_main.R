recode_data_main <-
function (data, dimens, alpha) {
  
  recoded_data <- array (0, c(dim(data)[1], dim(data)[2]))
  for (j in 1:(dim(data)[2]-1)) {
    v <- as.factor(data[,j])
    if (length(levels(v)) > dimens[j]) {
      stop(paste("The dimens vector does not agree with the data. For example, dimens[",j,"] must be increased to at least ", length(levels(v)), ".", sep = ""))
    }
    recoded_data[,j] <- as.double(v) - 1
  }

  v <- as.factor(data[,dim(data)[2]])
  if (length(levels(v)) != 2)
    stop ("Response must be binary")
  recoded_data[,dim(data)[2]] <- as.double(v) - 1
  
  recoded_data <- as.data.frame(recoded_data)
  colnames(recoded_data) <- colnames(data)  

  recoded_dimens <- dimens
  response <- recoded_data[,dim(data)[2]]

  for (i in 1:(dim(data)[2]-1)) {
    if (dimens[i] == 3) {
      candidates <- array (0, c(dim(data)[1],4))
      candidates[,1] <- data[,i]
      candidates[,2] <- ifelse (recoded_data[,i] == 0, 0, 1)
      candidates[,3] <- ifelse (recoded_data[,i] == 1, 0, 1)
      candidates[,4] <- ifelse (recoded_data[,i] == 2, 0, 1)
      temp <- optimal_coding (as.data.frame(cbind(candidates, response)), dimens = c(3,2,2,2,2), alpha = alpha)
      recoded_data[,i] <- temp[[1]]
      recoded_dimens[i] <- temp[[2]]
    }
  }
  
  return(list(recoded_data = recoded_data, recoded_dimens = recoded_dimens))

}
