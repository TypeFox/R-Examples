# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


ldclassify = function (data, means, covariance, template = NULL, posterior = 'no'){
  
  if (class(template) == 'template'){
    covariance = template$covariance
    means = template$means
  }
  
  data = as.matrix(data)
  means = as.matrix(means)
  covariance = as.matrix(covariance)
  
  distances = sapply (1:nrow(data), function (i){
    tmp = matrix(rep(data[i,], nrow(means)),nrow(means),ncol(means),byrow = TRUE)
    d = diag((tmp-means) %*% solve(covariance)%*% t(tmp-means))
  })
  winner = sapply (1:nrow(data), function (i){
    tmp = order(distances[,i])[1]
  })
  if (!is.null(rownames (means))) winner = as.factor (rownames(means)[winner])
  
  if (posterior=='winner'){
    tmppost = sapply (1:nrow(data), function (i){
      tmp = exp(-sort(distances[,i])[1]/2) / sum (exp(-distances[,i]/2))
    })
    winner = data.frame (winner, tmppost)
  }  
  if (posterior=='all'){
    tmppost = sapply (1:nrow(data), function (i){
      tmp = exp(-distances[,i]/2) / sum (exp(-distances[,i]/2))
    })
    winner = data.frame (winner, t(tmppost))
    colnames (winner)[-1] = labels
  }  
  return (winner)  
}

