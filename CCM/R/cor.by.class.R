cor.by.class <-
function(x,y, method = "pearson", use = "complete") {
  index = 1:length(y)
  index.by.class = split(index, y)
  K.list = list(NULL)
  for (i in 1:length(index.by.class)) {
    cat("calculating cor for: ", names(index.by.class)[i], "\n")
    K.list[[i]] = cor(x[,index.by.class[[i]]], method = method)
  }
  names = names(index.by.class)
  KK = list(NULL)
  for (i in 1:length(K.list)) {
    diag(K.list[[i]]) = NA
    K.list[[i]] = K.list[[i]][!is.na(K.list[[i]])]
    KK[[i]] = as.vector(K.list[[i]])
  }
  names(KK) = names
  return(K = KK)
}

