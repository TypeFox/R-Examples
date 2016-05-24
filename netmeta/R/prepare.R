prepare <- function(TE, seTE, treat1, treat2, studlab, tau=0){
  
  weights <- 1/(seTE^2 + tau^2)
  
  data <- data.frame(studlab,
                     treat1, treat2,
                     treat1.pos=NA, treat2.pos=NA,
                     TE, seTE, weights,
                     narms=NA, stringsAsFactors=FALSE)
  ##
  ## Ordering data set
  ##
  o <- order(data$studlab, data$treat1, data$treat2)
  ##
  data <- data[o,]
  ##
  ## Adapt numbers to treatment IDs
  ##
  names.treat <- sort(unique(c(data$treat1, data$treat2)))
  data$treat1.pos <- match(data$treat1, names.treat)
  data$treat2.pos <- match(data$treat2, names.treat)
  
  newdata <- data[1,][-1,]
  ##
  sl <- unique(data$studlab)
  ##
  ## Determining number of arms and adjusting weights of
  ## multi-arm studies
  ##
  for (s in sl){
    subgraph <- data[data$studlab==s,]
    subgraph$narms <- (1+sqrt(8*dim(subgraph)[1]+1))/2
    if (dim(subgraph)[1] > 1)
      subgraph$weights <- 1/multiarm(1/subgraph$weights)$v ## Reciprocal new weights
    newdata <- rbind(newdata, subgraph)
  }
  res <- newdata
  ##
  res$order <- o
  ##
  res
}
