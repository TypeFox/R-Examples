##' outliers
##' 
##' see function
##' 
##' 
##' @keywords internal
##' @export outliers
"outliers" <-  function(data, labels, threshold) {
  model <- train(data,labels)
  dist <- distance(data, model, labels, metric="mahal")
  good <- NULL
  bad  <- NULL
  ddd  <- NULL
  for(n in 1:length(model$label)){
    l <- model$label[n]
    which <- labels==l
    dd <- dist[which, n]
    meandd <- mean(dd)
    sdevdd <- sqrt(var(dd))
    
    cutoff <- sdevdd*threshold
    
    idx <- (1:length(labels))[which]
    tmp <- idx[dd<cutoff]
    good <- c(good, tmp)
    tmp <- idx[dd>=cutoff]
    bad  <- c(bad, tmp )
    ddd <- c(ddd, dist[which,n] )
  }
  list(good=good[!is.na(good)], bad=bad[!is.na(bad)], dist=ddd)
}
