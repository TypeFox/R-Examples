appnn <- function(sequences) UseMethod("appnn")

appnn.default <- function(sequences){
  strsqncs <- .aux$str2char(sequences)
  results <- list()
  for(i in 1:length(strsqncs)){
    results[[i]] <- list()
    ncdcsqnc <- .aux$encode(strsqncs[[i]],.aaindex,.mapping)
    rawdata <- .nn$classify(.nn, ncdcsqnc)
    results[[i]]$sequence <- sequences[i]
    results[[i]]$overall <- max(rawdata)
    results[[i]]$aminoacids <- .aux$aaprediction(rawdata)
    results[[i]]$hotspots <- .aux$hotspotsidentification(results[[i]]$aminoacids)
  }
  class(results) <- "appnn"
  return(results)
}

print.appnn <- function(x, ...){
  for(i in 1:length(x)){
    cat('\nPrediction: \n')
    print(x[[i]])
  }
}

plot.appnn <- function(x, indices, ...){
  par(mfrow=c(1,length(indices)))
  for(i in indices){
    plot(x[[i]]$aminoacids, xaxt = "n", main = x[[i]]$sequence, ylab="Propensity", xlab="Sequence",ylim=c(-0.5,1.5))
    axis(1, at=1:length(.aux$str2char(x[[i]]$sequence)[[1]]), labels=.aux$str2char(x[[i]]$sequence)[[1]])
    lines(x[[i]]$aminoacids)
    abline(h=0.5,col=2,lty=2)
  }
}
