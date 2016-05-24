topographicError <- function (sommap) {
  norm.data <- preprocessData(sommap$data, sommap$parameters$scaling)
  norm.proto <- preprocessProto(sommap$prototypes, sommap$parameters$scaling,
                                sommap$data)
  if (sommap$parameters$type=="numeric") {
    all.dist <- apply(norm.data, 1, function(x) {
      apply(norm.proto,1,function(y) sum((x-y)^2) )
    })
    ind.winner2 <- apply(all.dist,2,function(x) order(x)[2])
  } else if (sommap$parameters$type=="korresp") {

    nr <- nrow(sommap$data)
    nc <- ncol(sommap$data)
    all.dist.row <- apply(norm.data[1:nr,1:nc], 1, function(x) {
      apply(sommap$prototypes[,1:nc], 1, function(y) sum((x-y)^2) )
    })
    all.dist.col <- apply(norm.data[(nr+1):(nr+nc),(nc+1):(nr+nc)],
                          1, function(x) {
      apply(sommap$prototypes[,(nc+1):(nr+nc)], 1, function(y) sum((x-y)^2) )
    })
    ind.winner2 <- c(apply(all.dist.col,2,function(x) order(x)[2]),
                     apply(all.dist.row,2,function(x) order(x)[2]))
  } else if (sommap$parameters$type=="relational") {
    dist.1 <- -0.5*diag(norm.proto%*%norm.data%*%t(norm.proto))
    dist.2 <- norm.proto%*%norm.data
    all.dist <- sweep(dist.2, 1, dist.1, "+")
    ind.winner2 <- apply(all.dist, 2, function(x) order(x)[2])
  }
  res.error <- mean(!sapply(1:nrow(sommap$data), function(x) {
      is.element(ind.winner2[x], selectNei(sommap$clustering[x],
                                           sommap$parameters$the.grid, 1,
                                           radius.type="letremy",
                                           dist.type="maximum"))
  }))
  return(res.error)
}

quantizationError <- function(sommap) {
  norm.data <- preprocessData(sommap$data, sommap$parameters$scaling)
  norm.proto <- preprocessProto(sommap$prototypes, sommap$parameters$scaling,
                                sommap$data)
  if (sommap$parameters$type=="numeric") {
    quantization.error <- sum(apply((norm.data-
                                       norm.proto[sommap$clustering,])^2,
                                    1,sum))/nrow(norm.data)
  } else if (sommap$parameters$type=="korresp") {

    nr <- nrow(sommap$data)
    nc <- ncol(sommap$data)
    quantization.error <- sum(apply((norm.data[1:nr,1:nc]-
                                       sommap$prototypes[sommap$clustering[
                                         (nc+1):(nc+nr)],1:nc])^2,1,sum))
    quantization.error <- quantization.error +
      sum(apply((norm.data[(nr+1):(nr+nc),(nc+1):(nr+nc)]-
                   sommap$prototypes[sommap$clustering[1:nc],(nc+1):(nr+nc)])^2,
                1,sum))
    quantization.error <- quantization.error / (nc+nr)
  } else if (sommap$parameters$type=="relational") {
    clust.proto <- norm.proto[sommap$clustering,]
    quantization.error <- clust.proto%*%norm.data - 0.5*
      tcrossprod(diag(clust.proto%*%norm.data%*%t(clust.proto)),
                 rep(1,ncol(norm.data)))
    quantization.error <- sum(diag(quantization.error))/nrow(norm.data)
  }
    
  quantization.error
}

# main function
quality.somRes <- function(sommap, quality.type=c("all", "quantization",
                                                  "topographic"), ...) {
  quality.type <- match.arg(quality.type)
  switch(quality.type,
         "all"=list("topographic"=topographicError(sommap),
                    "quantization"=quantizationError(sommap)),
         "topographic"=topographicError(sommap),
         "quantization"=quantizationError(sommap)
         )
}

quality <- function(sommap, quality.type,...) {
  UseMethod("quality")
}
