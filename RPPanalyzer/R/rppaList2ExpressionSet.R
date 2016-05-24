`rppaList2ExpressionSet` <-
function(x) {
if (length(x)==3){
  ns <- x[[1]] ## normalized signal
  pd <- x[[2]] ## phenodata
  fd <- x[[3]] ## featuredata
  ## assayNames must equal fd
  rownames(fd) <- rownames(ns) <- 1:nrow(ns)
  } else {
         ns <- x[[1]] ## normalized signal
         pd <- x[[3]] ## phenodata
          fd <- x[[4]] ## featuredata
          ## assayNames must equal fd
          rownames(fd) <- rownames(ns) <- 1:nrow(ns)
          }
          
  ## pData
  pD <- new("AnnotatedDataFrame",
            data=as.data.frame(t(pd)))
  ## fData
  fD <- new("AnnotatedDataFrame",
               as.data.frame(fd))

  ## expression
  x <- new("ExpressionSet",
         exprs=ns,
         phenoData=pD,
         featureData=fD)
  return(x)
}

