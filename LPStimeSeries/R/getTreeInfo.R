getTreeInfo <- function(object, which.tree=1) {
  if (is.null(object$forest)) {
    stop("No forest component in ", deparse(substitute(object)))
  }
  if (which.tree > object$ntree) {
    stop("There are fewer than ", which.tree, "trees in the forest")
  }
  k <- which.tree
  if (object$type == "regression") {
      tree <- cbind(object$forest$leftDaughter[,k],
                    object$forest$rightDaughter[,k],
                    object$forest$bestvar[,k],
                    object$forest$splitType[,k],
                    object$forest$xbestsplit[,k],
                    object$forest$nodestatus[,k],
                    object$forest$nodedepth[,k],
                    object$forest$nodepred[,k])[1:object$forest$ndbigtree[k],]
  } 
  
  dimnames(tree) <- list(1:nrow(tree), c("left daughter", "right daughter",
                                         "split segment","split type", "split point",
                                         "status", "depth","prediction"))
  out <- list(segment.length=object$segment.length[k],target=object$target[k],
				target.type=object$target.type[k],tree=tree)
  
  return(out)
}
