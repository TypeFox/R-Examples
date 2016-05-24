node.age <- function(tree,tp=-Inf,t0=NULL) {
  if (! is.binary.tree(tree)) tree <- multi2di(tree)
  node.ids <- unique(tree$edge[,1])
  root <- node.ids[! node.ids %in% tree$edge[,2]]
  if (length(root) != 1) stop("Can't identify root edge.")
  if (is.null(t0)) {
    node.ages <- get.age(tree,root,0.0)
    t0 <- max(node.ages[,3])
    tp <- tp + t0
  }
  get.age(tree,root,0.0,tp=tp)
}

