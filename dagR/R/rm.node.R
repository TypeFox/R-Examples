rm.node <-
function(dag, node)
{ # this conveniently removes a node;
  # also removes all paths, if an arc is removed;
  # searchType and searchRes are generally removed,
  #  for simplicity;
    dag$cov.types <- dag$cov.types[-node]
    dag$x <- dag$x[-node]
    dag$y <- dag$y[-node]
    for (i in nrow(dag$arc):1) {
        if ((dag$arc[i, 1] == node) || (dag$arc[i, 2] == node)) {
            dag <- rm.arc(dag, i)
        }
    }
    dag$arc <- matrix(sapply(X = dag$arc, FUN = function(x) {
        if (x > node) x - 1 else x }), byrow = FALSE, ncol = 2)
    dag$names <- dag$names[-node]
    dag$adj <- dag$adj[dag$adj != node]
    dag$searchType <- NULL
    dag$searchRes <- NULL
    return(dag)
}


