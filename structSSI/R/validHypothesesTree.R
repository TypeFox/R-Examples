## We check to see that the tree that is being defined is valid.
## I.e., we check that the number of p-values, nodes in the tree,
## and hypotheses under consideration are consistent.

validHypothesesTree <- function(object){
    p.vals <- slot(object, "p.vals")
    tree <- slot(object, "tree")
    if(ncol(tree) != 2){
        return("Edgelist for '@tree' is invalid, does not two columns (parent -> child).")
    }
    if(nrow(p.vals) != nrow(tree)){
        return("Number of hypotheses in '@p.vals' does not match number of nodes in hypotheses tree in '@tree'")
    }
    TRUE
}
