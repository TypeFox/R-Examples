##' Add a semi-colon to end of tree string
##'
##' Check if tree string ends in semi-colon and if not add one.  This
##' is mostly done for compatibility with ape, which requires them.
##' @param x A character string or vector of character strings each
##' representing a tree in Newick format.
##' @return The same value, but with a semi-colon added to the end
##' of any strings which did not already end in semi-colons.
##' @export
##' @example inst/examples/fix-semicolon-tree.R
##' @author Melissa J. Hubisz
fix.semicolon.tree <- function(x) {
  n <- nchar(x)
  lastCh <- substr(x, n, n)
  f <- lastCh != ";"
  if (sum(f) > 0L)
    x[f] <- paste(x[f], ";", sep="")
  x
}
    


##' Read a tree from a file
##'
##' Reads a tree in newick format
##' @title Read a Newick Tree from a File
##' @param filename The file containing the tree.
##' @return a character string representing the tree in newick format
##' @keywords trees newick
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
##' @example inst/examples/read-newick-tree.R
read.newick.tree <- function(filename) {
  check.arg(filename, "filename", "character", null.OK=FALSE)
  tr <- .Call.rphast("rph_tree_read", filename)
  fix.semicolon.tree(tr)
}


##' Get the number of nodes in a tree
##' @title Number of Nodes in a Tree
##' @param tree A vector of character strings, each containing a newick tree
##' @return A numeric vector containing the number of nodes in each tree
##' @keywords trees
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
##' @example inst/examples/numnodes-tree.R
numnodes.tree <- function(tree) {
  check.arg(tree, "tree", "character", null.OK=FALSE, min.length=1,
            max.length=NULL)
  result <- integer(length(tree))
  for (i in 1:length(tree)) {
    result[i] <- .Call.rphast("rph_tree_numnodes", tree[i])
  }
  result
}


##' Get the number of leaves in a tree
##' @title Number of leaves in a Tree
##' @param tree A vector of character strings, each containing a newick tree
##' @return A numeric vector containing the number of leaves (species) in each tree
##' @keywords trees
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
numleaf.tree <- function(tree) {
  check.arg(tree, "tree", "character", null.OK=FALSE, min.length=1,
            max.length=NULL)
  result <- integer(length(tree))
  for (i in 1:length(tree)) {
    result[i] <- .Call.rphast("rph_tree_numnodes", tree[i])
  }
  (result+1)/2
}

##' Get the total length of the edges of a tree
##' @param tree A vector of character strings, each containing a newick tree
##' @return A numeric vector containing the total branchlength of each tree
##' @keywords trees
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
branchlength.tree <- function(tree) {
  check.arg(tree, "tree", "character", null.OK=FALSE, min.length=1,
            max.length=NULL)
  result <- numeric(length(tree))
  for (i in 1:length(tree))
    result[i] <- .Call.rphast("rph_tree_branchlen", tree[i])
  result
}


##' Get the distance from a node to the root of a tree
##' @param tree A vector of character strings, each containing a newick tree
##' @param node A vector of character strings, giving the node name to
##' use for each tree.  Will be recycled to the length of the first argument.
##' @return A numeric vector containing the distance from each given
##' node to the root of the corresponding tree.
##' @keywords trees
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
depth.tree <- function(tree, node) {
  check.arg(tree, "tree", "character", null.OK=FALSE, min.length=1,
            max.length=NULL)
  node <- rep(node, length.out = length(tree))
  result <- numeric(length(tree))
  for (i in 1:length(tree))
    result[i] <- .Call.rphast("rph_tree_depth", tree[i], node[i])
  result
}



##' Prune sequences from a file
##' @title Prune a Tree
##' @param tree A vector of character strings, each containing a newick tree
##' @param seqs The sequences to prune from the trees
##' @param all.but A logical value.  If false, prunes all the named sequences
##' from the tree.  If TRUE, prunes all sequences except the ones named.
##' @return a vector of character strings representing the pruned trees.
##' @keywords trees
##' @author Melissa J. Hubisz and Adam Siepel
##' @example inst/examples/prune-tree.R
##' @export
prune.tree <- function(tree, seqs, all.but=FALSE) {
  check.arg(tree, "tree", "character", null.OK=FALSE,
            min.length=1, max.length=NULL)
  check.arg(seqs, "seqs", "character", null.OK=FALSE, min.length=1,
            max.length=NULL)
  check.arg(all.but, "all.but", "logical", null.OK=FALSE)
  result <- character(length(tree))
  for (i in 1:length(tree)) {
    result[i] <- .Call.rphast("rph_tree_prune", tree[i], seqs, all.but)
  }
  fix.semicolon.tree(result)
}


##' Name ancestors of a tree
##' @title Name Ancestral Nodes
##' @param tree A vector of character strings, each containing a newick tree
##' @return A vector of character strings containing newick trees with all
##' ancestors named.
##' @keywords trees
##' @author Melissa J. Hubisz and Adam Siepel
##' @example inst/examples/name-ancestors.R
##' @export
name.ancestors <- function(tree) {
  check.arg(tree, "tree", "character", null.OK=FALSE,
            min.length=1, max.length=NULL)
  result <- character(length(tree))
  for (i in 1:length(tree)) {
    result[i] <- .Call.rphast("rph_tree_name_ancestors", tree[i])
  }
  fix.semicolon.tree(result)
}


##' Get a subtree
##' @title Subtree
##' @param tree A vector of character strings, each containing a newick tree
##' @param node A vector of character strings, each representing the name
##' of the node which will be the new root of the tree.  If node is shorter
##' than tree, values will be recycled, and a warning produced if \code{length(tree) \%\% length(node) != 0}
##' @param super.tree A vector of logical values.  If TRUE, then remove all
##' nodes which are descendants of node, rather than keeping them.
##' @return A vector of trees which have been pruned, removing all nodes
##' which are not descendants of the given node.
##' @keywords trees
##' @author Melissa J. Hubisz and Adam Siepel
##' @example inst/examples/subtree.R
##' @export
subtree <- function(tree, node, super.tree=FALSE) {
  check.arg(tree, "tree", "character", null.OK=FALSE,
            min.length=1, max.length=NULL)
  check.arg(node, "node", "character", null.OK=FALSE,
            min.length=1, max.length=length(tree))
  check.arg(super.tree, "super.tree", "logical", null.OK=FALSE,
            min.length=1, max.length=length(tree))
  if (length(tree) %% length(node) != 0)
    warning("number of trees is not multiple of number of given nodes")
  node <- rep(node, length.out = length(tree))
  super.tree <- rep(super.tree, length.out = length(tree))
  result <- character(length(tree))
  for (i in 1:length(tree)) {
    if (super.tree[i]) {
      result[i] <- .Call.rphast("rph_tree_supertree", tree[i], node[i])
    } else result[i] <- .Call.rphast("rph_tree_subtree", tree[i], node[i])
  }
  fix.semicolon.tree(result)
}



##' Rescale a tree
##' @title Scale a Tree or Subtree
##' @param tree A vector of character strings, each containing a newick tree
##' @param scale A vector of scale factors for each tree (will be recycled
##' as necessary if shorter than trees)
##' @param subtree If not NULL, scaling will be on subtree defined by the
##' named node.  Subtrees will be recycled as necessary if shorter than trees.
##' @param include.leading (Only applicable when subtree used) If \code{TRUE},
##' include the branch leading to the named node in the subtree.
##' @return A vector of trees whose branches have been scaled
##' @keywords trees
##' @author Melissa J. Hubisz and Adam Siepel
##' @example inst/examples/rescale-tree.R
##' @export
rescale.tree <- function(tree, scale, subtree=NULL, include.leading=FALSE) {
  check.arg(tree, "tree", "character", null.OK=FALSE,
            min.length=1, max.length=NULL)
  check.arg(scale, "scale", "numeric", null.OK=FALSE,
            min.length=1, max.length=length(tree))
  check.arg(subtree, "subtree", "character", null.OK=TRUE,
            min.length=1, max.length=length(tree))
  if (!is.null(subtree))
    check.arg(include.leading, "include.leading", "logical", null.OK=FALSE,
              min.length=1, max.length=length(tree))
  if (length(tree) %% length(scale) != 0)
    warning("number of trees is not multiple of number of given scales")
  if ((!is.null(subtree)) && length(tree) %% length(subtree) != 0)
    warning("number of trees is not multiple of number of given subtrees")
  if ((!is.null(subtree)) && length(tree) %% length(include.leading) != 0)
    warning("number of trees is not multiple of length of include.leading")
  if (is.null(subtree)) subtreeVal <-  NULL else subtreeIdx <- 1
  scale <- rep(scale, length.out=length(tree))
  if (!is.null(subtree)) {
    subtree <- rep(subtree, length.out=length(tree))
    include.leading <- rep(include.leading, length.out=length(tree))
  }
  result <- character(length(tree))
  for (i in 1:length(tree)) {
    if (!is.null(subtree)) {
      subtreeVal <- subtree[i]
      includeVal <- include.leading[i]
    } else includeVal <- NULL
    result[i] <- .Call.rphast("rph_tree_scale", tree[i], scale[i], subtreeVal,
                              includeVal)
  }
  fix.semicolon.tree(result)
}

##' Rename nodes of trees
##' @title Tree Node Renaming
##' @param tree A vector of character strings, each containing a newick tree
##' @param old.names A vector of current names to be substituted
##' @param new.names A vector of equal length to old.names giving the
##' substitutions
##' @return A vector of character strings, in which all nodes with names
##' given in old.names are replaced with values from new.names
##' @keywords trees
##' @author Melissa J. Hubisz and Adam Siepel
##' @example inst/examples/rename-tree.R
##' @export
rename.tree <- function(tree, old.names, new.names) {
  check.arg(tree, "tree", "character", null.OK=FALSE,
            min.length=1, max.length=NULL)
  check.arg(old.names, "old.names", "character", null.OK=FALSE,
            min.length=1, max.length=NULL)
  check.arg(new.names, "new.names", "character", null.OK=FALSE,
            min.length=length(old.names), max.length=length(old.names))
  fix.semicolon.tree(.Call.rphast("rph_tree_rename", tree, old.names, new.names))
}


##' Apply a label to some branches
##' @title Label tree branches
##' @param tree A vector of character strings, each containing a newick tree
##' @param branches A vector of character strings, indicating the branches
##' which should get the label.  The branch is named by the node which
##' descends from it.  If multiple trees and branches are given, all named
##' branches will be labelled in every tree.
##' @param label A single character string giving the label to apply
##' to the named branches.
##' @return A vector of character strings containing the modified trees;
##' the branches are labelled by appending a pound sign and the label
##' to the node name in the tree string.
##' @keywords trees
##' @author  Melissa J. Hubisz
##' @example inst/examples/label-branches.R
##' @export
label.branches <- function(tree, branches, label) {
  check.arg(tree, "tree", "character", null.OK=FALSE,
            min.length=1L, max.length=NULL)
  check.arg(branches, "branches", "character", null.OK=FALSE,
            min.length=1L, max.length=NULL)
  check.arg(label, "label", "character", null.OK=FALSE)
  .Call.rphast("rph_tree_label_branches", tree, branches, label)
}


##' Apply a label to a subtree
##' @title Label subtree
##' @param tree A vector of character strings, each containing a newick tree
##' @param node A character string, giving the node at the head of the subtree.
##' @param label A single character string giving the label to apply
##' to the branches in the subtree.
##' @param include.leading A logical value; if \code{TRUE}, include the
##' branch leading to the subtree in the labelled group; otherwise include
##' only descendants of the named node.
##' @return A vector of character strings containing the modified trees;
##' the branches are labelled by appending a pound sign and the label
##' to the node name in the tree string.
##' @keywords trees
##' @author  Melissa J. Hubisz
##' @example inst/examples/label-subtree.R
##' @export
label.subtree <- function(tree, node, label, include.leading=FALSE) {
  check.arg(tree, "tree", "character", null.OK=FALSE,
            min.length=1L, max.length=NULL)
  check.arg(node, "node", "character", null.OK=FALSE)
  check.arg(label, "label", "character", null.OK=FALSE)
  check.arg(include.leading, "include.leading", "logical", null.OK=FALSE)
  .Call.rphast("rph_tree_label_subtree", tree, node, include.leading, label)
}


##' Get a summary of a Newick-formatted tree, edge lengths, node names, etc
##' @param object A character string containing a newick tree
##' @param ... Not currently used (exists for S3 compatibility)
##' @return A data frame with a row for every node, containing columns:
##' branch length (tparent), distance to root (troot), name, label (if tree labels
##' present), and parent, rchild, lchild.
##' @method summary tree
##' @export summary.tree
##' @export
##' @example inst/examples/summary-tree.R
##' @author Melissa J. Hubisz and Adam Siepel
summary.tree <- function(object, ...) {
  check.arg(object, "object", "character", null.OK=FALSE)
  tree <- object
  names <- .Call.rphast("rph_tree_summary_nodenames", tree)
  t <- .Call.rphast("rph_tree_summary_len", tree)
  if (sum(t < 0) >= 1L) t[t < 0] <- NA
  troot <- .Call.rphast("rph_tree_summary_depth", tree)
  parent <- .Call.rphast("rph_tree_summary_parent", tree)
  if (sum(parent < 0) >= 1L) parent[parent < 0] <- NA
  lchild <- .Call.rphast("rph_tree_summary_lchild", tree)
  if (sum(lchild < 0) >= 1L) lchild[lchild < 0] <- NA
  rchild <- .Call.rphast("rph_tree_summary_rchild", tree)
  if (sum(rchild < 0) >= 1L) rchild[rchild < 0] <- NA
  rv <- data.frame(name=names, tparent=t, troot=troot, parent=parent, lchild=lchild, rchild=rchild)
  label <- .Call.rphast("rph_tree_summary_label", tree)
  if (!is.null(label)) {
    rv <- data.frame(rv, label=label)
  }
  rv
}


##' Get the names of a tree's leaf nodes
##' @param object A character string containing a newick tree
##' @param ... Not currently used
##' @return A character vector containing the names of the leaf nodes
##' @export
##' @author Melissa J. Hubisz and Adam Siepel
leafnames.tree <- function(object, ...) {
  x <- summary.tree(object, ...)
  as.character(x[is.na(x$lchild),"name"])
}
