# $Id: prune.bagging.R,v 1.2 2002/09/12 08:59:13 hothorn Exp $

prune.classbagg <- function(tree, cp=0.01,...)
{
  for(i in 1:length(tree$mtrees))
    tree$mtrees[[i]]$btree <- prune( tree$mtrees[[i]]$btree, cp=cp, ...)
  tree
}

prune.regbagg <- function(tree, cp=0.01,...)
{
  for(i in 1:length(tree$mtrees))
    tree$mtrees[[i]]$btree <- prune( tree$mtrees[[i]]$btree, cp=cp, ...)
  tree
}


prune.survbagg <- function(tree, cp=0.01,...)
{
  for(i in 1:length(tree$mtrees))
    tree$mtrees[[i]]$btree <- prune( tree$mtrees[[i]]$btree, cp=cp, ...)
  tree
}
