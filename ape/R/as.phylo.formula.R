## as.phylo.formula.R (2005-12-10)

##   Conversion from Taxonomy Variables to Phylogenetic Trees

## Copyright 2005 Julien Dutheil

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

as.phylo.formula <- function(x, data=parent.frame(), ...)
{
  # Testing formula syntax:
  err <- "Formula must be of the kind \"~A1/A2/.../An\"."
  if(length(x) != 2) stop(err)
  if(x[[1]] != "~") stop(err)
  f <- x[[2]]
  taxo <- list()
  while(length(f) == 3) {
    if(f[[1]] != "/") stop(err)
    if(!is.factor(data[[deparse(f[[3]])]])) stop(paste("Variable", deparse(f[[3]]), "must be a factor."))
    taxo[[deparse(f[[3]])]] <- data[[deparse(f[[3]])]]
    if(length(f) > 1) f <- f[[2]]
  }
  if(!is.factor(data[[deparse(f)]])) stop(paste("Variable", deparse(f), "must be a factor."))
  taxo[[deparse(f)]] <- data[[deparse(f)]]
  taxo.data <- as.data.frame(taxo)
  leaves.names <- as.character(taxo.data[,1])
  taxo.data[,1] <- 1:nrow(taxo.data)
  # Now builds the phylogeny:

  f.rec <- function(subtaxo) { # Recurrent utility function
    u <- ncol(subtaxo)
    levels <- unique(subtaxo[,u])
    if(u == 1) {
      if(length(levels) != nrow(subtaxo))
        warning("Error, leaves names are not unique.")
      return(as.character(subtaxo[,1]))
    }
    t <- character(length(levels))
    for(l in 1:length(levels)) {
      x <- f.rec(subtaxo[subtaxo[,u] == levels[l],][1:(u-1)])
      if(length(x) == 1) t[l] <- x
      else t[l] <- paste("(", paste(x, collapse=","), ")", sep="")
    }
    return(t)
  }
  string <- paste("(", paste(f.rec(taxo.data), collapse=","), ");", sep="")
  phy<-read.tree(text=string)
  phy$tip.label <- leaves.names[as.numeric(phy$tip.label)]
  return(phy)
}
