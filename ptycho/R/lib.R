# Convenience functions to avoid code repetition in the package
createXNames <- function(p) paste("X",seq_len(p),sep="")
createYNames <- function(q) paste("y",seq_len(q),sep="")

createGroupsSim <- function(G, p) {
  var2group <- rep(seq_len(G), each=ceiling(p/G))[seq_len(p)]
  names(var2group) <- createXNames(p)
  group2var <- NULL
  for (ngrp in seq_len(G)) {
    group2var[[ngrp]] <- which(var2group == ngrp)
    names(group2var)[ngrp] <- paste("g",ngrp,sep="")
  }
  list(var2group=var2group, group2var=group2var, sizes=laply(group2var, length))
}

# Not user level because we cannot distribute the actual data.
# ARGUMENTS
#   X  design matrix: rows correspond to subjects and columns correspond to
#      variants.  Must have column names that are the same as the row names in
#      var.detail; they will be carried through the entire process, making
#      analysis easier.
#   var.detail
#      data frame: Row names same as column names of X; must have columns MAF
#      and GENE
#   MAF.threshold
#      All variants in a gene that have var.detail$MAF less than this value will
#      be grouped together.  All other variants will be in groups by themselves.
#   is.sim is logical
#      TRUE: When simulating phenotypes from actual genotype data, I identify
#      common variants as being in groups by themselves and rare variants as
#      being in groups of size larger than one.  All groups that would consist
#      of one rare variant are omitted from group2var in the output.
#      FALSE: All groups are included in group2var.
# VALUE
#   List as output by createGroupsSim() except that if is.sim is TRUE then
#   length(group2var) is strictly less than max(var2group) and some variants are
#   not in unlist(group2var).
createGroupsActualX <- function(X, var.detail, MAF.threshold, is.sim) {
  p <- ncol(X)
  if (length(intersect(colnames(X), rownames(var.detail))) != p) {
    stop("createGroups: colnames(X) must equal rownames(var.detail)")
  }
  var2group <- rep(0, p)
  names(var2group) <- colnames(X)
  var.cmn <- with(var.detail, var.detail[MAF >= MAF.threshold, ])
  var.grp <- with(var.detail,
                  var.detail[MAF < MAF.threshold & !is.na(GENE)
                             & annotation %in% c("missense","nonsense"), ])
  genes.grp <- ddply(droplevels(var.grp), c("GENE"), nrow)
  genes.grp <- with(genes.grp, genes.grp[V1 > 1, "GENE"])
  var.grp <- droplevels(with(var.grp, var.grp[GENE %in% genes.grp, ]))
  # First assign variants in var.grp to groups.
  ngenes <- length(levels(var.grp$GENE))
  X.ind <- match(rownames(var.grp), colnames(X))
  var2group[X.ind] <- var.grp$GENE
  # Next assign variants in var.cmn to groups by themselves.
  X.ind <- match(rownames(var.cmn), colnames(X))
  ncommon <- nrow(var.cmn)
  var2group[X.ind] <- ngenes + seq_len(ncommon)
  # Finally assign all remaining variants to groups by themselves.
  var.bad <- (var2group == 0)
  var2group[var.bad] <- ngenes + ncommon + seq_len(sum(var.bad))
  # Create group2var
  group2var <- NULL
  ngrp2var <- if (is.sim) ngenes+ncommon else max(var2group)
  for (ngrp in seq_len(ngrp2var)) {
    group2var[[ngrp]] <- which(var2group == ngrp)
    names(group2var)[ngrp] <- (if (ngrp <= ngenes) levels(var.grp$GENE)[ngrp]
                               else colnames(X)[group2var[[ngrp]]])
  }
  list(var2group=var2group, group2var=group2var, sizes=laply(group2var, length))
}

createOrthogonalX <- function(n, p) {
  if (n < p) stop("createOrthogonalX: cannot have p > n")
  eye <- diag(nrow=p)
  X <- matrix(rep(eye, ceiling(n/p)), ncol=p, byrow=TRUE)
  colnames(X) <- createXNames(p)
  scale(X[seq_len(n),], center=FALSE)
}
