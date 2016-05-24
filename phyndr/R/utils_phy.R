split_genus <- function(str) {
  str_split <- strsplit(str, "[_ ]+")
  vcapply(str_split, "[[", 1L)
}

## Ape's version is too slow; this is around 10x faster.
is_monophyletic <- function(phy, tips) {
  if (length(tips) == 1L) {
    TRUE
  } else {
    tips_i <- match(tips, phy$tip.label)
    mrca <- mrca_tipset(phy, tips_i)
    desc <- get_descendants(mrca, phy, tips.only=TRUE)
    all(desc %in% tips_i)
  }
}

is_monophyletic_group <- function(group, table, phy) {
  is_monophyletic(phy, phy$tip.label[table == group])
}

## Find the things that make a group non-monophyletic
find_paraphyletic <- function(group, table, phy) {
  tips <- phy$tip.label[table == group]
  if (length(tips) == 1L) {
    character(0)
  } else {
    tips_i <- match(tips, phy$tip.label)
    mrca <- mrca_tipset(phy, tips_i)
    desc <- get_descendants(mrca, phy, tips.only=TRUE)
    setdiff(table[desc], group)
  }
}

## Find the largest clade in phy that includes the tip species `tip`
## but does not include any species listed in the vector `exclude`.
find_exclusive_clade <- function(tip, exclude, phy) {
  i <- match(tip, phy$tip.label)
  if (is.na(i)) {
    character(0)
  } else {
    ret <- list(species=tip, node=i)
    exclude <- setdiff(exclude, tip)

    repeat {
      i_parent <- phy$edge[match(i, phy$edge[, 2]), 1]
      if (is.na(i_parent)) {
        ## Hit the root:
        return(ret)
      }
      desc <- get_descendants_names(i_parent, phy)
      if (any(desc %in% exclude)) {
        return(ret)
      }

      i <- i_parent
      ret <- list(species=desc, node=i)
    }
  }
}
