##' phyndr for topo trees
##' @title phyndr for topo trees
##' @param phy An ape phylogeny.  The tree must be ultrametric and
##' this is enforced on entry.
##' @param data_species A vector of species names for which we have
##' trait data.  Species names in both the tree and in this vector
##' must be separated with underscores, not with spaces.
##' @param topology A topological tree used to determine patterns of
##' relatedness among species in \code{phy}.  The better resolved this
##' tree is and the better the overlap it has with \code{phy} and
##' \code{data_species} the more useful this will be.  Branch length
##' information is not used.
##'
##' @importFrom stats reorder setNames
##' @export
phyndr_topology <- function(phy, data_species, topology) {
  if (!is.ultrametric(phy)) {
    stop("phy must be ultrametric")
  }
  ## First, discard all species that are in the time tree but not in
  ## the topo tree and don't have data as these cannot be used.  I
  ## feel we could discard more than this (think of the case of a
  ## single data point) but perhaps that can't be done.
  to_drop <- setdiff(phy$tip.label, union(data_species, topology$tip.label))
  phy <- drop_tip(phy, to_drop)

  to_drop <- setdiff(topology$tip.label, union(data_species, phy$tip.label))
  topology <- drop_tip(topology, to_drop)

  ## Organise the time tree into postorder because we'll need that:
  phy <- reorder(phy, "postorder")

  candidates <- vector("list", max(phy$edge))
  desc <- vector("list", max(phy$edge))
  complete <- logical(length(candidates))

  ## Species that have data have themselves as possibilities:
  phy_data <- intersect(data_species, phy$tip.label)
  candidates[match(phy_data, phy$tip.label)] <- phy_data
  complete[match(phy_data, phy$tip.label)] <- TRUE

  desc[seq_along(phy$tip.label)] <- phy$tip.label

  ## Species without data are more complicated; we need to find the
  ## biggest clade in the topo tree that includes the tip that does
  ## not include any other species in the time tree.
  phy_nodata <- setdiff(phy$tip.label, data_species)

  tmp <- lapply(phy_nodata, find_exclusive_clade, phy$tip.label, topology)
  candidates[match(phy_nodata, phy$tip.label)] <-
    lapply(tmp, "[[", "species")

  nds <- unique(phy$edge[, 1])

  for (nd in nds) {
    nd_d <- phy$edge[phy$edge[, 1] == nd, 2]

    if (any(complete[nd_d])) {
      complete[[nd]] <- TRUE
    } else {
      desc[[nd]] <- unlist(desc[nd_d])
      ## Is it possible that we'll only hit a single species in the
      ## tree?  I hope not!  Might be worth checking though.
      topo_mrca <- mrca_tipset(topology, match(desc[[nd]], topology$tip.label))
      if (is.null(topo_mrca)) {
        ## Unlikely to be triggered...
        stop("Error reconciling trees")
      }
      topo_desc <- get_descendants_names(topo_mrca, topology)
      if (any(complete[match(topo_desc, phy$tip.label)], na.rm=TRUE)) {
        complete[[nd]] <- TRUE
      } else {
        candidates[[nd]] <- unlist(candidates[nd_d])
        candidates[nd_d] <- list(NULL)
      }
    }
  }

  candidates_data <- lapply(candidates, intersect, data_species)

  nc <- viapply(candidates_data, length)
  is_tip <- seq_along(nc) <= length(phy$tip.label)

  to_drop_tips <- phy$tip.label[nc == 0L & is_tip]

  ## Sort out the internal nodes:
  nodes <- nc > 1L & !is_tip
  nodes_keep <- vcapply(desc[nodes], "[[", 1L)
  nodes_drop <- unlist(lapply(desc[nodes], "[", -1L))
  nodes_cand <- setNames(candidates[nodes], nodes_keep)

  tips_cand <- setNames(candidates[is_tip], phy$tip.label)

  to_drop <- setdiff(union(nodes_drop, to_drop_tips), nodes_keep)

  phy2 <- drop_tip(phy, to_drop)

  clades <- tips_cand
  clades[nodes_keep] <- nodes_cand
  clades <- clades[phy2$tip.label]

  ## hmm.  I think that the n=1 cases here should just be straight up
  ## relabellings.
  is_single <- viapply(clades, length) == 1L
  clades_single <- unlist(clades[is_single])

  relabel <- clades_single[clades_single != names(clades_single)]
  if (length(relabel) > 0L) {
    phy2$tip.label[match(names(relabel), phy2$tip.label)] <- unname(relabel)
  }

  clades <- clades[!is_single]

  phy2$clades <- clades
  if (!inherits(phy2, "phyndr")) {
    class(phy2) <- c("phyndr", class(phy2))
  }
  phy2
}
