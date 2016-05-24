##' Taxonomic method of phyndr that works with generalised sets of
##' taxonomic information.  Requires a nested set of taxonomic classes
##' (e.g, genus, family, order, etc) but does not assume that these
##' classes are necessarily monophyletic.  \code{phyndr_genus} does
##' this for the genus level only but attempts to automatically detect
##' genera from the tree (assuming that tips are all genus/species
##' pairs, with genus and species separated by either an underscore or
##' space).
##'
##' The algorithm (not including recursion and substitute whatever
##' taxonomic level for genus).
##'
##' 1: drop genera if there is a species match for that genus, but
##' don't drop the actual matches.
##'
##' 2: work out which genera can be collapsed to a single tip due to
##' monophyly.
##'
##' 2a. For each genus, determine if they are monophyletic
##'
##' 2b. Then, for genera that are not monophyletic determine which can
##' be made monophyletic by dropping groups that are not represented
##' in the data set.
##'
##' 2c. Collapse each genus down, dropping only genera that are
##' required to achieve monophyly.
##'
##' 3. For each collapsed node we'll mangle the name to "genus::name"
##' and then make a list of suitable species in the tree and data set;
##' for now stored as \code{clades}.
##' @title Phyndr taxonomic
##' @param phy An ape phylogeny
##' @param data_species A vector of species names for which we have
##' trait data.  Species names in both the tree and in this vector
##' must be separated with underscores, not with spaces.
##' @param taxonomy A data.frame with taxonomic information.  Row
##' names must be present and must list every species in \code{phy}
##' and every species in \code{data_species}.  One or more columns
##' must be present; the first column is the lowest (finest) taxonomic
##' grouping and the last column is the highest (coarsest) taxonomic
##' grouping.  The names are arbitrary but will be used in creating
##' mangled names in the resulting phylogeny.
##'
##' @importFrom stats setNames
##' @export
phyndr_taxonomy <- function(phy, data_species, taxonomy) {
  if (!is.ultrametric(phy)) {
    stop("phy must be ultrametric")
  }
  if (any(is.na(taxonomy))) {
    stop("Missing values in taxonomy")
  }

  ## Might be worth doing an initial round of cleaning here.
  ## We only need species in the lookup that are in the taxonomy and
  ## in the tree.
  to_drop <- setdiff(phy$tip.label, union(data_species, rownames(taxonomy)))
  phy <- drop_tip(phy, to_drop)

  to_drop <- setdiff(rownames(taxonomy), union(data_species, phy$tip.label))
  taxonomy <- taxonomy[setdiff(rownames(taxonomy), to_drop), , drop=FALSE]

  to_drop <- setdiff(data_species, union(phy$tip.label, rownames(taxonomy)))
  data_species <- setdiff(data_species, to_drop)

  ## TODO: check
  ##   - taxonomy is a data.frame
  ##   - has row labels
  ##   - has unique column labels
  ##   - has at least one column
  ##   - has a tree structure?

  ## This is the recursive exit condition:
  if (ncol(taxonomy) < 1L || all(phy$tip.label %in% data_species)) {
    return(phyndr_taxonomy_cleanup(phy, data_species))
  }

  msg <- setdiff(phy$tip.label, rownames(taxonomy))
  if (length(msg) > 0L) {
    extra <- taxonomy[msg, , drop=FALSE]
    rownames(extra) <- msg
    for (i in names(extra)) {
      extra[[i]] <- sprintf("Unknown-%s-%s", i, msg)
    }
    taxonomy <- rbind(taxonomy, extra)
  }

  phy_g <- taxonomy[phy$tip.label, 1, drop=TRUE]
  dat_g <- taxonomy[data_species,  1, drop=TRUE]

  ## I don't think we want to run the whole way down the taxonomy
  ## straight away here.
  ##
  ## If there are unmatched species in the same genus as things that
  ## have data already, those get discarded from the tree as there's
  ## nothing that we can do.  We do keep other genera though, even if
  ## they're in a family that we have matches for so that implies that
  ## we don't roll back more than the first taxonomic grouping.
  match_s <- phy$tip.label %in% data_species

  ## On the second way around here we have to be more careful; most
  ## tips are actually ok by now.
  ## This one is wrong because it's against the wrong tree.  Here we
  ## need to know what group things belong to.  That means updating
  ## all the book-keeping so that's a pain.  Easiest way would be to
  ## append rows to the table and rematch.
  match_g <- phy_g %in% unique(phy_g[match_s])
  to_drop <- match_g & !match_s

  phy2 <- drop_tip(phy, phy$tip.label[to_drop])

  ## These are groups in the tree which did not match any species:
  phy_g_msg <- unique(phy_g[!match_g])

  ## Of these, these are the ones that have no data, so we don't want to
  ## go crazy with them.
  phy_g_msg_nodata <- setdiff(phy_g_msg, dat_g)
  ## and ones with data:
  phy_g_msg_data <- intersect(phy_g_msg, dat_g)

  ## Test missing genera to detemine which are monophyletic:
  ## TODO: This is really slow on very large trees, so that's not
  ## great.
  phy_g_msg_data_is_mono <-
    vlapply(phy_g_msg_data, is_monophyletic_group, phy_g, phy)
  phy_g_msg_data_mono <- phy_g_msg_data[phy_g_msg_data_is_mono]

  ## Then, of the ones that aren't, which can be fixed?
  check <- phy_g_msg_data[!phy_g_msg_data_is_mono]
  problems <- lapply(check, find_paraphyletic, phy_g, phy)
  can_fix <- vlapply(problems, function(x) all(x %in% phy_g_msg_nodata))
  phy_g_msg_data_fixable <- check[can_fix]

  ## Then, of the things with *no data* and which we don't drop above:
  ## * Genera that have no data but are in the tree (phy_g_msg_nodata)
  ## * not things that are going to be removed because of they make
  ##   other things paraphyletic.
  phy_g_msg_nodata_is_mono <- logical(length(phy_g_msg_nodata))
  names(phy_g_msg_nodata_is_mono) <- phy_g_msg_nodata
  check <- setdiff(phy_g_msg_nodata,
                   unlist(problems[can_fix]))
  phy_g_msg_nodata_is_mono[check] <-
    vlapply(check, is_monophyletic_group, phy_g, phy)
  phy_g_msg_nodata_mono <- phy_g_msg_nodata[phy_g_msg_nodata_is_mono]

  g_collapse <- c(phy_g_msg_data_mono,
                  phy_g_msg_data_fixable,
                  phy_g_msg_nodata_mono)
  g_drop <- unname(unlist(problems[can_fix]))

  tmp <- split(phy$tip.label, phy_g)
  to_drop <- c(unlist(lapply(tmp[g_collapse],
                             function(x) x[-1]), use.names=FALSE),
               unlist(tmp[g_drop], use.names=FALSE))

  relabel <- vcapply(tmp[g_collapse], function(x) x[[1]])
  names(relabel) <- sprintf("%s::%s", names(taxonomy)[[1]], names(relabel))

  taxonomy_extra <- taxonomy[match(g_collapse, taxonomy[[1]]), -1, drop=FALSE]
  rownames(taxonomy_extra) <- names(relabel)
  if (nrow(taxonomy_extra) > 0L) {
    taxonomy2 <- rbind(taxonomy[, -1, drop=FALSE], taxonomy_extra)
  } else {
    taxonomy2 <- taxonomy[, -1, drop=FALSE]
  }

  phy2 <- drop_tip(phy2, to_drop)
  phy2$tip.label[match(relabel, phy2$tip.label)] <- names(relabel)

  ## Now, assemble the list of plausible matches:
  tmp <- split(data_species, dat_g)
  clades <- setNames(tmp[g_collapse], names(relabel))

  data_species2 <- c(data_species,
                     names(clades)[viapply(clades, length) > 0L])

  phy2$clades <- c(phy2$clades, clades)
  phy2$clades <- phy2$clades[names(phy2$clades) %in% phy2$tip.label]

  if (!inherits(phy2, "phyndr")) {
    class(phy2) <- c("phyndr", class(phy2))
  }

  if (ncol(taxonomy_extra) == 0) {
    ## All done!
    phyndr_taxonomy_cleanup(phy2, data_species2)
  } else {
    ## Let's recurse!
    phyndr_taxonomy(phy2, data_species2, taxonomy2)
  }
}

##' @rdname phyndr_taxonomy
##' @export
phyndr_genus <- function(phy, data_species) {
  spp <- union(phy$tip.label, data_species)
  taxonomy <- data.frame(genus=split_genus(spp))
  rownames(taxonomy) <- spp
  phyndr_taxonomy(phy, data_species, taxonomy)
}

phyndr_taxonomy_cleanup <- function(phy, data_species) {
  to_drop <- setdiff(phy$tip.label, data_species)
  if (length(to_drop) > (length(phy$tip.label) - 2L)) {
    keep <- intersect(phy$tip.label, data_species)
    if (length(keep) == 0L) {
      stop("Zero overlap in data_species and phy")
    } else {
      stop("Only one species/clade in tree: ", keep)
    }
  }
  phy$clades <- phy$clades[names(phy$clades) %in% data_species]
  drop_tip(phy, to_drop)
}
