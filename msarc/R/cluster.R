################################################################################

# brain-damaged way to create a list with n copies of uniprot, named from "gos"
# In: list of N GO terms, uniprot ID
# Out: list of N copies of uniprot ID, named with the incoming GO terms
augmentGOlist <- function(gos,uniprot) {
  z <- lapply(gos,function(x) uniprot)
  names(z) <- gos
  return(z)
}

# In: named list: uniprot IDs to lists of GO terms
# Out: named list: GO terms to lists of Uniprot IDs
invertGOlist <- function(gos) {
  pile <- mapply(augmentGOlist,gos,names(gos))
  flat <- list()
  for (rmap in pile) {
    for (n in names(rmap)) {
      flat[[n]] = c(rmap[[n]],flat[[n]])
    }
  }
  return(flat)
}

# In: a list of GO terms.
# Out: the list of those terms plus all ancestors.
addOneAncestorSet <- function(terms) {
  res <- unique(unlist(c(terms,lapply(terms,function(x) GO.db::GOMFANCESTOR[[x]]))))
  return(res)
}

# ASSUMPTION!!!!! Go terms do not have ancestors in other hierarchies, e.g.
# a term in the "molecular function" tree will not have an ancestor in
# the "biological process" tree.  Must check this sometime...
# In: a list of lists of GO terms
# Out: the same list, but the sub-lists contain all ancestors of the original
#      set
addGoAncestors <- function(golist) {
  result <- lapply(golist,addOneAncestorSet)
  return(result)
}

# filter list of GO terms on GO tree name ("MF", "BP", "CC")
# In: list of GO objects (from eg2GO
# Out: list of GO ids.
filterGo <- function(goset,gocat) {
  flag <- sapply(goset,function(x) x$Ontology == gocat)
  filt <- goset[flag]
  goid <- sapply(filt,function(x) x$GOID)
  return(unique(goid))
}

# given a vector of Uniprot IDs and a bimap, return the IDs that are found
# as keys in the bimap.
filterUniprotInXref <- function(uniprots,bimap) {
  found <- uniprots %in% AnnotationDbi::mappedkeys(bimap)
  nf <- uniprots[!found]
  if (length(nf) > 0) {
    warning("missing Uniprot IDs: ",paste(nf,sep=",",collapse=","),call.=F)
  }
  return(uniprots[found])
}

# In: a single Uniprot ID,plus the two maps
# Out: a list of GO ids associated with this Uniprot ID
oneUniprot2go <- function(id,uniprot2eg,eg2go) {
  egmap = eg2go[uniprot2eg[[id]]]
  egs = AnnotationDbi::mappedkeys(egmap)
  gos = lapply(egs,function(x) eg2go[[x]])
  fgos = unlist(gos,recursive=F)
  return(fgos)
}

# note that "map" is a table that maps entrez gene IDs to Uniprot ids, from
# the "org.Hs.eg.db" package, or equivalent from another organism's package
#
# "gocat" is "MF", "BP", "CC"
#
# "uniprot" is a column of uniprot IDs.
#
# "org.Hs.egGO" converts Entrez IDs to GO terms
uniprot2go <- function(uniprot,eg2uniprot,eg2go,gocat) {
  uniprot2eg <- AnnotationDbi::revmap(eg2uniprot)
  uniprotFilt <- filterUniprotInXref(uniprot,uniprot2eg)
  gos = lapply(uniprotFilt,function(x) oneUniprot2go(x,uniprot2eg,eg2go))
  gfilt <- lapply(gos,filterGo,gocat)
  names(gfilt) <- uniprotFilt
  return(gfilt)
}

# In: list of lists of Uniprot IDs, named with GO terms
# Out: data frame of GO terms, counts of each category, GO descriptions
go2table <- function(gos) {
  ids <- names(gos)
  sizes <- sapply(gos,length)
  terms <- sapply(ids,function(x) GO.db::GOTERM[[x]]@Term)
  df <- data.frame(sizes,ids,terms,stringsAsFactors=F)
  colnames(df) = c('count','id','description')
  return(df)
}

# return T if a GO term's ancestor is in the list of candidates, F otherwise
parentInList <- function(id,candidates) {
  if (is.na(id)) {
    warning("Got NA in parentInList")
    return(F)
  }
  anc <- GO.db::GOMFANCESTOR[[id]]
  found = candidates %in% anc
  return(sum(found) > 0)
}

# return list of go terms which don't have ancestors in this list
topLevel <- function(gos) {
  tl <- sapply(gos,function(x) !parentInList(x,gos))
  return(tl)
}

# assign kids to parents.  If there are multiple parents for a child, pick
# the one with the smallest count...
assignKidsToParents <- function(parents,kids,counts) {
}

# return largest element of the list
listmax <- function(lst) {
  m <- lst[[1]]
  for (n in lst) {
    if (n > m) {
      m <- n
    }
  }
  return(m)
}

# get smallest element of ns, based on counts in "counts"
getSmallest <- function(ns,counts) {
  if (length(ns) == 1) {
    return(ns[[1]])
  }
  curcount <- listmax(counts)+1
  curid <- "zork"
  for (n in ns) {
    if (counts[[n]] < curcount) {
      curcount <- counts[[n]]
      curid <- n
    }
  }
  return(curid)
}

# In: a) an msarc object
#     b) eg2uniprot -- an AnnotationDBI map from entrez IDs to Uniprot IDs
#     c) eg2go -- an AnnotationDBI map from entrez IDs to GO terms
#     d) gocat -- which GO category to use (MF)
# Out: named list of GO terms, each entry being a list of included uniprot IDs.
msarc.findGOterms <- function(msarc,eg2uniprot=org.Hs.eg.db::org.Hs.egUNIPROT,eg2go=org.Hs.eg.db::org.Hs.egGO,minCount=10) {
  if (requireNamespace("GO.db",quietly=TRUE)) {
    gocat <- "MF"
    gos <- uniprot2go(msarc$data$uniprot,eg2uniprot,eg2go,gocat)
    msarc$go2uniBase <- invertGOlist(gos)
    nonzeroflags <- sapply(gos,function(x) length(x) > 0)
    nonzeroterms <- names(gos)[nonzeroflags]
    missing = msarc$data$uniprot[!(msarc$data$uniprot %in% nonzeroterms)]
    # augment with ancestors
    gosAll <- addGoAncestors(gos)
    # invert the list: go terms to uniprot IDs
    go2uni <- invertGOlist(gosAll)
    go2uni$all <- c(go2uni$all,missing)
    msarc$go2uni <- go2uni
    msarc$go2uniAll <- go2uni
    gotbl <- go2table(msarc$go2uni)
    colnames(gotbl) <- c('count','id','description')
    rownames(gotbl) <- gotbl$id
    msarc$gotbl <- gotbl[gotbl$count >= minCount,]
    return(msarc)
  } else {
    stop("Failed to load required package 'GO.db'.")
  }
}

msarc.getTerms <- function(msarc) {
  return(msarc$gotbl)
}

msarc.saveTerms <- function(msarc,filename='go_terms.txt') {
  write.table(msarc$gotbl,file=filename,quote=F,row.names=F,col.names=F,sep='\t')
}

setCandidates <- function(msarc) {
  tbl <- msarc$gotbl
  candidates <- as.list(tbl$id)
  counts <- as.list(tbl$count)
  names(counts) <- candidates
  msarc$go2uni <- msarc$go2uni[tbl$id]
  msarc$candidates <- candidates
  msarc$counts <- counts
  return(msarc)
}

msarc.loadTerms <- function(msarc,filename='go_terms.txt') {
  tbl = read.table(filename,header=F,stringsAsFactors=F,sep='\t')
  colnames(tbl) <- c('count','id','description')
  rownames(tbl) <- tbl$id
  msarc$gotbl <- tbl
  msarc <- setCandidates(msarc)
  return(msarc)
}

msarc.filterTerms <- function(msarc,keepers) {
  if (class(keepers) == "data.frame") {
    klist <- keepers$id
  } else {
    klist <- keepers
  }
  msarc$gotbl <- msarc$gotbl[klist,]
  msarc$gotbl <- msarc$gotbl[!is.na(msarc$gotbl$id),]
  msarc <- setCandidates(msarc)
  return(msarc)
}

makeHierarchy <- function(terms,counts) {
  if (length(terms) == 0) {
    return(terms)
  }
  tl = topLevel(terms)
  parents = terms[tl]
  kids = terms[!tl]
  kidsets = list()
  for (p in parents) {
    kidsets[[p]] = list()
  }
  for (k in kids) {
    candidates = GO.db::GOMFANCESTOR[[k]]
    ps = parents[parents %in% candidates]
    p = getSmallest(ps,counts)
    kidsets[[p]] <- c(kidsets[[p]],k)
  }
  tree <- list()
  for (p in parents) {
    tree[[p]] = makeHierarchy(kidsets[[p]],counts)
  }
  return(tree)
}
