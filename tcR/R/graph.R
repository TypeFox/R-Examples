######### MUTATION NETWORK MANAGING ##########


#' Make mutation network for the given repertoire.
#' 
#' @description
#' Mutation network (or a mutation graph) is a graph with vertices representing nucleotide or in-frame amino acid sequences (out-of-frame amino acid sequences
#' will automatically filtered out) and edges are connecting pairs of sequences with hamming distance or edit distance between them
#' no more than specified in the \code{.max.errors} function parameter.
#' 
#' @param .data Either character vector of sequences, data frame with \code{.label.col}
#' or shared repertoire (result from the \code{shared.repertoire} function) constructed based on \code{.label.col}.
#' @param .method Either "hamm" (for hamming distance) or "lev" (for edit distance). Passed to the \code{find.similar.sequences} function.
#' @param .max.errors Passed to the \code{find.similar.sequences} function.
#' @param .label.col Name of the column with CDR3 sequences (vertex labels).
#' @param .seg.col Name of the column with V gene segments.
#' @param .prob.col Name of the column with clonotype probability.
#' 
#' @return Mutation network, i.e. igraph object with input sequences as vertices labels, ???
#' 
#' @seealso \link{shared.repertoire}, \link{find.similar.sequences}, \link{set.people.vector}, \link{get.people.names}
#' 
#' @examples
#' \dontrun{
#' data(twb)
#' twb.shared <- shared.repertoire(twb)
#' G <- mutation.network(twb.shared)
#' get.people.names(G, 300, T)  # "Subj.A|Subj.B"
#' get.people.names(G, 300, F)  # list(c("Subj.A", "Subj.B"))
#' }
mutation.network <- function (.data, .method = c('hamm', 'lev'), .max.errors = 1,
                                   .label.col = 'CDR3.amino.acid.sequence', .seg.col = 'V.gene', .prob.col = 'Probability') {
  # Make vertices and edges.
  if (has.class(.data, 'character')) {
    .data <- data.frame(A = .data, stringsAsFactors = F)
    colnames(.data) <- .label.col
  }
  G <- graph.empty(n = nrow(.data), directed=F)
  G <- add.edges(G, t(find.similar.sequences(.data[[.label.col]], .method = .method[1], .max.errors = .max.errors)))
  G <- simplify(G)
  
  # Every label is a sequence.
  G <- set.vertex.attribute(G, 'label', V(G), .data[[.label.col]])
  
  # Add V-segments to vertices.
  if (.seg.col %in% colnames(.data)) {
    G <- set.vertex.attribute(G, 'vseg', V(G), .data[[.seg.col]])
  } else {
    G <- set.vertex.attribute(G, 'vseg', V(G), 'nosegment')
  }
  
  # Set sequences' indices as in the given shared repertoire.
  G <- set.vertex.attribute(G, 'repind', V(G), 1:vcount(G))
  
  # Set probabilities for sequences.
  if (.prob.col %in% colnames(.data)) {
    G <- set.vertex.attribute(G, 'prob', V(G), .data[[.prob.col]])
  } else {
    G <- set.vertex.attribute(G, 'prob', V(G), rep.int(-1, vcount(G)))
  }
  
  # Add people to vertices.
  if ('People' %in% colnames(.data)) {
    attr(G, 'people') <- colnames(shared.matrix(.data))
    G <- set.people.vector(G, .data)
  } else {
    attr(G, 'people') <- "Individual"
    G <- set.vertex.attribute(G, 'people', V(G), 1)
    G <- set.vertex.attribute(G, 'npeople', V(G), 1)
  }
  
  G
}


#' Set and get attributes of a mutation network related to source people.
#' 
#' @aliases set.people.vector get.people.names
#' 
#' @description
#' Set vertice attributes 'people' and 'npeople' for every vertex in the given graph.
#' Attribute 'people' is a binary string indicating in which repertoire sequence are
#' found. Attribute 'npeople' is a integer indicating number of repertoires, in which
#' this sequence has been found.
#' 
#' @usage
#' set.people.vector(.G, .shared.rep)
#' 
#' get.people.names(.G, .V = V(.G), .paste = T)
#' 
#' @param .G Mutation network.
#' @param .shared.rep Shared repertoire.
#' @param .V Indices of vertices.
#' @param .paste If TRUE than concatenate people names to one string, else get a character vector of names.
#' 
#' @return New graph with 'people' and 'npeople' vertex attributes or character vector of length .V or list of length .V.
#' 
#' @examples
#' \dontrun{
#' data(twb)
#' twb.shared <- shared.repertoire(twb)
#' G <- mutation.network(twb.shared)
#' get.people.names(G, 300, T)  # "Subj.A|Subj.B"
#' get.people.names(G, 300, F)  # list(c("Subj.A", "Subj.B"))
#' }
set.people.vector <- function (.G, .shared.rep) {
  .shared.rep[is.na(.shared.rep)] <- 0
  .G <- set.vertex.attribute(.G, 'people', V(.G),
                             apply(as.matrix(.shared.rep[, -(1:(match('People', colnames(.shared.rep))))]),
                                   1,
                                   function (row) { paste0(as.integer(row > 0), collapse='') }))
  .G <- set.vertex.attribute(.G, 'npeople', V(.G),
                             apply(as.matrix(.shared.rep[, -(1:(match('People', colnames(.shared.rep))))]),
                                   1,
                                   function (row) { sum(row > 0) }))
}

get.people.names <- function (.G, .V = V(.G), .paste = T) {
  ppl <- attr(.G, 'people')
  if (!.paste) {
    lapply(strsplit(get.vertex.attribute(.G, 'people', .V), '', fixed=T, useBytes=T), function (l) {
      ppl[l == '1']
    })
  } else {
    sapply(strsplit(get.vertex.attribute(.G, 'people', .V), '', fixed=T, useBytes=T), function (l) {
      paste0(ppl[l == '1'], collapse='|')
    }, USE.NAMES = F)
  }
} 


#' Set group attribute for vertices of a mutation network
#' 
#' @aliases set.group.vector get.group.names
#' 
#' @description
#' asdasd
#' 
#' @usage
#' set.group.vector(.G, .attr.name, .groups)
#' 
#' get.group.names(.G, .attr.name, .V = V(.G), .paste = T)
#' 
#' @param .G Mutation network.
#' @param .attr.name Name of the new vertex attribute.
#' @param .V Indices of vertices.
#' @param .groups List with integer vector with indices of subjects for each group.
#' @param .paste if T then return character string with concatenated group names, else return list with character vectors
#' with group names.
#' 
#' @return igraph object with new vertex attribute \code{.attr.name} with binary strings for \code{set.group.vector}.
#' Return character vector for \code{get.group.names}.
#' 
#' @examples
#' \dontrun{
#' data(twb)
#' twb.shared <- shared.repertoire(twb)
#' G <- mutation.network(twb.shared)
#' G <- set.group.vector(G, "twins", list(A = c(1,2), B = c(3,4)))  # <= refactor this
#' get.group.names(G, "twins", 1)       # "A|B"
#' get.group.names(G, "twins", 300)     # "A"
#' get.group.names(G, "twins", 1, F)    # list(c("A", "B"))
#' get.group.names(G, "twins", 300, F)  # list(c("A"))
#' # Because we have only two groups, we can assign more readable attribute.
#' V(G)$twin.names <- get.group.names(G, "twins")
#' V(G)$twin.names[1]  # "A|B"
#' V(G)$twin.names[300]  # "A"
#' }
set.group.vector <- function (.G, .attr.name, .groups) {  
  d <- get.vertex.attribute(.G, 'people', V(.G))
  d <- do.call(rbind, lapply(strsplit(d, '', T, useBytes = T), as.integer))
  
  group.vec <- rep('nogroup', times = max(unlist(.groups)))
  for (gr.i in 1:length(.groups)) {
    for (elem in .groups[[gr.i]]) {
      group.vec[elem] <- names(.groups)[gr.i]
    }
  }
  grnames <- names(.groups)
  .G <- set.vertex.attribute(.G, .attr.name, V(.G),
                            apply(d, 1, function (row) { 
                              paste0(as.integer(grnames %in% sort(group.vec[row > 0])), collapse='')
                              }))
  
  attr(.G, .attr.name) <- names(.groups)
  .G
}

get.group.names <- function (.G, .attr.name, .V = V(.G), .paste = T) {
  grs <- attr(.G, .attr.name)
  if (!.paste) {
    lapply(strsplit(get.vertex.attribute(.G, .attr.name, .V), '', fixed=T, useBytes=T), function (l) {
      grs[l == '1']
    })
  } else {
    sapply(strsplit(get.vertex.attribute(.G, .attr.name, .V), '', fixed=T, useBytes=T), function (l) {
      paste0(grs[l == '1'], collapse='|')
    }, USE.NAMES = F)
  }
}


#' Get vertex neighbours.
#' 
#' @description
#' Get all properties of neighbour vertices in a mutation network of specific vertices.
#' 
#' @param .G Mutation network.
#' @param .V Indices of vertices for which return neighbours.
#' @param .order Neighbours of which order return.
#' 
#' @return List of length \code{.V} with data frames with vertex properties. First row in each data frame
#' is the vertex for which neighbours was returned.
#' 
#' @examples
#' \dontrun{
#' data(twb)
#' twb.shared <- shared.repertoire(twb)
#' G <- mutation.network(twb.shared)
#' head(mutated.neighbours(G, 1)[[1]])
#' #           label             vseg repind prob people npeople
#' # 1 CASSDRDTGELFF          TRBV6-4      1   -1   1111       4
#' # 2 CASSDSDTGELFF          TRBV6-4     69   -1   1100       2
#' # 3 CASSYRDTGELFF TRBV6-3, TRBV6-2    315   -1   1001       2
#' # 4 CASKDRDTGELFF TRBV6-3, TRBV6-2   2584   -1   0100       1
#' # 5 CASSDGDTGELFF          TRBV6-4   5653   -1   0010       1
#' # 6 CASSDRETGELFF          TRBV6-4   5950   -1   0100       1
#' }
mutated.neighbours <- function (.G, .V, .order = 1) {
  neis <- neighborhood(.G, .order, .V, mode = 'all')
  lapply(neis, function (l) { 
    res <- as.data.frame(lapply(list.vertex.attributes(.G), function (vattr) { 
      get.vertex.attribute(.G, vattr, l) } ),
      stringsAsFactors = F)
    colnames(res) <- list.vertex.attributes(.G)
    res
    } )
}


# srcppl.distribution <- function (.G) {  
#   one.count <- sapply(strsplit(V(.G)$people, '', fixed = T, useBytes = T), function (x) sum(x == '1'))
# #   mean(V(.G)$npeople)
# mean(one.count)
# }
# 
# 
# neippl.distribution <- function (.G, .exclude.zeros = F) {
#   if (.exclude.zeros) {
#     .G <- induced.subgraph(.G, degree(.G) != 0)
#   }
#   
#   one.count <- sapply(strsplit(V(.G)$people, '', fixed = T, useBytes = T), function (x) sum(x == '1'))
# #   one.count <- V(.G)$npeople
#   quantile(sapply(neighborhood(.G, 1), function (x) {
#     if (length(x) == 1) {
#       0
#     } else {
#       mean(one.count[x[-1]]) / (length(x) - 1)
#     }
#   }), prob = c(.025, .975))
# }
# 
# 
# pplvar.distribution <- function (.G, .exclude.zeros = F) {
#   if (.exclude.zeros) {
#     .G <- induced.subgraph(.G, degree(.G) != 0)
#   }
#   
#   ppl.inds <- lapply(strsplit(V(.G)$people, '', fixed = T, useBytes = T), function (x) which(x == '1'))
#   c1 <- quantile(sapply(neighborhood(.G, 1), function (x) {
#     if (length(x) == 1) {
#       0
#     } else {
# #       length(unique(unlist(ppl.inds[x[-1]]))) / length(x[-1])
#       length(unique( unlist(ppl.inds[x[-1]]) [ !(unlist(ppl.inds[x[-1]]) %in% unlist(ppl.inds[x[1]])) ] ))
#     }
#   }), prob = c(.25, .75))
# 
#   c2 <- mean(sapply(neighborhood(.G, 1), function (x) {
#     if (length(x) == 1) {
#       0
#     } else {
# #       length(unique(unlist(ppl.inds[x[-1]]))) / length(x[-1])
#       length(unique( unlist(ppl.inds[x[-1]]) [ !(unlist(ppl.inds[x[-1]]) %in% unlist(ppl.inds[x[1]])) ] ))
#     }
#   }))
# 
#   c(c1[1], Mean = c2, c1[2])
# }