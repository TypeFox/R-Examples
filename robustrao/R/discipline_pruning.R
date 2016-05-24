#install.packages("igraph", repos="http://cran.rstudio.com/")

#library("igraph")

#' Set of permissible disciplines for redistribution.
#'
#' Computes the set of disciplines to which uncategorized references can be redistributed.
#' This set is computed taking into account the mutual similarities of the already referenced disciplines, as explained in Calatrava et al. (2016).
#' This function allows to set a tolerance of similarity that only permits similar disciplines to participate in the redistribution process.
#' Therefore, it avoids redistributions that include very dissimilar and improbable disciplines.
#'
#' @param r A logical vector indicating which disciplines are referenced by the current document.
#' Its length is equal to the total number of disciplines.
#' @param tolerance A real number in the interval [0,1].
#' This argument modulates the similarity between disciplines with which the strictness of the pruning of unlikely disciplines is controlled.
#' A value of 0 allows all disciplines to participate in the redistribution process.
#' A value of 1 permits no tolerance.
#' This argument is optional and leaving it unspecified deactivates tolerances.
#' @param similarity A positive semi-definite matrix that encodes the similarity between disciplines, as explained in Porter and Rafols (2009).
#' The dimensions of this matrix are \emph{n} x \emph{n}, being \emph{n} the total number of disciplines.
#' The number of rows and the number of columns of this matrix needs to be equal to the length of \code{r}.
#' The self-similarities (i.e. the diagonal elements) have to be 1.
#' @return A logical vector indicating to which disciplines a reference redistribution is permissible.
#' @examples
#' #Load data
#' data(pubdata1)
#'
#' #Get counts of citations of one of the publications in the dataset
#' counts <- pd1.count.matrix[,1]
#'
#' #Get logical vector indicating which disciplines are referenced by the publication
#' logic.disciplines <- counts > 0
#'
#' PruneDisciplines(logic.disciplines, 0.233, pd1.similarity)
#' @references
#' Calatrava Moreno, M. C., Auzinger, T. and Werthner, H. (2016) On the uncertainty of interdisciplinarity measurements due to incomplete bibliographic data. Scientometrics. DOI:10.1007/s11192-016-1842-4
#'
#' Porter, A. and Rafols, I. (2009) Is science becoming more interdisciplinary? Measuring and mapping six research fields over time. Scientometrics, Vol. 81, No. 3 (719-745). DOI:10.1007/s11192-008-2197-2
#' @import igraph
#' @export
PruneDisciplines <- function(r,
                             tolerance = 1,
                             similarity) {

  n <- length(r)
  # Error handling.
  if (n < 1 || n != nrow(similarity) || n != ncol(similarity)) {
    stop("Arguments 'r' and 'similarity' have incompatible sizes.")
  }
  if (is.nan(tolerance) || tolerance < 0 || tolerance > 1) {
    stop("Argument 'tolerance' is out of range.")
  }
  if (any(is.nan(similarity)) || any(similarity < 0) || any(similarity > 1)) {
    stop("Elements of 'similarity' are out of range.")
  }
  if (any(diag(similarity) != 1)) {
    stop("Elements of the diagonal of 'similarity' are not 1.")
  }

  # Check if at least one discipline is referenced.
  if (all(!r)) {
    return(!logical(length = n)) # All disciplines are permissible if none is referenced.
  }
  # Check for total tolerance.
  if (tolerance == 0) {
    return(!logical(length = n)) # All disciplines are permissible.
  }

  # Take the referenced subset of the similarity matrix.
  referenced.discipline.count <- sum(r)
  referenced.disciplines      <- which(r)
  unreferenced.disciplines    <- which(!r)
  referenced.part.similarity <- similarity[-unreferenced.disciplines, -unreferenced.disciplines, drop = FALSE]
  stopifnot(nrow(referenced.part.similarity) == referenced.discipline.count)

  # Compute the minimum spanning tree of the subset when interpreted as a weighted undirected graph.
  if (referenced.discipline.count == 1) {
    spanning.tree <- matrix(1, 1, 1) # 'minimum.spanning.tree' would yield matrix(0, 1, 1).
  } else {
    g <- graph.adjacency(1 - referenced.part.similarity, mode = "undirected", weighted = TRUE)
    spanning.tree <- minimum.spanning.tree(g)
    spanning.tree <- get.adjacency(spanning.tree, type = "both", sparse = FALSE)
  }

  # Get the minimal similarity between a disciplines and all connected ones.
  spanned.similarity <- referenced.part.similarity
  spanned.similarity[spanning.tree == 0] <- Inf
  min.similarity.per.discipline <- apply(spanned.similarity, 1, min)
  stopifnot(length(min.similarity.per.discipline) == referenced.discipline.count)

  # Mark disciplines that are sufficiently similar.
  permissible.disciplines <- logical(length = n) # Initialize as all FALSE.
  for (k in 1:referenced.discipline.count) {
    all.similarities.per.discipline <- similarity[referenced.disciplines[k], ]
    permissible.disciplines[all.similarities.per.discipline >= tolerance * min.similarity.per.discipline[k]] <- TRUE
  }

  return(permissible.disciplines)
}
