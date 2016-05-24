#install.packages("gmp", repos="http://cran.rstudio.com/")
#install.packages("iterpc", repos="http://cran.rstudio.com/")
#install.packages("quadprog", repos="http://cran.rstudio.com/")

#library(gmp)
#library(iterpc)
#library(quadprog)

#' Lower bound of the uncertainty interval of the Rao-Stirling diversity index.
#'
#' This function computes the lower bound of the uncertainty interval of the Rao-Stirling diversity index, as explained in Calatrava et al. (2016).
#' The computation involves the redistribution of uncategorized references to various disciplines.
#' In order to avoid improbable redistributions of disciplines, a set of permissible disciplines for redistribution can be defined.
#' Furthermore, the number of disciplines redistributed to uncategorized references can be limited.
#'
#' @param known.ref.counts A vector of positive integers. Each element represents the count of references to each discipline.
#' @param uncat.ref.count A positive integer denoting the number of references that are not categorized into disciplines.
#' @param permissible.disciplines A logical vector denoting to which disciplines uncategorized references can be distributed.
#' Its length needs to be equal to the length of \code{known.ref.counts}.
#' This argument is optional and leaving it unspecified or supplying NULL permits redistribution to all disciplines.
#' @param redistribution.limit A positive integer that limits the number of disciplines that each uncategorized reference can be redistributed to.
#' This argument is optional and leaving it unspecified permits redistribution to all disciplines at once.
#' @param similarity A positive semi-definite matrix that encodes the similarity between disciplines, as explained in Porter and Rafols (2009).
#' The dimensions of this matrix are \emph{n} x \emph{n}, being \emph{n} the total number of disciplines.
#' The self-similarities (i.e. the diagonal elements) have to be 1.
#' @param max.batch.size A positive integer that sets the size of the batch of candidates that is computed at once.
#' This positive value determines the quantity of allocated memory and has to be reduced if corresponding errors arise.
#' This argument is optional and leaving it unspecified sets it to a default value.
#' @return The lower bound of the uncertainty interval of the Rao-Stirling diversity index.
#' @section Warning:
#' This function solves a computationally intensive optimization problem. In order to reduce the search space it is recommended to provide the function with the vector of permissible disciplines and redistribution limit.
#' When very dissimilar disciplines are referenced by the categorized references, a warning message is displayed to inform the user.
#' Such cases require longer computation times.
#' The dataset \code{\link{pubdata2}} contains an example of a publication that requires intensive computation in order to calculate the uncertainty interval of the Rao-Stirling diversity index.
#' @examples
#' ##EXAMPLE 1
#' #Load data
#' data(pubdata1)
#'
#' #Get counts of citations of one of the publications in the dataset
#' counts <- pd1.count.matrix[,1]
#'
#' #Get number of uncategorized references in the publication
#' uncat <- pd1.uncat.refs[1]
#'
#' #Get vector of permissible disciplines.
#' logic.disciplines <- counts > 0
#' permissible <- PruneDisciplines(logic.disciplines, 0.233, pd1.similarity)
#'
#' LowerIndexBound(counts, uncat, pd1.similarity, permissible)
#'
#' @references
#' Calatrava Moreno, M. C., Auzinger, T. and Werthner, H. (2016) On the uncertainty of interdisciplinarity measurements due to incomplete bibliographic data. Scientometrics. DOI:10.1007/s11192-016-1842-4
#'
#' Porter, A. and Rafols, I. (2009) Is science becoming more interdisciplinary? Measuring and mapping six research fields over time. Scientometrics, Vol. 81, No. 3 (719-745). DOI:10.1007/s11192-008-2197-2
#' @import gmp iterpc
#' @export
LowerIndexBound <- function(known.ref.counts,
                                   uncat.ref.count,
                                   similarity,
                                   permissible.disciplines = NULL,
                                   redistribution.limit = 4,
                                   max.batch.size = 131072) {

  n <- length(known.ref.counts)
  # Error handling.
  if (n < 1 || n != nrow(similarity) || n != ncol(similarity)) {
    stop("Arguments 'known.ref.counts' and 'similarity' have incompatible sizes.")
  }
  if (any(is.nan(known.ref.counts)) || any(known.ref.counts < 0)) {
    stop("Elements of 'known.ref.counts' are out of range.")
  }
  if (is.nan(uncat.ref.count) || uncat.ref.count < 0) {
    stop("Argument 'uncat.ref.count' is out of range.")
  }
  if (sum(known.ref.counts) == 0 && uncat.ref.count == 0) {
    stop("Elements of 'known.ref.counts' and argument 'uncat.ref.count' cannot simultaneously be 0.")
  }
  if (!is.null(permissible.disciplines)) {
    if (!is.logical(permissible.disciplines)) {
      stop("Argument 'permissible.disciplines' is not a logical vector.")
    } else if (length(permissible.disciplines) != n) {
      stop("Arguments 'permissible.disciplines' and 'known.ref.counts' have incompatible sizes.")
    } else if (any(known.ref.counts[!permissible.disciplines] > 0)) {
      stop("Argument 'known.ref.counts' must not contain positive values for disciplines that are excluded from redistribution by 'permissible.disciplines'.")
    }
  } else {
    permissible.disciplines <- !logical(length = n); # Allow all disciplines.
  }
  if (is.nan(redistribution.limit) || redistribution.limit < 0) {
    stop("Argument 'redistribution.limit' is out of range.")
  }
  if (any(is.nan(similarity)) || any(similarity < 0) || any(similarity > 1)) {
    stop("Elements of 'similarity' are out of range.")
  }
  if (any(diag(similarity) != 1)) {
    stop("Elements of the diagonal of 'similarity' are not 1.")
  }
  if (is.nan(max.batch.size) || max.batch.size < 1) {
    stop("Argument 'max.batch.size' is out of range.")
  }

  # Early out.
  if (uncat.ref.count == 0) {
    return (RaoCounts(known.ref.counts, similarity)) # No unknown references.
  }
  if (sum(permissible.disciplines) == 0) {
    warning("No disciplines are permissible for redistribution. No discipline can be assigned to the uncategorized references.")
    return (RaoCounts(known.ref.counts, similarity))
  }
  if (redistribution.limit == 0) {
    warning("The redistribution limit is set to 0 disciplines. Disciplines cannot be assigned to the uncategorized references.")
    return (RaoCounts(known.ref.counts, similarity))
  }

  # Compute the number of hypercube vertices that have to be checked.
  permissible.discipline.count <- sum(permissible.disciplines)
  candidate.vertex.count <- as.bigz(0)
  for (hypercube.vertex.norm in 1:min(n-1, redistribution.limit, permissible.discipline.count)) {
    # Norms 0 and 'n' do not need to be checked.
    # The 1-norm cannot exceed the number of dimensions, which are given by the permissible discipline count.
    # For a given 1-norm of the vertex distance, compute in how many ways ones can be chosen to sum to this amount.
    candidate.vertex.count <- candidate.vertex.count + chooseZ(permissible.discipline.count, hypercube.vertex.norm)
  }

  # Restrict to permissible disciplines.
  c   <- known.ref.counts[permissible.disciplines]
  sim <- similarity[permissible.disciplines, permissible.disciplines, drop = FALSE]
  u   <- uncat.ref.count
  k   <- redistribution.limit

  # Iterate over all permissible vertex 1-norms to identify the smallest possible index.
  min.index <- Inf
  large.candidate.count.reported <- FALSE
  for (hypercube.vertex.norm in 1:min(n-1, k, permissible.discipline.count)) {
    # Norms 0 and 'n' do not need to be checked.
    # The 1-norm cannot exceed the number of dimensions, which are given by the permissible discipline count.
    # Start with norm 1 to enable the early out at the end of the loop body.
    redistributed.reference.counts.norm <- sum(c) + hypercube.vertex.norm * u
    redistributed.reference.counts.norm.square <-
      redistributed.reference.counts.norm * redistributed.reference.counts.norm

    # For a given 1-norm of the vertex distance, enumerate the multiset permutations of 0 and 1 coordinates.
    vertex.iterator <- iterpc(c(permissible.discipline.count - hypercube.vertex.norm, hypercube.vertex.norm),
                              labels = c(0, u),
                              ordered = TRUE)
    remaining.vertex.count <- getlength(vertex.iterator)

    # Iterate over all vertices of the current 1-norm.
    while (remaining.vertex.count > 0) {
      # Perform iteration in batches.
      vertex.batch <- getnext(vertex.iterator, min(remaining.vertex.count, max.batch.size), drop = FALSE)

      # Compute index for the current vertex batch
      vertex.batch <- t(vertex.batch)
      redistributed.reference.counts.batch <- vertex.batch + c # Offset hypercube vertices by known reference count.
      bilinear.form.componentwise.unnormalized.batch <-
        redistributed.reference.counts.batch * (sim %*% redistributed.reference.counts.batch)
      bilinear.form.unnormalized.batch <- colSums(bilinear.form.componentwise.unnormalized.batch)
      bilinear.form.normalized.batch <- bilinear.form.unnormalized.batch / redistributed.reference.counts.norm.square
      min.index.batch <- 1 - max(bilinear.form.normalized.batch)
      min.index <- min(min.index, min.index.batch) # Update preliminary minimal index.
      remaining.vertex.count <- remaining.vertex.count - max.batch.size # Update remaining vertices.
    }
    # Early out
    if (min.index == 0) {
      break
    }
    # If the early out does not apply, report a potentially lengthy computation.
    # An error is thrown for computations that would take an exceedingly long time to compute.
    if (!large.candidate.count.reported) {
      if (candidate.vertex.count > 2^30) {
        stop("Too many candidates detected (", as.character(candidate.vertex.count), " candidates).")
      } else if (candidate.vertex.count > 2^20) {
        warning("Large number of candidates detected (", as.character(candidate.vertex.count), " candidates).",
                immediate. = TRUE)
      }
      large.candidate.count.reported <- TRUE
    }
  }

  return(min.index)
}


#' Upper bound of the uncertainty interval of the Rao-Stirling diversity index.
#'
#' This function computes the upper bound of the uncertainty interval of the Rao-Stirling diversity index, as explained in Calatrava et al. (2016).
#' The computation involves the redistribution of uncategorized references to various disciplines.
#' In order to avoid improbable redistributions of disciplines, a set of permissible disciplines for redistribution can be defined.
#' Furthermore, the number of disciplines redistributed to uncategorized references can be limited.
#'
#' @param known.ref.counts A vector of positive integers. Each element represents the count of references to each discipline.
#' @param uncat.ref.count A positive integer denoting the number of references that are not categorized into disciplines.
#' @param permissible.disciplines A logical vector denoting to which disciplines uncategorized references can be distributed.
#' Its length needs to be equal to the length of \code{known.ref.counts}.
#' This argument is optional and leaving it unspecified or supplying NULL permits redistribution to all disciplines.
#' @param redistribution.limit A positive integer that limits the number of disciplines that each uncategorized reference can have redistributed.
#' This argument is optional and leaving it unspecified will set the redistribution.limit to default.
#' @param similarity A positive semi-definite matrix that encodes the similarity between disciplines, as explained in Porter and Rafols (2009).
#' The dimensions of this matrix are \emph{n} x \emph{n}, being \emph{n} the total number of disciplines.
#' The self-similarities (i.e. the diagonal elements) have to be 1.
#' @return The upper bound of the uncertainty interval of the Rao-Stirling diversity index.
#' @examples
#' #Load data
#' data(pubdata1)
#'
#' #Get counts of citations of one of the publications in the dataset
#' counts <- pd1.count.matrix[,1]
#'
#' #Get number of uncategorized references in the publication
#' uncat <- pd1.uncat.refs[1]
#'
#' #Get vector of permissible disciplines.
#' logic.disciplines <- counts > 0
#' permissible <- PruneDisciplines(logic.disciplines, 0.233, pd1.similarity)
#'
#' UpperIndexBound(counts, uncat, pd1.similarity, permissible)
#'
#' @references
#' Calatrava Moreno, M.C., Auzinger, T. and Werthner, H. (2016) On the uncertainty of interdisciplinarity measurements due to incomplete bibliographic data. Scientometrics. DOI:10.1007/s11192-016-1842-4
#'
#' Porter, A. and Rafols, I. (2009) Is science becoming more interdisciplinary? Measuring and mapping six research fields over time. Scientometrics, Vol. 81, No. 3 (719-745). DOI:10.1007/s11192-008-2197-2
#' @import quadprog
#' @export
UpperIndexBound <- function(known.ref.counts,
                                   uncat.ref.count,
                                   similarity,
                                   permissible.disciplines = NULL,
                                   redistribution.limit = 4) {

  n <- length(known.ref.counts)
  # Error handling.
  if (n < 1 || n != nrow(similarity) || n != ncol(similarity)) {
    stop("Arguments 'known.ref.counts' and 'similarity' have incompatible sizes.")
  }
  if (any(is.nan(known.ref.counts)) || any(known.ref.counts < 0)) {
    stop("Elements of 'known.ref.counts' are out of range.")
  }
  if (is.nan(uncat.ref.count) || uncat.ref.count < 0) {
    stop("Argument 'uncat.ref.count' is out of range.")
  }
  if (sum(known.ref.counts) == 0 && uncat.ref.count == 0) {
    stop("Elements of 'known.ref.counts' and argument 'uncat.ref.count' cannot simultaneously be 0.")
  }
  if (!is.null(permissible.disciplines)) {
    if (!is.logical(permissible.disciplines)) {
      stop("Argument 'permissible.disciplines' is not a logical vector.")
    } else if (length(permissible.disciplines) != n) {
      stop("Arguments 'permissible.disciplines' and 'known.ref.counts' have incompatible sizes.")
    } else if (any(known.ref.counts[!permissible.disciplines] > 0)) {
      stop("Argument 'known.ref.counts' must not contain positive values for disciplines that are excluded from redistribution by 'permissible.disciplines'.")
    }
  } else {
    permissible.disciplines <- !logical(length = n); # Allow all disciplines.
  }
  if (is.nan(redistribution.limit) || redistribution.limit < 0) {
    stop("Argument 'redistribution.limit' is out of range.")
  }
  if (any(is.nan(similarity)) || any(similarity < 0) || any(similarity > 1)) {
    stop("Elements of 'similarity' are out of range.")
  }
  if (any(diag(similarity) != 1)) {
    stop("Elements of the diagonal of 'similarity' are not 1.")
  }

  # Early out.
  if (uncat.ref.count == 0) {
    return (RaoCounts(known.ref.counts, similarity)) # No unknown references.
  }
  if (sum(permissible.disciplines) == 0) {
    warning("No disciplines are permissible for redistribution. Disciplines cannot be assigned to uncategorized references.")
    return (RaoCounts(known.ref.counts, similarity))
  }
  if (redistribution.limit == 0) {
    warning("The redistribution limit is set to 0 disciplines. Disciplines cannot be assigned to uncategorized references.")
    return (RaoCounts(known.ref.counts, similarity))
  }

  # Restrict to permissible disciplines.
  c   <- known.ref.counts[permissible.disciplines]
  sim <- similarity[permissible.disciplines, permissible.disciplines, drop = FALSE]

  # Start building quadratic programming problem.
  lc <- length(c)
  u  <- uncat.ref.count
  k  <- redistribution.limit

  # Preallocate (in)equality constraint structures.
  maximal.constraint.count <- 1 + 2*lc + lc*lc
  bvec <- numeric(length = maximal.constraint.count)
  Amat <- matrix(data = 0, nrow = maximal.constraint.count, ncol = lc)

  # Enforce unit-1-norm constraint.
  current.row <- 1
  bvec[current.row]   <- 1
  Amat[current.row, ] <- rep(1, lc)
  current.row <- current.row + 1

  # Enforce limited redistribution constraint.
  # In the absence of a specified limit, this enforces positivity.
  bvec[current.row:(current.row + lc - 1)]   <- c / (sum(c) + k * u)
  Amat[current.row:(current.row + lc - 1), ] <- diag(lc)
  current.row <- current.row + lc

  # Enforce transformed hypercube constraints.
  for (i in 1:lc) {
    for (j in 1:lc) {
      if (i != j && c[j] != 0) {
        Amat[current.row, i] <- -c[j]
        Amat[current.row, j] <-  c[i] + u
        current.row <- current.row + 1
      }
    }
  }
  actual.row.count <- current.row - 1

  # Delete unused rows.
  length(bvec) <- actual.row.count
  Amat <- Amat[1:actual.row.count, ]

  # Finished building quadratic programming problem.
  # Solve quadratic programming problem.
  upper.index.bound.permissible.proportions <- solve.QP(
    Dmat = 2 * sim,
    dvec = rep(0, length(c)),
    Amat = t(Amat),
    bvec = bvec,
    meq  = 1,
    factorized = FALSE)$solution

  # Clamp solution to [0,1]
  if (min(upper.index.bound.permissible.proportions) < 0) {
    warning("Proportion below 0 clamped (", toString(min(upper.index.bound.permissible.proportions)), ")")
  }
  upper.index.bound.permissible.proportions[upper.index.bound.permissible.proportions < 0] <- 0
  if (max(upper.index.bound.permissible.proportions) > 1) {
    warning("Proportion above 1 clamped (", toString(max(upper.index.bound.permissible.proportions)), ")")
  }
  upper.index.bound.permissible.proportions[upper.index.bound.permissible.proportions > 1] <- 1

  # Add non-permitted disciplines
  upper.index.bound.proportions <- rep(0, n)
  upper.index.bound.proportions[which(permissible.disciplines)] <- upper.index.bound.permissible.proportions

  return(RaoProportions(upper.index.bound.proportions, similarity))
}
