sb.check.args <- function(operation, nattributes, nlevels, nalternatives, 
  attribute.names, design, generators, effect, interactions, determinant,
  nblocks, seed) {
# Arguments: see sb.design()


# operation
  if (isTRUE(operation != "construct" & operation != "check")) {
    stop("'operation' is incorrect")
  }


# nattributes, nlevels, and nalternatives
  if (nattributes < 2 | nattributes > 20) {
    stop("'nattributes' must be integer values from 2 to 20")
  }
  if(any(nlevels < 2 | nlevels > 20)) {
    stop("'nlevels' must be integer values from 2 to 20")
  }
  if (!isTRUE(all.equal(nattributes, length(nlevels)))) {
    stop("'nattributes' must be equal to the length of 'nlevels'")
  }
  if (nalternatives < 2 | nalternatives > 20) {
    stop("'nalternatives' must be integer values from 2 to 20")
  }


# attribute.names
  if (!is.null(attribute.names)) {
    if (!is.list(attribute.names)) {
      stop("'attribute.names' must be list")
    }
    if (!isTRUE(all.equal(length(attribute.names), nattributes))) {
      stop("length of 'attribute.names' must be equal to 'nattributes'")
    }
    if (!isTRUE(all(sapply(attribute.names, length) == nlevels))) {
      stop("number of levels in each attribute in 'attribute.names' must be equal to 'nlevels' corresponding to the attribute")
    }
  }


# design
  if (operation == "check" & is.null(design)) {
    stop("'design' must be assigned when 'operation' is set as 'check'")
  }
  if (!is.null(design)) {
    if (!is.matrix(design)) {
      stop("'design' must be matrix")
    }
    storage.mode(design) <- "integer"
    if (!isTRUE(all.equal(min(design), 0))) {
        stop("minimum of all values in each column in 'design' must be 0")
    }
    if (operation == "check") {
      tmp.nlevels <- rep(nlevels, times = nalternatives)
    } else {
      tmp.nlevels <- nlevels
    }
    nlevels.list1 <- vector("list", ncol(design))
    for (i in 1:ncol(design)) {
      nlevels.list1[[i]] <- sort(unique(design[, i], )) 
    }
    nlevels.list2 <- lapply(tmp.nlevels, function(x) c(0:(x - 1)))
    if (!isTRUE(all.equal(nlevels.list1, nlevels.list2, check.attributes = FALSE))) {
      stop("values in each column in 'design' must correspond to those created from 'nlevels'")
    }
  }


# generators
  if (operation == "construct") {
    if (is.matrix(generators)) {
      ncol.generators <- ncol(generators)
    } else {
      ncol.generators <- length(generators)
    }
    if (!isTRUE(all.equal(ncol.generators, ((nalternatives - 1) * nattributes)))) {
      stop("length or number of columns in 'generators' must be equal to ('nalternatives' - 1) * 'nattributes'")
    }
    if (is.matrix(generators)) {
      if (!isTRUE(all(generators >= 0))) {
        stop("elements of 'generators' must be non-negative integers")
      }
      max.nlevels <- matrix(rep(nlevels, times = nrow(generators)), nrow = nrow(generators), byrow = TRUE) - 1
      if (!isTRUE(all(generators <= max.nlevels))) {
        stop("maximum of all values in each element of 'generators' must be ('nlevels' - 1) corresponding to the attribute")
      }
      generators.matrix <- matrix(as.vector(t(generators)), ncol = nattributes, byrow = TRUE)
      zero.matrix <- matrix(0, nrow = nrow(generators.matrix), ncol = ncol(generators.matrix))
      if (isTRUE(any(apply((generators.matrix == zero.matrix), 1, all)))) {
        stop("a set of generator with only elements of 0 is not permitted")
      }
    } else {
      generators.matrix <- matrix(generators, ncol = nattributes, byrow = TRUE)
      zero.matrix <- matrix(0, nrow = nrow(generators.matrix), ncol = ncol(generators.matrix))
      if (isTRUE(any(apply((generators.matrix == zero.matrix), 1, all)))) {
        stop("a set of generator with only elements of 0 is not permitted")
      }
      if (!isTRUE(all(generators >= 0))) {
        stop("each element of 'generators' must be non-negative integers")
      }
      if (!isTRUE(all(generators <= (nlevels - 1)))) {
        stop("maximum of all values in each element of 'generators' must be ('nlevels' - 1) corresponding to the attribute")
      }
    }
  }


# effect, interactions, and determinant
  if (isTRUE(effect != "main" & effect != "mplusall" & effect != "mplussome")) {
    stop("'effect' is incorrect")
  }
  if (effect == "main") {
    if (!isTRUE(is.null(determinant) & is.null(interactions))) {
      stop("'determinant' and 'interactions' must be NULL when 'effect' is set as 'main'")
    }
  } else {
    if (!is.null(determinant)) {
      if (!is.character(determinant)) {
        stop("'determinant' must be set as a character")
      }
      if (!isTRUE(type.convert(determinant) <= 1 & type.convert(determinant) >= 0)) {
        stop("'determinant' must be in the range [0, 1]")
      }
    }
  }
  if (effect == "mplusall") {
    if (!is.null(interactions)) {
      stop("'interactions' must be NULL when 'effect' is set as 'mplusall'")
    }
  }
  if (effect == "mplussome") {
    if (is.list(interactions)) {
      if (!isTRUE(all.equal(sapply(interactions, length), rep(2, length(interactions))))) {
        stop("each component of 'interactions' must be a vector with two elements")
      }
      if (isTRUE(any(sapply(interactions, function(x) identical(x[1], x[2]))))) {
        stop("values in each component of 'interactions' must differ")
      }
      if (!isTRUE((min(unlist(interactions)) >= 1) & (max(unlist(interactions)) <= nattributes))) {
        stop("values in each component of 'interactions' must be in the range [1, 'nattributes']")
      }
      if (isTRUE(any(duplicated(do.call("rbind", interactions), MARGINE = 1)))) {
        stop("components of 'interactions' must differ")
      }
    } else {
      if (!isTRUE(all.equal(length(interactions), 2))) {
        stop("length of 'interactions' must be 2")
      }
      if (isTRUE(all.equal(interactions[1], interactions[2]))) {
        stop("values in 'interactions' must differ")
      }
      if (!isTRUE((min(interactions) >= 1) & (max(interactions) <= nattributes))) {
        stop("values in 'interactions' must be in range [1, 'nattributes']")
      }
    }
  }


# nblocks
  if (operation == "check" & !isTRUE(all.equal(nblocks, 1))) {
    stop("'nblocks' must be defalut setting when 'operation' is set as 'check'")
  }
  if (operation == "construct") {
    if (!isTRUE(all.equal(nrow(design) %% nblocks, 0))) {
      stop("'nblocks' must be divisors of number of rows of 'design'")
    }
  }

}