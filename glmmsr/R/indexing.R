
find_dim_sub <- function(x, var, d = NULL){
  if(is.atomic(x) || is.name(x)) {
    d
  } else if (is.call(x)) {
    if(identical(x[[1]], quote(`[`)) && identical(as.character(x[[2]]), var)) {
      if(length(d) > 0L){
        if(d != (length(x) - 2)) {
          stop(paste0("\'", var, "\' is indexed inconsistently"), call. = FALSE)
        }
      } else{
        d <- length(x) - 2
      }
    }
    for(i in seq_along(x)) {
      d <- find_dim_sub(x[[i]], var, d)
    }
  } else if(is.pairlist(x)) {
    for(i in seq_along(x)) {
      d <- find_dim_sub(x[[i]], var, d)
    }
  } else{
    stop("Don't know how to handle type ", typeof(x), call. = FALSE)
  }
  d
}

find_indices <- function(x, i, var) {
  if(is.atomic(x) || is.name(x)) {
    NULL
  } else if (is.call(x)) {
    if(identical(x[[1]], quote(`[`)) && identical(as.character(x[[2]]), var)) {
      return(as.character(x[[i+2]]))
    }
    unique(unlist(lapply(x, find_indices, i = i, var = var)))
  } else if(is.pairlist(x)) {
    unique(unlist(lapply(x, find_indices, i = i, var = var)))
  } else{
    stop("Don't know how to handle type ", typeof(x), call. = FALSE)
  }
}

# find all variables with ith index equal to index
find_vars_with_index <- function(x, i, index) {
  if(is.atomic(x) || is.name(x)) {
    NULL
  } else if (is.call(x)) {
    if(identical(x[[1]], quote(`[`)) && (length(x) > i + 1) &&
         identical(as.character(x[[i+2]]), index)) {
      return(as.character(x[[2]]))
    }
    unique(unlist(lapply(x, find_vars_with_index, i = i, index = index)))
  } else if(is.pairlist(x)) {
    unique(unlist(lapply(x, find_vars_with_index, i = i, index = index)))
  } else{
    stop("Don't know how to handle type ", typeof(x), call. = FALSE)
  }
}

find_indices_subform <- function(sub, data) {
  subvar <- sub$subvar
  subform <- sub$subform
  subexpr <- sub$subexpr

  d <- find_dim_sub(subexpr, subvar)

  # check that we get same dimension from subform
  if(find_dim_sub(subform, subvar) != d) {
    stop(paste0("\'", subvar, "\' is indexed inconsistently"), call. = FALSE)
  }

  # find the indices used to index subvar in subexpr
  indices_subexpr <- list()
  indices_subform <- list()
  for(i in 1:d){
    indices_subexpr[[i]] <- find_indices(subexpr, i, subvar)
    # check all subexpr indices are in data
    in_data <- vapply(indices_subexpr[[i]], exists, TRUE, where = data)
    if(any(!in_data)) {
      stop(paste0("Can't find indexing variables: ",
                  indices_subexpr[[i]][!in_data]))
    }
    indices_i <- mget(indices_subexpr[[i]], as.environment(data))
    index <- find_indices(subform, i, subvar)
    if(exists(index, data)) {
      # coerce to factor with common levels
      levels_i <- get(index, data)
      indices_i <- lapply(indices_i, factor, levels = levels_i)
      # replace in data, so have factors with correct levels
      for(j in seq_along(indices_i)) {
        data[[names(indices_i)[j]]] <- indices_i[[j]]
      }
    } else if(all(vapply(indices_i, is.factor, TRUE))) {
      levels_i_across_indices <- lapply(indices_i, levels)
      levels_i <- levels_i_across_indices[[1]]
      index_i <- factor(levels_i, levels = levels_i)
      # check that indexing the same across all indices
      if(any(!vapply(levels_i_across_indices, identical, TRUE, y = levels_i))) {
        stop("Indexing factors must have identical levels", call. = FALSE)
      }
    } else if(all(vapply(indices_i, is.numeric, TRUE))) {
      index_i <- sort(unique(Reduce(c, indices_i)))
    } else {
      stop("Should use factor or numeric to index variables", call. = FALSE)
    }
    #data <- drop_unused_levels(index, i, as.numeric(levels_i), subform, data)
    indices_subform[[i]] <- find_index_i(index, i, index_i, subform, data)
    names(indices_subform)[i] <- index
  }
  return(list(data = data, indices_subform = indices_subform))
}


find_index_i <- function(index, i, index_i, subform, data) {
  vars <- find_vars_with_index(subform, i, index)
  in_data <- vapply(vars, exists, TRUE, where = data)
  vars <- vars[in_data]
  x_vars <- data[vars]
  # for each var, we subset vars to keep only levels
  if(length(x_vars) > 0) {
    levels_i_max <- max(vapply(x_vars, find_dim_in_dir, 1L, i = i))
  } else {
    levels_i_max <- length(index_i)
  }
  if(is.factor(index_i)) {
    if(length(index_i) != levels_i_max) {
      stop(paste0("Indexing factor ", index, " has wrong number of levels"),
           call. = FALSE)
    } else{
      return(index_i)
    }
  } else {
    return(1:levels_i_max)
  }
}

matrix_indexing <- function(x, indices) {
  if(is.atomic(x) || is.name(x)) {
    x
  } else if (is.call(x)) {
    if(identical(x[[1]], quote(`[`)) && length(x) > 3) {
      y <- x[1:2]
      #indices[, as.list(x[-c(1,2)]), drop = FALSE]
      y[[3]] <- as.call(c(list(quote(cbind)), as.list(x[-c(1,2)])))
      return(y)
    } else{
      as.call(lapply(x, matrix_indexing, indices = indices))
    }
  } else if(is.pairlist(x)) {
    as.pairlist(lapply(x, matrix_indexing, indices = indices))

  } else{
    stop("Don't know how to handle type ", typeof(x), call. = FALSE)
  }
}

flatten_formula <- function(formula, indices){
  # drop indexing in LHS of formula
  lhs_var <- find_subvar(formula, char = FALSE)
  tf <- terms(formula)
  tf[[2]] <- lhs_var
  formula_flat <- formula(tf)

  d <- length(indices)

  if(d > 1) {
    # each time we see `[`, extract relevant columns of index_matrix
    # and combine with cbind
    # for now, insist all indexing is done with indices in indices

    formula_flat <- matrix_indexing(formula_flat, indices)
  }
  formula_flat
}


multi_to_flat <- function(multi_id, indices) {
  l <- vapply(indices, length, 1L)
  h <- c(1, cumprod(l)[-length(l)])
  unname(1 + colSums((multi_id-1)*h))
}

flatten_vars <- function(var_names, indices, data) {
  multi_id <- Reduce(rbind, lapply(data[var_names], as.numeric))
  multi_to_flat(multi_id, indices)
}

