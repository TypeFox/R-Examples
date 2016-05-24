#' Add a feature or parameter to a model
#' @export
#' @param e1 The Model to which the feature/parameter should be added
#' @param e2 The feature/parameter to add
#' @return The extended model
#' @keywords internal
"+.coalmodelpart" <- function(e1, e2) {
  e2name <- deparse(substitute(e2)) # Passed through for error messages

  if (is.model(e1)) return(add_to_model(e2, e1, e2name))

  if (is_partial_model(e1)) {
    e1[[length(e1) + 1]] <- e2
    return(e1)
  }

  partial_model <- list(e1, e2)
  class(partial_model) <- c("partial_model", "coalmodelpart")
  partial_model
}

is_partial_model <- function(x) inherits(x, "partial_model")

#' @export
print.partial_model <- function(x, ...) {
  cat("Partial coal_model with ", length(x), " components:\n")
  for (part in x) {
    cat("* ")
    print(part)
  }
}

add_to_model <- function(x, model, x_name) UseMethod("add_to_model")
add_to_model.default <- function(x, model, x_name) {
  stop("Can not add `", x_name, "` to model")
}


add_to_model.parameter <- function(par, model, par_name) model


add_to_model.named_par <- function(par, model, par_name) {
  if (par$get_name() %in% get_par_names(model))
    stop("There is already a parameter with name ", par_name)

  model$parameter[[length(model$parameter) + 1]] <- par

  model$id <- get_id()
  model
}


add_to_model.variation_par <- function(par, model, par_name) {
  model <- add_variation(model)
  for (par in par$get_base_par()) model <- model + par
  model$id <- get_id()
  model
}


add_to_model.feature <- function(feat, model, feat_name) {
  # Check that the population in the feature exists
  pop <- feat$get_population()
  if (!is.null(pop)) {
    if (pop != "all" && !all(pop %in% get_populations(model))) {
      stop("Invalid population in ", feat_name)
    }
  }

  # Execute the features checks
  feat$check(model)

  # Add the parameters in the feature
  for (para in feat$get_parameters()) model <- model + para

  # Add the feature itself
  model$features[[length(model$features) + 1]] <- feat

  model$id <- get_id()
  model
}


add_to_model.locus <- function(locus, model, locus_name) {
  model$loci[[length(model$loci) + 1]] <- locus
  model$id <- get_id()
  model
}


add_to_model.partial_model <- function(partial_model, model, partical_name) {
  for (part in partial_model) {
    model <- model + part
  }
  model
}
