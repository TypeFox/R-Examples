find_subexpr <- function(subvar) {
  parse(text = substr(subvar, 5, nchar(subvar) - 1))[[1]]
}



split_formula <- function(formula) {
  tf <- terms(formula, specials = "Sub")
  var_sub <- attr(tf, "specials")$Sub
  if(length(var_sub) > 0L) {
    # variables numbered starting at LHS.
    # assuming that there is exactly one LHS variable
    # could make this more general
    var_sub_rhs <- var_sub - 1L
    fac <- attr(tf, "factors")
    rhs_vars <- attr(fac, "dimnames")[[2]]
    if(length(rhs_vars) > length(var_sub_rhs)){
      attr(tf, "intercept") <- 1
      tf_no_sub <- drop.terms(tf, var_sub_rhs, keep.response = TRUE)
      attr(tf_no_sub, "intercept") <- 0
      form_no_sub <- fixFormulaObject(tf_no_sub)
    } else{
      form_no_sub <- update.formula(formula, . ~ 0)
    }
    subexprs <- lapply(rhs_vars[var_sub_rhs], find_subexpr)
  }
  else{
    form_no_sub <- formula
    subexprs <- NULL
  }
  return(list(form_no_sub = form_no_sub, subexprs = subexprs))
}

#drop all indexing from an expression
drop_index <- function(x) {
  if(is.atomic(x) || is.name(x)){
    x
  } else if (is.call(x)) {
    if(identical(x[[1]], quote(`[`))) {
      drop_index(x[[2]])
    } else {
      as.call(lapply(x, drop_index))
    }
  } else if (is.pairlist(x)) {
    as.pairlist(lapply(x, drop_index))
  } else{
    stop("Don't know how to handle type ", typeof(x), call. = FALSE)
  }
}

find_subvar <- function(subform, char = TRUE) {
  tf <- terms(subform)
  lhs <- tf[[2]]
  subvar <- drop_index(lhs)
  if(!is.name(subvar)) {
    stop("Should have only a single variable on LHS of subformula",
         call. = FALSE)
  }
  if(char) {
    res <- as.character(subvar)
  } else{
    res <- subvar
  }
  res
}

# find the vars (from vars), involved in x
find_vars <- function(x, vars) {
  if(is.atomic(x)) {
    NULL
  } else if (is.name(x)) {
    if(is.element(as.character(x), vars)){
      as.character(x)
    } else {
      NULL
    }
  } else if (is.call(x) || is.pairlist(x)) {
      unique(unlist(lapply(x, find_vars, vars = vars)))
  } else{
    stop("Don't know how to handle type ", typeof(x), call. = FALSE)
  }
}

add_subexpr <- function(subform, subexprs, which_subvars) {
  subvar <- find_subvar(subform)
  # find those subexpr involving (only) subvar
  which_subexpr <- which(vapply(which_subvars, function(x) {x[1] == subvar}, TRUE))
  if(length(which_subexpr) == 0L){
    warning(paste0("No subexpressions involving \'", subvar, "\'"), call. = FALSE)
    return(NULL)
  } else if(length(which_subexpr) > 1L) {
    stop(paste0("Multiple subexpressions involving ", subvar))
  }
  subexpr <- subexprs[[which_subexpr]]
  return(list(subvar = subvar, subform = subform, subexpr = subexpr))
}

match_subform_subexpr <- function(subforms, subexprs, data) {
  subvars <- vapply(subforms, find_subvar, "test")

  # find the names of the subvars involved in each subexpr
  which_subvars <- lapply(subexprs, find_vars, vars = subvars)
  # check how many subvars involved in each subexpr
  n_subvars <- vapply(which_subvars, length, 1L)
  if(any(n_subvars) > 1L) {
    stop("Each Sub(.) should only involve a single substituted variable",
      .call = FALSE)
  } else if (any(n_subvars) == 0L) {
    stop("Each Sub(.) should involve a substituted variable", .call = FALSE)
  }
  lapply(subforms, add_subexpr, subexprs = subexprs, which_subvars = which_subvars)
}


parse_sub <- function(sub, data, family)
{
  subvar <- sub$subvar
  subform <- sub$subform
  subexpr <- sub$subexpr

  indices_subform_tot <- find_indices_subform(sub, data)
  data <- indices_subform_tot$data
  indices_subform <- indices_subform_tot$indices_subform

  indices_expand <- expand.grid(indices_subform)
  subform_flat <- flatten_formula(subform, indices_expand)

  # find fake data for the subvar, of the correct length
  subvar_data <- list(rep(0, NROW(indices_expand)))
  names(subvar_data) <- subvar

  data_subform <- c(as.list(indices_expand), subvar_data, as.list(data))

  modfr_subform <- parse_subformula(subform_flat, data_subform)
  modfr_subform_list <- list(modfr_subform)
  names(modfr_subform_list) <- subvar

  # now need to substitute model frames into subexpr, in place of subvar
  # use a version of subexpr expecting model frames
  subexpr_mod <- modify_subexpr(subexpr, subvar)
  to_flatten <- extract_to_flatten(subexpr, subvar)
  subset_vars_flat <- lapply(to_flatten, flatten_vars,
                             indices = indices_subform,
                             data = data)

  names(subset_vars_flat) <- lapply(to_flatten, paste, collapse = "_")

  # we want to replace subvar with modfr_subform
  data_subexpr <- c(modfr_subform_list, subset_vars_flat, as.list(data))

  eval(subexpr_mod, envir = data_subexpr)
}

#' Parse a formula (and possibly subformulas)
#'
#' @inheritParams glmm
#' @export
find_modfr_glmm <- function (formula, subformula = NULL, data = NULL,
                             family = gaussian, weights = NULL, offset = NULL)
{
  if(is.list(subformula) || length(subformula) == 0L){
    subforms <- subformula
  } else {
    subforms <- list(subformula)
  }
  formula_split <- split_formula(formula)
  form_no_sub <- formula_split$form_no_sub
  subexprs <- formula_split$subexprs
  modfr_no_sub <- parse_formula(form_no_sub, data = data,
                                family = family, weights = weights,
                                off = offset)

  if(length(subexprs) == 0L) {
    return(modfr_no_sub)
  }
  subs <- match_subform_subexpr(subforms, subexprs, data)

  modfr_list <- lapply(subs, parse_sub, data = data, family = family)

  modfr <- attach_subframes(modfr_no_sub, modfr_list)
  modfr$formula <- formula
  modfr$subformula <- subformula

  modfr
}

