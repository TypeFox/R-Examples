# This file contains all the helper funcitons that the outer exported
# functions utilize. I try to have a pr at the start of the name for all
# the private functions.
#
# Author: max
###############################################################################

#' Looks for unique rowname match without grep
#'
#' Since a rowname may contain characters reserved by regular
#' expressions I've found it easier to deal with the rowname
#' finding by just checking for matching strings at the beginning
#' of the name while at the same time excluding names that have the
#' same stem, i.e. DM and DM_COMP will cause an issue since DM will
#' match both rows.
#'
#' @param rnames A vector with the rownames that are looked for
#' @param vn The variable name that is of interest
#' @param vars A vector with all the names and the potentially competing names
#' @return integer A vector containing the position of the matches
#'
#' TODO: remove this function in favor of the more powerful prMapVariable2Name
#' @keywords internal
prFindRownameMatches <- function(rnames, vn, vars){
  # Find the beginning of the string that matches exactly to the var. name
  name_stub <- substr(rnames, 1, nchar(vn))
  matches <- which(name_stub == vn)
  
  # Since the beginning of the name may not be unique we need to
  # check for other "competing matches"
  # TODO: make this fix more elegant
  vars_name_stub <- substr(vars, 1, nchar(vn))
  if (sum(vars_name_stub == vn) > 1){
    competing_vars <- vars[vars != vn &
                             vars_name_stub == vn]
    
    competing_matches <- NULL
    for(comp_vn in competing_vars){
      competing_name_stub <- substr(rnames, 1, nchar(comp_vn))
      competing_matches <-
        c(competing_matches,
          which(competing_name_stub == comp_vn))
    }
    
    # Clean out competing matches
    matches <- matches[!matches %in% competing_matches]
  }
  
  return(matches)
}

#' Get model outcome
#'
#' Uses the model to extract the outcome variable. Throws
#' error if unable to find the outcome.
#'
#' @param model The fitted model
#' @param mf The dataset that the model is fitted to - if missing it
#'  uses the \code{\link[stats]{model.frame}} dataset. This can cause
#'  length issues as there may be variables that are excluded from the
#'  model for different reasons.
#' @return vector
#'
#' @keywords internal
prExtractOutcomeFromModel <- function(model, mf){
  if (missing(mf)){
    mf <- model.frame(model)
    outcome <- mf[,names(mf) == deparse(as.formula(model)[[2]])]
  }else{
    outcome <- eval(as.formula(model)[[2]], envir = mf)
  }
  if (is.null(outcome))
    stop("Could not identify the outcome: ", deparse(as.formula(model)[[2]]),
         " among the model.frame variables: '", paste(names(mf), collapse="', '"),"'")
  
  # Only use the status when used for survival::Surv objects
  if (inherits(outcome, "Surv"))
    return(outcome[,"status"])
  
  return(outcome)
}

#' Get model data.frame
#'
#' Returns the raw variables from the original data
#' frame using the \code{\link[stats]{get_all_vars}}
#' but with the twist that it also performs any associated
#' subsetting based on the model's \code{subset} argument.
#'
#' @param x The fitted model.
#' @return data.frame
#'
#' @keywords internal
prGetModelData <- function(x){
  # Extract the variable names
  true_vars <- all.vars(as.formula(x))
  
  # Get the environment of the formula
  env <- environment(as.formula(x))
  data <- eval(x$call$data,
               envir = env)
  
  # The data frame without the
  mf <- get_all_vars(as.formula(x),
                     data=data)
  
  if (!is.null(x$call$subset)){
    if (!is.null(data)){
      # As we don't know if the subsetting argument
      # contained data from the data frame or the environment
      # we need this additional check
      mf <- tryCatch(mf[eval(x$call$subset,
                             envir = data,
                             enclos = env), ],
                     error = function(e){
                       stop("Could not deduce the correct subset argument when extracting the data. ", e)
                     })
    }else{
      mf <- mf[eval(x$call$subset,
                    envir=env), ]
    }
  }
  
  return(mf)
}

#' Get the models variables
#'
#' This function extract the modelled variables. Any interaction
#' terms are removed as those should already be represented by
#' the individual terms.
#'
#' @param model A model fit
#' @param remove_splines If splines, etc. should be cleaned
#'  from the variables as these no longer are "pure" variables
#' @param remove_interaction_vars If interaction variables are
#'  not interesting then these should be removed. Often in
#'  the case of \code{\link{printCrudeAndAdjustedModel}} it is impossible
#'  to properly show interaction variables and it's better to show
#'  these in a separate table
#' @param add_intercept Adds the intercept if it exists
#' @return vector with names
#'
#' @importFrom stringr str_split
#' @importFrom stringr str_trim
#' @keywords internal
prGetModelVariables <- function(model,
                                remove_splines = TRUE,
                                remove_interaction_vars=FALSE,
                                add_intercept = FALSE){
  # We need the call names in order to identify
  # - interactions
  # - functions such as splines, I()
  if (inherits(model, "nlme")){
    vars <- attr(model$fixDF$terms, "names")
  }else{
    vars <- attr(model$terms, "term.labels")
  }
  
  strata <- NULL
  if (any(grepl("^strat[a]{0,1}\\(", vars))){
    strata <- vars[grep("^strat[a]{0,1}\\(", vars)]
    vars <- vars[-grep("^strat[a]{0,1}\\(", vars)]
  }
  
  cluster <- NULL
  if (any(grepl("^cluster{0,1}\\(", vars))){
    cluster <- vars[grep("^cluster{0,1}\\(", vars)]
    vars <- vars[-grep("^cluster{0,1}\\(", vars)]
  }
  # Fix for bug in cph
  if (is.null(cluster) &&
      inherits(model, "cph")){
    alt_terms <- stringr::str_trim(strsplit(deparse(model$call$formula[[3]]),
                                            "+", fixed = TRUE)[[1]])
    if (any(grepl("^cluster{0,1}\\(", alt_terms))){
      cluster <- alt_terms[grep("^cluster{0,1}\\(", alt_terms)]
    }
  }
  
  # Remove I() as these are not true variables
  unwanted_vars <- grep("^I\\(.*$", vars)
  if (length(unwanted_vars) > 0){
    attr(vars, "I() removed") <- vars[unwanted_vars]
    vars <- vars[-unwanted_vars]
  }
  
  pat <- "^[[:alpha:]\\.]+[^(]+\\(.*$"
  fn_vars <- grep(pat, vars)
  if(length(fn_vars) > 0){
    if (remove_splines){
      # Remove splines and other functions
      attr(vars, "functions removed") <- vars[fn_vars]
      vars <- vars[-fn_vars]
    }else{
      # Cleane the variable names into proper names
      # the assumption here is that the real variable
      # name is the first one in the parameters
      pat <- "^[[:alpha:]\\.]+.*\\(([^,)]+).*$"
      vars[fn_vars] <- sub(pat, "\\1", vars[fn_vars])
    }
  }
  
  # Remove interaction terms as these are not variables
  int_term <- "^.+:.+$"
  in_vars <- grep(int_term, vars)
  if (length(in_vars) > 0){
    if (remove_interaction_vars){
      in_vn <- unlist(str_split(vars[in_vars], ":"),
                      use.names = FALSE)
      in_vars <- unique(c(in_vars, which(vars %in% in_vn)))
    }
    attr(vars, "interactions removed") <- vars[in_vars]
    vars <- vars[-in_vars]
  }
  
  if (add_intercept &&
      grepl("intercept", names(coef(model))[1], ignore.case = TRUE)){
    vars <- c(names(coef(model))[1],
              vars)
  }
  
  clean_vars <- unique(vars)
  attributes(clean_vars) <- attributes(vars)
  if (!is.null(strata))
    attr(clean_vars, "strata") <- strata
  if (!is.null(cluster))
    attr(clean_vars, "cluster") <- cluster
  
  return(clean_vars)
}

#' Get statistics according to the type
#'
#' A simple function applied by the \code{\link[Gmisc]{getDescriptionStatsBy}}
#' for the total column. This function is also used by \code{\link{printCrudeAndAdjustedModel}}
#' in case of a basic linear regression is asked for a raw stat column
#'
#' @param x The variable that we want the statistics for
#' @param show_perc If this is a factor/proportion variable then we
#'  might want to show the percentages
#' @param html If the output should be in html or LaTeX formatting
#' @param digits Number of decimal digits
#' @param numbers_first If number is to be prior to the percentage
#' @param useNA If missing should be included
#' @param show_all_values This is by default false as for instance if there is
#'  no missing and there is only one variable then it is most sane to only show
#'  one option as the other one will just be a complement to the first. For instance
#'  sex - if you know gender then automatically you know the distribution of the
#'  other sex as it's 100 \% - other \%.
#' @param continuous_fn A function for describing continuous variables
#'  defaults to \code{\link{describeMean}}
#' @param prop_fn A function for describing proportions, defaults to
#'  the factor function
#' @param factor_fn A function for describing factors, defaults to
#'  \code{\link{describeFactors}}
#' @param percentage_sign If you want to suppress the percentage sign you
#'  can set this variable to FALSE. You can also choose something else that
#'  the default \% if you so wish by setting this variable.
#' @return A matrix or a vector depending on the settings
#'
#' TODO: Use the Gmisc function instead of this copy
#'
#' @importFrom Gmisc describeMean
#' @importFrom Gmisc describeFactors
#' @keywords internal
prGetStatistics <- function(x,
                            show_perc = FALSE,
                            html = TRUE,
                            digits = 1,
                            numbers_first = TRUE,
                            useNA = "no",
                            show_all_values = FALSE,
                            continuous_fn = describeMean,
                            factor_fn = describeFactors,
                            prop_fn = factor_fn,
                            percentage_sign = percentage_sign)
{
  useNA <- prConvertShowMissing(useNA)
  if (is.factor(x) ||
      is.logical(x) ||
      is.character(x)){
    if (length(unique(x)) == 2){
      if (show_perc){
        total_table <- prop_fn(x,
                               html=html,
                               digits=digits,
                               number_first=numbers_first,
                               useNA = useNA,
                               percentage_sign = percentage_sign)
      }else{
        total_table <- table(x, useNA=useNA)
        names(total_table)[is.na(names(total_table))] <- "Missing"
        # Choose only the reference level
        # Note: Currently references are required
        if (show_all_values == FALSE && FALSE)
          total_table <- total_table[names(total_table) %in% c(levels(x)[1], "Missing")]
      }
      
    } else {
      if (show_perc)
        total_table <- factor_fn(x,
                                 html=html,
                                 digits=digits,
                                 number_first=numbers_first,
                                 useNA = useNA,
                                 percentage_sign = percentage_sign)
      else{
        total_table <- table(x, useNA=useNA)
        names(total_table)[is.na(names(total_table))] <- "Missing"
      }
    }
  }else{
    total_table <- continuous_fn(x,
                                 html=html, digits=digits,
                                 number_first=numbers_first,
                                 useNA = useNA)
    
    # If a continuous variable has two rows then it's assumed that the second is the missing
    if (length(total_table) == 2 &&
        show_perc == FALSE)
      total_table[2] <- sum(is.na(x))
  }
  return(total_table)
}

#' Gets the boundaries for a survival fit
#'
#' @param fit A survival model of either competing risk regression or cox regression type
#' @param conf.int The interval of interest 0-1, see levels in confint()
#' @param exp If the value should be in exponential form (default)
#' @return A matrix with the columns:
#' \item{beta}{The estimated coefficient}
#' \item{p_val}{P-value}
#' \item{low}{The lower confidence interval}
#' \item{high}{The upper confidence interval}
#' \item{order}{A column that later can be used in ordering}
#'
#' @keywords internal
prGetFpDataFromSurvivalFit <- function (fit,
                                        conf.int = 0.95,
                                        exp      = TRUE){
  # Get the p-value, I use the method in the
  # print.cph -> prModFit from the rms package
  Z <- coef(fit)/sqrt(diag(fit$var))
  p_val <- signif(1 - pchisq(Z^2, 1), 5)
  order <- rep(-1, length(beta))
  ci <- confint(fit, level=conf.int)
  
  if (exp){
    ret_matrix <- cbind(
      beta=exp(coef(fit)),
      p_val=p_val,
      low=exp(ci[,1]),
      high=exp(ci[,2]),
      order=order)
  }else{
    ret_matrix <- cbind(
      beta=coef(fit),
      p_val=p_val,
      low=ci[,1],
      high=ci[,2],
      order=order)
  }
  
  # Set the names of the rows
  rownames(ret_matrix) <- names(fit$coef)
  
  return(ret_matrix)
}

#' Gets the boundaries for a GLM fit that is poisson or quasipoisson based
#'
#' @param glm.fit A regression model
#' @param conf.int The interval of interest 0-1, see levels in confint()
#' @param exp If the value should be in exponential form (default)
#' @return A matrix with the columns:
#' \item{beta}{The estimated coefficient}
#' \item{p_val}{P-value}
#' \item{low}{The lower confidence interval}
#' \item{high}{The upper confidence interval}
#' \item{order}{A column that later can be used in ordering}
#'
#' @keywords internal
prGetFpDataFromGlmFit <- function(glm.fit,
                                  conf.int = 0.95,
                                  exp      = TRUE){
  summary_glm <- summary.glm(glm.fit)
  
  # Extract the summary values of interest
  summary_se <- summary_glm$coefficients[,colnames(summary_glm$coefficients) == "Std. Error"]
  if ("quasipoisson" %in% glm.fit$family){
    summary_p_val <- summary_glm$coefficients[,colnames(summary_glm$coefficients) == "Pr(>|t|)"]
  }else if ("poisson" %in% glm.fit$family){
    summary_p_val <- summary_glm$coefficients[,colnames(summary_glm$coefficients) == "Pr(>|z|)"]
  }else{
    stop("Type of analysis not prepared!")
  }
  
  order = rep(-1, length(glm.fit$coefficients))
  ci <- confint(glm.fit, level=conf.int)
  
  if (exp){
    ret_matrix <- cbind(
      beta=exp(coef(glm.fit)),
      p_val=summary_p_val,
      low=exp(ci[,1]),
      high=exp(ci[,2]),
      order=order)
  }else{
    ret_matrix <- cbind(
      beta=coef(glm.fit),
      p_val=summary_p_val,
      low=ci[,1],
      high=ci[,2],
      order=order)
  }
  
  # Set the names of the rows
  rownames(ret_matrix) <- names(glm.fit$coefficients)
  
  # Remove the intercept
  ret_matrix <- ret_matrix[names(glm.fit$coefficients) != "(Intercept)", ]
  
  return(ret_matrix)
}


#' Gets the confidence interval, p-values,
#' coefficients from a survival object
#'
#' @param model_fit A regression fit from CRR, coxph, cph object
#' @param conf.int The interval of interest 0-1, see levels in confint()
#' @param exp If the value should be in exponential form (default)
#' @return A matrix with the columns:
#' \item{beta}{The estimated coefficient}
#' \item{p_val}{P-value}
#' \item{low}{The lower confidence interval}
#' \item{high}{The upper confidence interval}
#' \item{order}{A column that later can be used in ordering}
#'
#' @keywords internal
prGetFpDataFromFit <- function(model_fit,
                               conf.int = 0.95,
                               exp = TRUE){
  # Get the estimates, confidence intervals and the p_values
  if (any(class(model_fit) %in% "coxph") ||
      any(class(model_fit) %in% "crr")){
    sd <- prGetFpDataFromSurvivalFit(fit = model_fit, conf.int = conf.int, exp = exp)
  } else if (any(class(model_fit) %in% "glm")){
    sd <- prGetFpDataFromGlmFit(glm.fit = model_fit, conf.int = conf.int, exp = exp)
  } else {
    stop(paste("Unknown fit class type:", class(model_fit)))
  }
  
  return(sd)
}

#' A functuon for converting a useNA variable
#'
#' The variable is suppose to be directly compatible with
#' table(..., useNA=useNA). It throughs an error
#' if not compatible
#'
#' @param useNA Boolean or "no", "ifany", "always"
#' @return string
#'
#' @keywords internal
prConvertShowMissing <- function(useNA){
  if (useNA == FALSE || useNA == "no")
    useNA <- "no"
  else if (useNA == TRUE)
    useNA <- "ifany"
  
  if (!useNA %in% c("no", "ifany", "always"))
    stop(sprintf("You have set an invalid option for useNA variable, '%s' ,it should be boolean or one of the options: no, ifany or always.", useNA))
  
  return(useNA)
}

#' A function that tries to resolve what variable corresponds to what row
#'
#' As both the \code{\link{getCrudeAndAdjustedModelData}} and the
#' \code{\link{printCrudeAndAdjustedModel}} need to now exactly
#' what name from the \code{\link[stats]{coef}}/\code{\link[rms]{summary.rms}}
#' correspond to we for generalizeability this rather elaborate function.
#'
#' @param var_names The variable names that are saught after
#' @param available_names The names that are available to search through
#' @param data The data set that is saught after
#' @param force_match Whether all variables need to be identified or not.
#'  E.g. you may only want to use some variables and already pruned the
#'  \code{available_names} and therefore wont have matches. This is the
#'  case when \code{\link{getCrudeAndAdjustedModelData}} has been used together
#'  with the \code{var_select} argument.
#' @return \code{list} Returns a list with each element has the corresponding
#'  variable name and a subsequent list with the parameters \code{no_rows}
#'  and \code{location} indiciting the number of rows corresponding to that
#'  element and where those rows are located. For factors the list also contains
#'  \code{lvls} and \code{no_lvls}.
#' @keywords internal
#' @import utils
prMapVariable2Name <- function(var_names, available_names,
                               data, force_match = TRUE){
  if (any(duplicated(available_names)))
    stop("You have non-unique names. You probably need to adjust",
         " (1) variable names or (2) factor labels.")
  
  # Start with figuring out how many rows each variable
  var_data <- list()
  for (name in var_names){
    if (grepl("intercept", name, ignore.case = TRUE)){
      var_data[[name]] <-
        list(no_rows = 1)
    }else if (is.factor(data[,name])){
      var_data[[name]] <-
        list(lvls = levels(data[,name]))
      # Sometimes due to subsetting some factors don't exist
      # we therefore need to remove those not actually in the dataset
      var_data[[name]]$lvls <-
        var_data[[name]]$lvls[var_data[[name]]$lvls %in%
                                as.character(unique(data[, name][!is.na(data[, name])]))]
      var_data[[name]][["no_lvls"]] <- length(var_data[[name]]$lvls)
      var_data[[name]][["no_rows"]] <- length(var_data[[name]]$lvls) - 1
    }else{
      var_data[[name]] <-
        list(no_rows = 1)
    }
  }

  # A function for stripping the name and the additional information
  # from the available name in order to get the cleanest form
  getResidualCharacters <- function(search, conflicting_name){
    residual_chars <- substring(conflicting_name, nchar(search) + 1)
    if (!is.null(var_data[[search]]$lvls)){
      best_resid <- residual_chars

      for (lvl in var_data[[search]]$lvls){
        new_resid <- sub(lvl, "", residual_chars,
                         fixed = TRUE)
        if (nchar(new_resid) < nchar(best_resid)){
          best_resid <- new_resid
          if (nchar(new_resid) == 0)
            break;
        }
      }
      residual_chars <- best_resid
    }
    return(residual_chars)
  }

  matched_names <- c()
  matched_numbers <- c()
  org_available_names <- available_names
  # Start with simple non-factored variables as these should give a single-line match
  # then continue with the longest named variable
  for (name in var_names[order(sapply(var_data, function(x) is.null(x$lvls)),
                               nchar(var_names), decreasing = TRUE)]){
    matches <- which(name == substr(available_names, 1, nchar(name)))
    if (length(matches) == 0){
      if (force_match)
        stop("Sorry but the function could not find a match for '", name , "'",
             " among any of the available names: '", paste(org_available_names,
                                                           collapse="', '") ,"'")
    }else if(length(matches) == 1){
      if (var_data[[name]]$no_rows != 1)
        stop("Expected more than one match for varible '", name, "'",
             " the only positive match was '", available_names[matches], "'")

    }else if (length(var_names) > length(matched_names) + 1){
      if (is.null(var_data[[name]]$lvls) &&
            sum(name == available_names) == 1){
        # Check if the searched for variable is a non-factor variable
        # if so then match if there is a perfect match

        matches <- which(name == available_names)

      }else if (length(var_names) > length(matched_names) + 1){

        # Check that there is no conflicting match
        conflicting_vars <- var_names[var_names != name &
                                        !var_names %in% matched_names]
        possible_conflicts <- c()
        for (conf_var in conflicting_vars){
          possible_conflicts <-
            union(possible_conflicts,
                   which(substr(available_names, 1, nchar(conflicting_vars)) %in%
                           conflicting_vars))
        }
        conflicts <- intersect(possible_conflicts, matches)
        if (length(conflicts) > 0){

          conflicting_vars <- conflicting_vars[sapply(conflicting_vars,
                 function(search)
                   any(search == substr(available_names, 1, nchar(search))))]

          for (conflict in conflicts){
            # We will try to find a better match that leaves fewer "residual characters"
            # than what we started with
            start_res_chars <- getResidualCharacters(name, available_names[conflict])

            best_match <- NULL
            best_conf_name <- NULL
            for (conf_name in conflicting_vars){
              resid_chars <- getResidualCharacters(conf_name, available_names[conflict])
              if (is.null(best_match) ||
                    nchar(best_match) > nchar(resid_chars)){
                best_match <- resid_chars
                best_conf_name <- conf_name
              }
            }

            if (nchar(start_res_chars) == nchar(best_match)){
              stop("The software can't decide which name belongs to which variable.",
                   " The variable that is searched for is '", name, "'",
                   " and there is a conflict with the variable '", best_conf_name ,"'.",
                   " The best match for '", name, "' leaves: '", start_res_chars, "'",
                   " while the conflict '", best_conf_name ,"' leaves: '", best_match ,"'",
                   " when trying to match the name: '", available_names[conflict] ,"'")

            }else if(nchar(start_res_chars) > nchar(best_match)){
              # Now remove the matched row if we actually found a better match
              matches <- matches[matches != conflict]
            }
          }
        }
      }
      if (length(matches) == 0){
          stop("Could not identify the rows corresponding to the variable '", name ,"'",
               " this could possibly be to similarity between different variable names",
               " and factor levels. Try to make sure that all variable names are unique",
               " the variables that are currently looked for are:",
               " '", paste(var_names,
                           collapse="', '"),
               "'.")
      }
    }

    # Check that multiple matches are continuous, everything else is suspicious
    if (length(matches) > 1){
      matches <- matches[order(matches)]
      if (any(1 != tail(matches, length(matches) - 1) -
                head(matches, length(matches) -1)))
        stop("The variable '", name, "' failed to provide an adequate",
             " consequent number of matches, the names matched are located at:",
             " '", paste(matches, collapse="', '"), "'")
    }

    # Since we remove the matched names we need to look back at the original and
    # find the exact match in order to deduce the true number
    true_matches <- which(org_available_names %in%
                            available_names[matches])
    # Avoid accidentally rematching
    true_matches <- setdiff(true_matches, matched_numbers)
    var_data[[name]][["location"]] <- true_matches
    # Update the loop vars
    if (length(matches) > 0)
      available_names <- available_names[-matches]

    matched_names <- c(matched_names, name)
    matched_numbers <- c(matched_numbers, true_matches)

    if (length(var_data[[name]][["location"]]) == 0 &
          !force_match){
      # Remove variable as it is not available
      var_data[[name]] <- NULL
    }else if (length(var_data[[name]][["location"]]) !=
          var_data[[name]][["no_rows"]]){
      warning("Expected the variable '", name ,"'",
              " to contain '",var_data[[name]][["no_rows"]],"' no. rows",
              " but got '", length(var_data[[name]][["location"]]), "' no. rows.")
      var_data[[name]][["no_rows"]] <- length(var_data[[name]][["location"]])
    }
  }

  return(var_data)
}

#' Runs an \code{fastDoCall} within the environment of the model
#'
#' Sometimes the function can't find some of the variables that
#' were available when running the original variable. This function
#' uses the \code{\link[stats]{as.formula}} together with
#' \code{\link[base]{environment}} in order to get the environment
#' that the original code used.
#'
#' @param model The model used
#' @param what The function or non-empty character string used for
#'  \code{\link[Gmisc]{fastDoCall}}
#' @param ... Additional arguments passed to the function
#' @keywords internal
prEnvModelCall <- function(model, what, ...){
  call_lst <- list(object = model)
  dots <- list(...)
  if (length(dots) > 0){
    for(i in 1:length(dots)){
      if (!is.null(names(dots)[i])){
        call_lst[[names(dots)[i]]] <- dots[[i]]
      }else{
        call_lst <- c(call_lst,
                      dots[[i]])
      }
    }
  }
  model_env <- new.env(parent=environment(as.formula(model)))
  model_env$what <- what
  model_env$call_lst <- call_lst
  fastDoCall(what, call_lst,
             envir = model_env)
}