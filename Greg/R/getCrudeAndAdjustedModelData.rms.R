#' @export
#' @importFrom rms bootcov
#' @importFrom rms robcov
#' @keywords internal
getCrudeAndAdjustedModelData.rms <- function(model, 
                                             level=.95, 
                                             remove_interaction_vars = TRUE, 
                                             remove_strata = FALSE,
                                             remove_cluster = FALSE,
                                             var_select,
                                             ...){
  
  var_names <- prGetModelVariables(model, 
                                   remove_splines = TRUE,
                                   add_intercept = FALSE,
                                   remove_interaction_vars = remove_interaction_vars)
  if (length(var_names) == 0)
    stop("You have no variables that can be displayed as adjusted/unadjusted",
         " since they all are part of an interaction, spline or strata.")

  if (grepl("intercept", coef(model)[1], ignore.case = TRUE))
    message("The rms-funcitons allow no easy way of extracting the intercept",
            ", it is thefore omitted from the output. Use regular lm/glm",
            " if the intercept is of interest")
  
  df <- prGetModelData(model)
  
  # Get the adjusted variables
  adjusted <- prCaRmsGetCoefAndCI(model, level = level,
                                  vn = var_names, data = df)
  
  if (!missing(var_select)){
    var_rows <- 
      prMapVariable2Name(var_names = var_names,
                         available_names = rownames(adjusted),
                         data = df)
    greps <- 
      prCaSelectAndOrderVars(names = names(var_rows),
                             order = var_select, 
                             ok2skip = TRUE)
    
    if (length(greps) == 0)
      stop("The variables that you have tried to select:",
           " '", paste(var_select, collapse="', '"), "'",
           " do not seem to exist - there is no match",
           " for the names: '", paste(names(var_rows), collapse="', '"), "'")
    
    var_rows <- var_rows[sort(greps)]
    var_names <- local({
      tmp <- var_names[var_names %in% names(var_rows)]
      copyAllNewAttributes(from = var_names, to = tmp)
    })
    
    adjusted <- adjusted[unlist(lapply(var_rows,
                                       function(x) x$location), 
                                use.names = FALSE), ,drop=FALSE]
  }
  
  unadjusted <- c()
  
  for(variable in var_names){
    
    # We should keep any strata information when running the models
    # TODO: Add the nlmn | options
    vars_4_frml <- variable
    if (!is.null(attr(var_names, "strata")) &&
          !remove_strata)
      vars_4_frml <- c(vars_4_frml, attr(var_names, "strata"))
    
    if (!is.null(attr(var_names, "cluster")) &&
          !remove_cluster)
      vars_4_frml <- c(vars_4_frml, attr(var_names, "cluster"))
    
    frml_4_single_var <- 
      paste(".~", paste(vars_4_frml,collapse="+"))
            
    # Run the same model but with only one variable
    model_only1 <- prEnvModelCall(model, update, frml_4_single_var)

    if ("boot.coef" %in% names(model)){
      # Redo bootstrap if this was a bootstrapped
      # rms model by rerunning the bootcov part
      model_only1 <- bootcov(model_only1, B=model$B, coef.reps="boot.Coef" %in% names(model))
    }else if (!is.null(attr(model, "robust"))){
      model_only1 <- robcov_alt(model_only1, type=attr(model, "robust"))
    }else if (!is.null(model$orig.var)){
      warning("Using the regular robcov function instead of the robcov_alt is not called for.",
        " If you get this message and you haven't used robcov then you may have a bug.", 
        " Try to set the model$org.var <- NULL and see if it works better.")
      model_only1 <- robcov(model_only1)
    }
    
    # Get the coefficients processed with some advanced
    # round part()
    new_vars <- prCaRmsGetCoefAndCI(model_only1, level = level,
                                    vn = variable, data = df)
    
    # Add them to the previous
    unadjusted <- rbind(unadjusted, new_vars)
  }
  
  # Fix if ordering got mixed up because the summary call
  if (any(rownames(adjusted) != rownames(unadjusted))){
    ua_order <- c()
    a_order <- c()
    ua_map <- prMapVariable2Name(var_names = var_names, 
                                 available_names = rownames(unadjusted),
                                 data = df)
    a_map <- prMapVariable2Name(var_names = var_names, 
                                available_names = rownames(adjusted),
                                data = df)
    for(variable in var_names){
      ua_order <- c(ua_order, 
                    ua_map[[variable]]$location)
      a_order <- c(a_order, 
                   a_map[[variable]]$location)
    }
    if (length(ua_order) != length(a_order)){
      stop ("An error happend when fetching the",
            "adjusted and unadjusted variables")
    }

    # Now reorder the matrices so they match
    unadjusted <- unadjusted[ua_order,]
    adjusted <- adjusted[a_order,]
  }
  
  # If just one variable it's not a proper matrix
  if (is.null(dim(adjusted))){
    both <- matrix(c(unadjusted, adjusted), nrow=1)
  }else{
    both <- cbind(unadjusted, adjusted)
  }
  
  levels_str <- c(sprintf("%.1f %%", 100*(1-level)/2),
    sprintf("%.1f %%", 100*(level + (1-level)/2)))
  
  colnames(both) <- c("Crude", levels_str,
                      "Adjusted", levels_str)
  
  attr(both, "model") <- model
  class(both) <- c("getCrudeAndAdjustedModelData", class(both))
  return(both)
}