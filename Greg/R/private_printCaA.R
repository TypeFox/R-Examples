#' Adds the ordering, references, and descriptions
#' 
#' This is a wrapper function around some more basic functions that
#' \code{\link{printCrudeAndAdjustedModel}} uses.
#'
#' @param x The main value matrix from the \code{\link{prCaPrepareCrudeAndAdjusted}}
#' @param model The model
#' @param order A vector A vector with regular expressions for each group.
#' @param var_order The output from the \code{\link{prMapVariable2Name}}
#' @param ds The dataset from the model
#'  
#' @return The reordered groups as a matrix
#' @family printCrudeAndAdjusted functions
#' @inheritParams printCrudeAndAdjustedModel
#' @keywords internal
prCaReorderReferenceDescribe <- function (
  x,
  model,
  order,
  var_order, 
  add_references, 
  add_references_pos, 
  reference_zero_effect, 
  ds, 
  desc_column, 
  desc_args, 
  use_labels)
{
  
  reordered_groups <- prCaReorder(mtrx2reorder = x,
                                  var_order = var_order,
                                  order = order)
  
  if (!missing(add_references)){
    var_order <- attr(reordered_groups, "var_order")

    if(add_references == TRUE){
      reordered_groups <- 
        prCaAddRefAndStat(model = model,
                          var_order = var_order,
                          add_references = add_references,
                          add_references_pos = add_references_pos,
                          reference_zero_effect = reference_zero_effect, 
                          values = reordered_groups, 
                          ds = ds,
                          desc_column = desc_column, 
                          desc_args = desc_args,
                          use_labels = use_labels)
    }else if (length(add_references) == 
                length(attr(reordered_groups, "greps"))){
      
      if (desc_column)
        warning("The descriptive column works so far only",
                " when used with automated references")
      
      reordered_groups <- 
        prCaAddUserReferences(reordered_groups = reordered_groups, 
                              var_order = var_order, 
                              add_references= add_references, 
                              add_references_pos = add_references_pos, 
                              reference_zero_effect = reference_zero_effect)
    }
  } 
  
  return(reordered_groups)
}

#' Function for retrieving the imputation arguments
#' 
#' @param impute_args The imputation arguments from \code{\link{printCrudeAndAdjustedModel}}
#'  function call.
#' @param output_mtrx The reordered groups matrix (a nx4 matrix)
#'  that have been prepared in for the \code{\link{printCrudeAndAdjustedModel}} 
#'  function. It is important that the references
#'  if any have been added.
#' @param model The imputation model. Currently only \code{\link[Hmisc]{fit.mult.impute}}
#'  is supported by the function.
#' @param data The data that has been used for generating the model.
#' 
#' @return \code{matrix} Returns a matrix with the requested columns
#' @family printCrudeAndAdjusted functions
#' @keywords internal
prCaGetImputationCols <- function(impute_args,
                                  output_mtrx,
                                  model,
                                  data){
  # Check if the reqquested imputation information has been implemented
  # if not return NULL
  if (!any(names(impute_args) %in% c("coef_change", 
                                     "variance.inflation"))){
    warning("You have specified imputation arguments currently not available.",
            " You can only specify coef_change or variance.inflation, see docs.")
    return(NULL)
  }
  
  impute_cols <- NULL
  custom_sprintf <- function(sp_str, variable){
    digits_search <- regexpr("%.[0-9]+f", sp_str)
    if(digits_search == -1)
      stop("Your output string for sprintf does not seem to contain a %.[0-9]+f",
           " regex compatible string, your output string is: ", sp_str)
    
    round_2 <- as.numeric(substr(sp_str, digits_search + 2, 
                                 digits_search + 2 + attr(digits_search, "match.length") - 4))
    n <- names(variable)
    # It's a little tricky to actually remove the -0.0 section...
    rounded_variable <- as.numeric(as.character(round(variable, digits = round_2)))
    out <- sprintf(sp_str, rounded_variable)
    names(out) <- n
    return(out)
  }
  
  # Compare the coefficients from the imputation with the 
  # original coefficients. This can give an idea to the 
  # direction of the imputed results, i.e. if the effect
  # size increases or decreases when adding the observations 
  # with missing data
  if (!is.null(impute_args$coef_change)){
    
    # Call the fitter without the imputed data in order
    # to get the original coefficients
    raw_call_lst <- list(formula=formula(model),
                         data=data)
    for (param_name in names(model$call)){
      # A few parameters are only used by the fit.mult.impute
      # and should not be forwarderd to the fitter
      if (!param_name %in% c("", "xtrans", "fitter", "n.impute",
                             "data", "formula", "fit.reps", 
                             "dtrans", "derived", "vcovOpts",
                             "pr")){
        raw_call_lst[[param_name]] <- model$call[[param_name]]
      }
    }
    non_imputed_fit <- 
      do.call(as.character(model$call$fitter),
              raw_call_lst)
    
    diff <- coef(model) - coef(non_imputed_fit)
    name <- "Coefficient change"
    
    
    # Do exp() if the variabels should be presented in that format
    # for that specific function type
    antilog <- FALSE
    if (inherits(model, "coxph") ||
          (!is.null(model$family) &&
             !is.null(model$family$link) &&
             grepl("^log", model$family$link))){
      antilog <- TRUE
    }
    
    if (!is.list(impute_args$coef_change)){
      change <- diff
      out_str <- "%.1f"
    }else{
      if (is.logical(impute_args$coef_change$antilog)){
        antilog <- impute_args$coef_change$antilog
      }
      
      out_str <- impute_args$coef_change$out_str
      
      if(is.null(impute_args$coef_change$type) ||
           tolower(impute_args$coef_change$type) %in% c("ratio",
                                                        "/")){
        if (is.null(out_str)){
          if (!is.numeric(impute_args$coef_change$digits)){
            out_str <- "%.2f"
          }else{
            out_str <- paste0("%.",
                              impute_args$coef_change$digits,
                              "f")
          }
        }
        change <- diff/abs(coef(non_imputed_fit))
        
        if (is.character(impute_args$coef_change$name)){
          name <- impute_args$coef_change$name
        }
        
      }else if (tolower(impute_args$coef_change$type) %in% c("%", 
                                                       "percent", 
                                                       "percentages")){
        if (is.null(out_str)){
          if (!is.numeric(impute_args$coef_change$digits)){
            out_str <- "%.0f%%"
          }else{
            out_str <- paste0("%.",
                              impute_args$coef_change$digits,
                              "f%%")
          }
        }
        change <- diff/abs(coef(non_imputed_fit))*100
        
        if (is.character(impute_args$coef_change$name)){
          name <- impute_args$coef_change$name
        }
      }else if(!is.character(impute_args$coef_change$type) ||
                 tolower(impute_args$coef_change$type) %in% c("abs",
                                                              "absolute",
                                                              "diff",
                                                              "difference",
                                                              "-")){
        if (is.null(out_str)){
          if (!is.numeric(impute_args$coef_change$digits)){
            out_str <- "%.2f"
          }else{
            out_str <- paste0("%.",
                              impute_args$coef_change$digits,
                              "f")
          }
        }
        change <- diff
        
        if (is.character(impute_args$coef_change$name)){
          name <- impute_args$coef_change$name
        }
        
      }else{
        stop("The requested type '", impute_args$coef_change$type ,"' for coef_reps",
             " is not yet implemented. Currently only percent or ratio is available.")
      }
    }
    
    if (antilog){
      change <- exp(change)
    }
    change <- custom_sprintf(out_str, change)
    impute_cols <- cbind(impute_cols,
                         change)
    
    colnames(impute_cols)[ncol(impute_cols)] <- name
  }
  
  if (!is.null(impute_args$variance.inflation)){
    
    name <- "Variance change"
    
    if (!is.list(impute_args$variance.inflation)){
      inflation <- custom_sprintf("%.2f", model$variance.inflation.impute)
    }else{
      if (is.character(impute_args$variance.inflation$name)) 
        name <- impute_args$variance.inflation$name
      
      if (is.character(impute_args$variance.inflation$type)){
        type <- tolower(impute_args$variance.inflation$type)
      }else{
        type <- "raw"
      }
      
      if (is.character(impute_args$variance.inflation$name)){
        name <- impute_args$variance.inflation$name
      }
      
      out_str <- impute_args$variance.inflation$out_str
      
      if (type %in% c("%", "percent", "percentages")){
        if (is.null(out_str)){
          if (!is.numeric(impute_args$variance.inflation$digits)){
            out_str <- "%.0f%%"
          }else{
            out_str <- paste0("%.",
                              impute_args$variance.inflation$digits,
                              "f%%")
          }
        }
        
        inflation <- custom_sprintf(out_str, model$variance.inflation.impute*100)
      }else if (type %in% c("raw", "ratio")){
        if (is.null(out_str)){
          if (!is.numeric(impute_args$variance.inflation.impute$digits)){
            out_str <- "%.2f"
          }else{
            out_str <- paste0("%.",
                              impute_args$variance.inflation.impute$digits,
                              "f")
          }
        }
        
        inflation <- custom_sprintf(out_str, model$variance.inflation.impute)
      }else{
        stop("The requested type '", impute_args$coef_change$type ,"'",
             " for variance.inflation.type",
             " is not yet implemented. Currently only percent or raw is available.")
      }
    }
    
    impute_cols <- cbind(impute_cols,
                         inflation)
    colnames(impute_cols)[ncol(impute_cols)] <- name
  }
  
  # To add additional imputation information we need to know
  # what rows contain the reference information
  reference_rows <- which(output_mtrx[,4] == "ref")
  
  if (length(reference_rows) > 0){
    for (row_no in reference_rows){
      impute_cols <-
        insertRowAndKeepAttr(impute_cols, r = row_no, 
                             v = rep("", times=ncol(impute_cols)))
      
    }
  }
  
  # The rms version does not provide the intercept and therefore the
  # intercept should be removed
  if (inherits(model, "rms") &&
        !grepl("intercept", x = rownames(output_mtrx)[1],  ignore.case = TRUE) &&
        grepl("intercept", x = rownames(impute_cols)[1],  ignore.case = TRUE)){
    impute_cols <- impute_cols[-1,,drop=FALSE]
  }
  
  return(impute_cols)
}

#' Add reference according to the model
#' 
#' This is of course for factored variables and not in general.
#'  
#' @param var_order The output from the \code{\link{prMapVariable2Name}}
#' @param values The values that are to be outputted
#' @param ds The dataset
#' @return \code{list}
#' 
#' @inheritParams printCrudeAndAdjustedModel
#' @importFrom Gmisc copyAllNewAttributes
#' @family printCrudeAndAdjusted functions
#' @keywords internal
prCaAddRefAndStat <- function(model, 
                              var_order, 
                              add_references, # Only called when present
                              add_references_pos,
                              reference_zero_effect, 
                              values, 
                              ds,
                              desc_column, 
                              desc_args,
                              use_labels){
  
  outcome <- NULL
  if (desc_column){
    stats <- list()
    
    # Get the original data
    outcome <- prExtractOutcomeFromModel(model, ds)
    if (is.matrix(outcome) && 
          "coxph" %in% class(model)){
      # Get the left part of the formula
      outcome <- outcome[,"status"]
    }
  }
  
  stats <- list()
  
  for(vn in names(var_order))
  {
    if (desc_column && !is.null(outcome)){
      stats[[vn]] <- prCaGetVnStats(model = model,
                                    vn = vn, 
                                    outcome = outcome,
                                    ds = ds,
                                    add_references = add_references,
                                    add_references_pos = add_references_pos,
                                    desc_args = desc_args)
    }
    
    # Add the refrence to the values matrix if it is a factor variable with a ref.
    if (!is.null(var_order[[vn]]$lvls)){
      values <- prCaAddReference(vn = vn, 
                                 var_order = var_order, 
                                 values = values,
                                 add_references_pos = add_references_pos,
                                 reference_zero_effect = reference_zero_effect,
                                 ds = ds,
                                 use_labels = use_labels)
      var_order <- attr(values, "var_order")
      
      # Update locations after the changed variable rows
      current_var_no <- which(vn == names(var_order))
      if (current_var_no < length(var_order)){
        for(i in (current_var_no + 1):length(var_order)){
          var_order[[i]]$location <- var_order[[i]]$location + 1
        }
      }
    }
  }
  
  
  if (desc_column){
    desc_mtrx <- matrix("-", 
                        ncol=NCOL(stats[[1]]), 
                        nrow=NROW(values))
    rownames(desc_mtrx) <- rownames(values)
    
    # Should probably make sure we're always dealing
    # with a matrix but this is a quick fix for now
    # TODO: fix consistent matrix handling
    getRows <- function(x){
      ifelse(is.matrix(x), nrow(x), length(x))
    }
    getRownames <- function(x){
      if(is.matrix(x))
        rownames(x)
      else
        names(x)
    }
    getValue <- function(x, rn){
      if(is.matrix(x))
        x[rn,,drop=FALSE] 
      else
        x[rn]
    }
    for(vn in names(var_order)){
      if (!is.null(var_order[[vn]]$lvls)){
        existing_labels <- rownames(values)[var_order[[vn]]$location]
        if (any(!existing_labels %in% rownames(stats[[vn]])))
          stop(paste("A few labels from factor", vn, "weren't found in the stats:",
                     paste(existing_labels[!existing_labels %in% rownames(stats[[vn]])], 
                           collapse=", ")))
        
        # Add the stats to the desc
        # This slightly complicated structure is to make sure that
        # the descriptive value corresponds to the regression row
        for(rn in existing_labels){
          # Find that row within the group
          group_rownames <- rownames(desc_mtrx[var_order[[vn]]$location, ,drop=FALSE])
          row_within_group <- which(rn == group_rownames)
          if (length(row_within_group) > 1)
            stop("There are more than one occurrence within group ", vn,
                 "that have the value: '", rn, "'\n",
                 "The rownames in that group are: '", paste(group_rownames, "', '"), "'")
          else if (length(row_within_group) == 0)
            stop("There was no match within group ", vn,
                 "for the value: '", rn, "'\n",
                 "The rownames in that group are: '", paste(group_rownames, "', '"), "'")
          
          # Set the value of that row
          desc_mtrx[head(var_order[[vn]]$location, 1) + 
                      row_within_group - 1, ] <- getValue(stats[[vn]], rn)
        }
        
        # There are more values in the stats than in the 
        # regression, this is probably due to missing values,
        # these will be added to the current group last
        if (getRows(stats[[vn]]) > var_order[[vn]]$no_rows){
          rows_2_add <- getRownames(stats[[vn]])[!getRownames(stats[[vn]]) %in% 
                                                   existing_labels]
          for (i in 1:length(rows_2_add)){
            rn <- rows_2_add[i]
            values <- insertRowAndKeepAttr(values, 
                                           r = tail(var_order[[vn]]$location, 1) + 1, 
                                           v = rep("-", length.out=NCOL(desc_mtrx)), 
                                           rName = rn)
            desc_mtrx <- insertRowAndKeepAttr(desc_mtrx, 
                                              r = tail(var_order[[vn]]$location, 1) + 1,
                                              v = getValue(stats[[vn]], rn), 
                                              rName = rn)
          }
          
          var_no <- which(vn == names(var_order))
          var_order[[var_no]]$no_rows <- var_order[[var_no]]$no_rows + length(rows_2_add)
          
          # Update the value_order
          for(ii in var_no:length(var_order)){
            start_pos <- 0
            if (ii > 1){
              start_pos <- tail(var_order[[ii - 1]]$location, 1) - 1
            }
            var_order[[ii]]$location <- start_pos + 1:var_order[[ii]]$no_rows
          }
        }
      }else{
        # This occurrs if the element is logical and you have
        # a TRUE/FALSE situation
        
        if (var_order[[vn]]$no_rows != NROW(stats[[vn]])){
          if (var_order[[vn]]$no_rows == 1 &&
                NROW(stats[[vn]]) == 2 &&
                is.logical(ds[,vn])){
            stats[[vn]] <- stats[[vn]]["TRUE",]
          }else{
            stop("The description statistics did not work for '", vn, "'",
                 " it returned ", NROW(stats[[vn]]), " rows",
                 " while expecting ", var_order[[vn]]$no_rows, " row(s).",
                 " The rowlabels returned are:",
                 " '", paste(rownames(stats[[vn]]), collapse="', '"), "'")
          }
        }
        desc_mtrx[var_order[[vn]]$location, ] <- stats[[vn]]
      }
    }
    
    values <- copyAllNewAttributes(from = values, 
                                   to = cbind(desc_mtrx, values)) 
    if (NCOL(desc_mtrx) == 1 && colnames(values)[1] == "")
      colnames(values)[1] <- desc_args$colnames[1]
    else if(all(colnames(values)[1:2] == ""))
      colnames(values)[1:2] <- desc_args$colnames
  }
  
  attr(values, "var_order") <- var_order
  return(values)
}

#' Adds a reference to value matrix
#' 
#' @param vn Variable name 
#' @param var_order The output from the \code{\link{prMapVariable2Name}}
#' @param values The value matrix
#' @param ds The data set
#' @return \code{matrix} A matrix with rgroup and n.rgroup attributes 
#' 
#' @inheritParams printCrudeAndAdjustedModel
#' @family printCrudeAndAdjusted functions
#' @keywords internal
prCaAddReference <- function(vn, 
                             var_order, 
                             values, 
                             add_references_pos, 
                             reference_zero_effect, 
                             ds, 
                             use_labels){
  ref_value <- rep(c(reference_zero_effect, "ref"), times=2)
  
  reference <- NULL
  rms_format <- FALSE
  # The rms package generates rownames with factor name:reference factor
  # and it is therefore a good idea to find the refreence by checking
  # which one is at the end
  for (f_name in var_order[[vn]]$lvls){
    # The substr is just to avoid having to check for regular expression
    # characters within f_name and escaping them
    beginning <- substr(rownames(values)[var_order[[vn]]$location], 
                        1,
                        nchar(rownames(values)[var_order[[vn]]$location])-
                          nchar(f_name))
    end <- substring(rownames(values)[var_order[[vn]]$location], 
                     nchar(rownames(values)[var_order[[vn]]$location])-
                       nchar(f_name) + 1)
    if (all(grepl(":$", beginning)) && 
          all(f_name %in% end)){
      reference <- f_name
      rms_format <- TRUE
      break
    }
  }
  
  if (is.null(reference)){
    # TODO: Could probably be extended to other regression models but needs testing
    used_factors <- gsub("^[ ]{0,1}[=-]{0,1}{0,1}", "", 
                         substring(rownames(values)[var_order[[vn]]$location], nchar(vn) + 1))
    
    # Fetch the reference level, would probably work just as well with a levels()[1]
    reference <- var_order[[vn]]$lvls[!var_order[[vn]]$lvls %in% used_factors]
    if (length(reference) != 1)
      stop("Error occurred in looking for reference,",
           " found ", length(reference), " reference categories",
           " while expecting 1, out of these factors:",
           "\n '", paste(var_order[[vn]]$lvls, collapse="'', '"), "'",
           ifelse(length(reference) > 1, 
                  sprintf(" \n The refrences found: '%s'",
                          paste(reference, collapse="', '")),
                  sprintf(" \n The rownames that have searched: '%s'", 
                          paste(rownames(values)[var_order[[vn]]$location], collapse="', '"))))
  }else{
    used_factors <- var_order[[vn]]$lvls[reference != var_order[[vn]]$lvls]
  }
  
  clean_rn = rownames(values)
  for (uf_name in used_factors){
    if (rms_format){
      clean_rn <- gsub("^.* - (.*):.*$", "\\1", rownames(values))
      r_no <- which(clean_rn == uf_name)
    }else{
      r_no <- grep(uf_name, clean_rn, fixed=TRUE)
    }
    
    if (!any(r_no %in% var_order[[vn]]$location))
      stop("Could not find rowname with factor ", uf_name, 
           " among any of the row names: ", paste(clean_rn, collapse=", "))
    
    r_no <- r_no[r_no %in% var_order[[vn]]$location]
    if (length(r_no) > 1)
      stop("Multiple rows matched the factor ", uf_name,
           " from the available: ", paste(clean_rn[var_order[[vn]]$location], collapse=", "))
    
    # Remove the main label as that goes into the attr(values, "rgroup")
    rownames(values)[r_no] <- uf_name
  }
  
  offset <- ifelse(vn %in% names(add_references_pos),
                   add_references_pos[[vn]] - 1,
                   0)
  if (offset > length(var_order[[vn]]$lvls) - 1 ||
        offset < 0){
    warning("You have a reference position '", add_references_pos[[vn]], "'",
            " that is outside the number of levels, '", offset + 1, "'",
            " is not among, 1 to ", length(var_order[[vn]]$lvls),
            ". This will therefore be ignored")
    offset <- 0
  }
  
  values <- insertRowAndKeepAttr(values, 
                                 var_order[[vn]]$location[1] + offset, 
                                 ref_value,
                                 rName=reference)

  var_order[[vn]] <- 
    within(var_order[[vn]],
           {location <- c(location[1] + offset, 
                         location[location <= offset], 
                         location[location > offset] + 1)
            no_rows <- length(location)
            reference_pos <- offset
            })
  attr(values, "var_order") <- var_order
  
  return(values)
}

#' Gets the variable stats
#' 
#' @param model The model
#' @param vn The variable name
#' @param outcome The outcome vector
#' @param ds The dataset
#'  
#' @return \code{matrix} A matrix from \code{\link[Gmisc]{getDescriptionStatsBy}} or
#'  \code{\link{prGetStatistics}}
#' @importFrom Gmisc getDescriptionStatsBy
#' @family printCrudeAndAdjusted functions
#' @inheritParams printCrudeAndAdjustedModel
#' @keywords internal
prCaGetVnStats <- function(model,
                           vn, 
                           outcome, 
                           ds,
                           add_references, # Only called when present
                           add_references_pos,
                           desc_args){
  # TODO: add some option of handling missing from the model, a second/third column
  # TODO: add handling for logical values
  
  # If there is a binomial outcome variable then 
  # it makes sense to have two columns, the overall
  # and the event data.
  if (any(class(model) %in% c("lrm", "coxph")) ||
        ("glm" %in% class(model) &&
           model$family$family == "binomial")){
    if (grepl("intercept", vn, ignore.case = TRUE)){
      desc_mtrx <- matrix("-", 
                          ncol = length(unique(outcome[is.na(outcome) == FALSE])),
                          nrow = 1)
    }else{
      desc_mtrx <- 
        getDescriptionStatsBy(x=ds[is.na(outcome) == FALSE,vn], 
                              by=outcome[is.na(outcome) == FALSE],
                              hrzl_prop = TRUE,
                              digits = desc_args$digits,
                              continuous_fn = desc_args$continuous_fn,
                              prop_fn = desc_args$prop_fn,
                              factor_fn = desc_args$factor_fn,
                              show_all_values = !missing(add_references),
                              useNA = desc_args$useNA,
                              add_total_col = TRUE,
                              total_col_show_perc = desc_args$show_tot_perc, 
                              html = TRUE)
    }
    
    # Don't select the no-event alternative as this is usually
    # not interesting since we have the total column
    desc_mtrx <- desc_mtrx[,c(1,3),drop=FALSE]
    colnames(desc_mtrx) <- desc_args$colnames
  }else{
    if (grepl("intercept", vn, ignore.case = TRUE)){
      desc_mtrx <- matrix("-", 
                          ncol = 1,
                          nrow = 1)
    }else{
      desc_mtrx <- 
        prGetStatistics(x=ds[is.na(outcome) == FALSE,vn],  
                        show_perc = desc_args$show_tot_perc, 
                        html = TRUE,
                        digits = desc_args$digits,
                        continuous_fn = desc_args$continuous_fn,
                        prop_fn = desc_args$prop_fn,
                        factor_fn = desc_args$factor_fn,
                        useNA = desc_args$useNA)
    }
  }
  
  if (!is.matrix(desc_mtrx)){
    rn <- names(desc_mtrx)
    desc_mtrx <- matrix(desc_mtrx, ncol=1)
    if (!is.null(rn))
      rownames(desc_mtrx) <- rn
  }
  
  # As the first element in a factor is always the
  # reference then we need to move it to the wanted
  # position
  if (!grepl("intercept", vn, ignore.case = TRUE) &&
        is.factor(ds[,vn]) && 
        vn %in% names(add_references_pos) &&
        add_references_pos[[vn]] != 1){
    
    if (nrow(desc_mtrx) == 2){
      if (add_references_pos[[vn]] == 2)
        desc_mtrx <- desc_mtrx[c(2,1), ]
    }else if (nrow(desc_mtrx) > 2){
      if (add_references_pos[[vn]] == nrow(desc_mtrx)){
        desc_mtrx <- desc_mtrx[c(2:nrow(desc_mtrx), 1), ] 
      }else{
        desc_mtrx <- desc_mtrx[c(2:add_references_pos[[vn]], 
                                 1,
                                 (add_references_pos[[vn]]+1):nrow(desc_mtrx)), ]
      }
    }
  }
  
  return(desc_mtrx)
}

#' Gets the labelled rowname if it exists
#' 
#' Looks for matches inside factors if rowname
#' contains the name of the column. 
#' 
#' @param vn The variable name 
#' @param use_labels If labels should be used
#' @param dataset The dataset
#' @return \code{string} The rowname 
#' 
#' @importFrom Hmisc label
#' @family printCrudeAndAdjusted functions
#' @keywords internal
prCaGetRowname <- function(vn, use_labels, dataset){
  vn <- as.character(vn)
  if(vn %in% colnames(dataset) &&
       use_labels && 
       label(dataset[,vn]) != ""){
    return(label(dataset[,vn]))
  }else if (any(vn == colnames(dataset))){
    # An exact match means that there is no factor information
    # for this row and we should be able to return this row
    return(vn)
  }else if (grepl("intercept", vn, ignore.case = TRUE)){
    return(vn)
  }
  
  # Check if this is actually a factor and return that factors name
  colno_containing_name <-
    unlist(lapply(colnames(dataset), 
                  function(x) grepl(x, vn, fixed=TRUE)))
  if (sum(colno_containing_name) == 1){
    # Remove the column name from the beginning of the vn 
    # as this may otherwise cause a search conflict if
    # the name consists the searched labels
    cn <- colnames(dataset)[colno_containing_name]
    if (cn == substr(vn, 1, nchar(cn))){
      vn <- substring(vn, nchar(cn) + 1)
    }
    
    lvls <- levels(dataset[,colno_containing_name])
    matching_lvl <- unlist(lapply(lvls, function(x) grepl(x, vn, fixed=TRUE)), 
                           use.names = FALSE)
    if (sum(matching_lvl) == 1)
      return(lvls[matching_lvl])
    
    # The rms-package returns the levels so that the reference appears after
    # each level and thus may give the appearance that the levels appears
    # several times. We therefore need to remove that information.
    if (grepl(sprintf("%s$", lvls[1]), vn)){
      vn <- gsub(sprintf("%s$", lvls[1]), "", vn)
      matching_lvl <- unlist(lapply(lvls, function(x) grepl(x, vn, fixed=TRUE)), 
                             use.names = FALSE)
      if (sum(matching_lvl) == 1)
        return(lvls[matching_lvl])
    }
    
    
    warning("Could not identify the rowname '", vn, "'",
            " from the factor variable '", colnames(dataset)[colno_containing_name], "'",
            " that has the factors: '", paste(lvls, collapse="', '"), "'",
            " The rowname with label therefore defaults to ", vn)
  }
  
  return(vn)
}

#' Sets the rownames of the reordered_groups
#' 
#' @param reordered_groups The value matrix that needs refrences 
#' @param var_order The output from the \code{\link{prMapVariable2Name}}
#' @param rownames.fn A rowname function for tailoring names
#' @param use_labels Whether to use labels or not
#' @param ds The model data set
#'  
#' @return \code{matrix} Returns the reordered_groups
#' @family printCrudeAndAdjusted functions
#' @keywords internal
prCaSetRownames <- function (reordered_groups,
                             var_order,
                             rowname.fn, 
                             use_labels, 
                             ds) {
  if (!missing(rowname.fn)){
    if (is.character(rowname.fn))
      rowname.fn <- get(rowname.fn)
    
    rn <- list()
    for (name in rownames(reordered_groups)){
      new_name <- rowname.fn(name)
      if (new_name == name)
        new_name <- prCaGetRowname(vn = name, use_labels = use_labels, dataset = ds)
      rn <- append(rn, new_name)
    }
  }else{
    rn <- rownames(reordered_groups)
    for (name in names(var_order)){
      # Only change the names of variables if they are not factors
      # Factors already have variable names assigned to them
      if (is.null(var_order[[name]]$lvls) &&
            var_order[[name]]$no_rows == 1){
        rn[var_order[[name]]$location] <- 
          prCaGetRowname(vn = name, 
                         use_labels = use_labels, 
                         dataset = ds)
      }else if (!is.null(var_order[[name]]$lvls)){
        # This is one easy way to remove the reference
        new_names <- tail(var_order[[name]]$lvls, 
                          var_order[[name]]$no_rows) 
        # Sanity check - we don't want to rename rows 
        # without knowing that there is some substance to it
        name_check_ok <- TRUE
        for (i in 1:length(new_names)){
          cn <- rn[var_order[[name]]$location[i]]
          # Remove the variable name
          cn <- gsub(name, "", cn, fixed = TRUE)
          if (!grepl(new_names[i], cn, fixed = TRUE)){
            name_check_ok <- FALSE
            warning("Tried to set the rownames for the variable's '", name ,"'",
                    " to the name of its level '", new_names[i] ,"'",
                    " but could not match the level to the row '", name ,"'",
                    " after subtracting the variable name, i.e. leaving '", cn ,"'")
            break;
          }
        }
        
        if (name_check_ok)
          rn[var_order[[name]]$location] <- new_names
      }
    }
  }
  rownames(reordered_groups) <- unlist(rn, use.names = FALSE)
  return(reordered_groups)
}

#' Re-order variables
#' 
#' @param names The names of the variables 
#' @param order The order regular expression
#' @param ok2skip If you have the intercept then
#'  it should be ok for the function to skip that
#'  variable if it isn't found among the variable list
#'  
#' @return \code{vector} A vector containing the greps 
#' @keywords intrenal
#' @family printCrudeAndAdjusted functions
#' @keywords internal
prCaSelectAndOrderVars <- function(names, order, ok2skip = FALSE){
  greps <- c()
  for (r_expr in order) {
    # Find the names that matches
    matches <- grep(r_expr, names)
    if (length(matches) == 0 & !ok2skip){
      stop("You have a strange selection order,",
           "this could be due to that you try to select a factor level",
           " and not the full variable.",
           " Re-arranging factors should be done in the factor() function and not here.",
           " Anyway the expression '", r_expr, "' was not found",
           " in these variable names: '",
           paste(names, collapse="', '"), "'")
    }else if (length(matches) > 0){
      # Avoid reselecting
      new_vars <- setdiff(matches, unlist(greps, 
                                          use.names = FALSE))
      if (length(new_vars) > 0){
        greps <- append(greps, list(new_vars))
      }
    }
  }
  return(unlist(greps, 
                use.names = FALSE))
}


#' Reorder according to the requested variables
#' 
#' Uses the \code{\link{prCaSelectAndOrderVars}} for finding the 
#' orders according to the \code{order} argument.
#' 
#' @param mtrx2reorder The matrix to reorder
#' @param var_order The variables representing different rows 
#'  \code{\link{prMapVariable2Name}}
#' @param order A vector of strings used for \code{\link{prCaSelectAndOrderVars}}
#' 
#' @return \code{matrix} Returns the \code{mtrx2reorder} rearranged with the 
#'  attribute "greps" for the greps from \code{\link{prCaSelectAndOrderVars}}
#'  and the attribute "var_order" for the new var_order
#' @family printCrudeAndAdjusted functions
#' @keywords internal
prCaReorder <- function (mtrx2reorder, var_order, order) {
  greps <- prCaSelectAndOrderVars(names = names(var_order), 
                                  order = order,
                                  ok2skip = TRUE)
  var_order <- var_order[greps]
  row_reorder <- c()
  for (i in 1:length(var_order)){
    last_pos <- length(row_reorder)
    row_reorder <- c(row_reorder,
                     var_order[[i]]$location)
    var_order[[i]]$location <- last_pos + 1:var_order[[i]]$no_rows
  }
  
  reordered_groups <- mtrx2reorder[row_reorder, ,drop=FALSE]
  attr(reordered_groups, "greps") <- greps
  attr(reordered_groups, "var_order") <- var_order
  
  # Just some warnings/checks
  if (any(!rownames(reordered_groups) %in% 
            rownames(mtrx2reorder))){
    groups_not_fount <-
      rownames(reordered_groups)[!rownames(reordered_groups) %in% 
                                   rownames(mtrx2reorder)]
    stop("An error occurred when reordering, there are now more",
         " variables than initially found, the following new",
         " vars exist: '", paste(groups_not_fount, collapse="', '"), "'")
    
  }else if (any(!rownames(mtrx2reorder) %in% 
                  rownames(reordered_groups))){
    rows_not_used <-
      rownames(mtrx2reorder)[!rownames(mtrx2reorder) %in% 
                               rownames(reordered_groups)]
    warning("Not all variables selected from the model when re-ordering,",
            " the following were not included:",
            " '", paste(rows_not_used, collapse="', '"), "'")
  }
  
  return(reordered_groups)
}

#' Adds references
#'
#' @param reordered_groups The value matrix that needs refrences 
#' @param var_order The output from the \code{\link{prMapVariable2Name}}
#' 
#' @return \code{matrix} The \code{reordered_groups} with references and the
#'  attribute "var_order" in order to keep track of no. of variables per row.
#' @family printCrudeAndAdjusted functions
#' @inheritParams printCrudeAndAdjustedModel
#' @keywords internal
prCaAddUserReferences <- function (reordered_groups, 
                                   var_order, 
                                   add_references, # Only called when present
                                   add_references_pos, 
                                   reference_zero_effect) {
  for(i in 1:length(var_order)){
    # Add reference if it's not empty
    if (length(add_references) > 1 &&
          is.na(add_references[i]) == FALSE){
      
      within_pos <- ifelse(add_references[i] %in% add_references_pos, 
                           add_references_pos[add_references[i]], 0)
      reordered_groups <- 
        insertRowAndKeepAttr(reordered_groups, 
                             head(var_order[[i]]$location, 1) + within_pos, 
                             rep(c(reference_zero_effect, "ref"), times=2),  
                             rName=add_references[i])
      
      var_order[[i]]$no_rows <- 
        var_order[[i]]$no_rows + 1
      
      for (ii in i:length(var_order)){
        start_pos <- 0
        if (ii > 1){
          start_pos <- tail(var_order[[ii-1]]$location, 1)
        }
        var_order[[ii]]$location <- start_pos + 1:var_order[[ii]]$no_rows
      }
    }
  }
  
  attr(reordered_groups, "var_order") <- var_order
  
  return(reordered_groups)
}

#' Prettify the text
#' 
#' Sets the number of digits, formats the confidence interval and 
#' changes the number of cols into 4 where the upper and lower CI 
#' meet in one string column
#' 
#' @param x The value matrix from getCrudeAndAdjusted 
#' @param ci_lim The limits of the confidence interval
#' @param digits The number of decimal digits to use
#' @param sprintf_ci_str The \code{\link{sprintf}} code for the confidence interval
#' 
#' @return \code{matrix} A string matrix with the values formated
#' @family printCrudeAndAdjusted functions
#' @keywords internal
prCaPrepareCrudeAndAdjusted <- function(x, ci_lim, digits, sprintf_ci_str){
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  # Just to make sure that it gives 1.0 and
  # not 1 if digits = 1, in cases where a
  # adding another decimal that is used
  # since everyone is so hyped about p-val < 0.05
  format_number <- function(x){
    if (length(grep("[0-9]", as.character(x))) == 0)
      return(x)
    
    # Matrix forces all values to be either char or
    # numeric and therefore we need to convert them
    # to numeric
    x <- as.numeric(x)
    
    if (is.wholenumber(x) && x > 100)
      return(as.character(x))
    
    # Remove -0.0 effect
    x <- 
      round(x, digits) %>%
      as.character %>%
      as.numeric
    return(sprintf(sprintf("%%0.%df", digits), x))
  }
  
  # A way to set the min/max of the confidence interval
  # according to the parameters,
  # round to appropriate number of digits and 
  # format into a string as specified in the sprintf string
  format_ci <- function(ci){
    # If there is an NA then we don't know
    # much about the data
    if (any(is.na(ci))){
      return("-")
    }
    
    upper <- max(ci)
    if (upper > max(ci_lim))
      upper <- sprintf("&gt; %s", format_number(max(ci_lim)))
    else
      upper <- format_number(upper)
    
    lower <- min(ci)
    if (lower < min(ci_lim))
      lower <- sprintf("&gt; %s", format_number(min(ci_lim)))
    else
      lower <- format_number(lower)
    
    return(sprintf(sprintf_ci_str, lower, upper))
  }
  
  values <- cbind(
    tapply(x[,1, drop=FALSE], 1:NROW(x), FUN = format_number),
    apply(x[,2:3, drop=FALSE], MARGIN=1, FUN=format_ci),
    tapply(x[,4, drop=FALSE], 1:NROW(x), FUN = format_number),
    apply(x[,5:6, drop=FALSE], MARGIN=1, FUN=format_ci))
  
  colnames(values) <- c(
    colnames(x)[1], 
    sprintf("%s to %s", colnames(x)[2], colnames(x)[3]),
    colnames(x)[4], 
    sprintf("%s to %s", colnames(x)[5], colnames(x)[6]))
  
  rownames(values) <- rownames(x)
  
  return(values)
}    