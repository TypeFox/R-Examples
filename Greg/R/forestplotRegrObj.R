#' Forest plot for multiple models
#' 
#' Plot different model fits with similar variables in order to 
#' compare the model's estimates and confidence intervals. Each
#' model is represented by a separate line on top of eachother
#' and are therefore ideal for comparing different models. This 
#' extra appealing when you have lots of variables included in
#' the models.
#' 
#' @param regr.obj A regression model object. It should be of coxph, crr or glm class.
#'   Warning: The glm is not fully tested.
#' @param skip.variables Which variables to use. The variables should be the 
#'   names of the fit output and not the true output names if you're using
#'   the rowname_translate_function. 
#' @param add.empty_row Add empty rows. This can either be a 
#'   vector or a list. 
#'   When you have a vector the number indicates the row number where 
#'   the empty row should be added, the format is: c(3, 5).
#'   If you give a list you have the option of specifying the name of the
#'   row, the format is: \code{list(list(3, "my rowname"), list(5, "my other rowname"))}.
#'   The rows will be added at the 3rd row and 5th row from the original
#'   position. Ie you don't have take into account that the 5:th row will be
#'   at the 6:th position after adding the 3rd row.
#' @param order.regexps A regexp vector that searches for matches along the original
#'   rownames and reorders according to those. 
#' @param order.addrows If there are ordered groups then often you want empty rows
#'   that separate the different groups. Set this to true if you want to add these
#'   empty rows between groups.
#' @param get_box_size A function for extracting the box sizes
#' 
#' @example inst/examples/forestplotRegrObj_example.R
#' 
#' @importFrom Gmisc insertRowAndKeepAttr
#' 
#' @inheritParams forestplotCombineRegrObj
#' @inheritParams forestplot::forestplot
#' 
#' @family \code{\link[forestplot]{forestplot}} wrappers
#' @rdname forestplotRegrObj
#' @import forestplot
#' @export
forestplotRegrObj <- function(  
  regr.obj, 
  skip.variables,
  add.empty_row,
  order.regexps,
  order.addrows,
  box.default.size,
  rowname.fn,
  xlab,
  xlog,
  exp,
  estimate.txt     = xlab,
  zero,
  get_box_size = fpBoxSize,
  ...)
{
  # Treat always as multiple regression object fits
  # is.list might seem more intuitive but it turns out
  # that objects are lists but lists are not objects
  if (is.object(regr.obj) == TRUE)
    regr.obj = list(regr.obj)
  
  # Initiate some standard values if the user
  # hasn't supplied any
  if (missing(xlab)){
    if (isFitCoxPH(regr.obj[[1]])){
        xlab <- "Hazard Ratio"
    }else if(isFitLogit(regr.obj[[1]])){
      xlab <- "Odds Ratio"
    }
  }
  
  if (missing(estimate.txt)){
    if (isFitCoxPH(regr.obj[[1]])){
      estimate.txt <- "HR"
    }else if(isFitLogit(regr.obj[[1]])){
      estimate.txt <- "OR"
    }
  }
  
  if (missing(xlog)){
    if (isFitCoxPH(regr.obj[[1]]) ||
          isFitLogit(regr.obj[[1]])){
      xlog <- TRUE
      if (missing(zero))
        zero <- 1
      if (missing(exp))
        exp <- TRUE
    }else{
      xlog <- FALSE
      if (missing(zero))
        zero <- 0
      if (missing(exp))
        exp <- FALSE
    }
  }else{
    if (xlog){
      if (missing(zero))
        zero <- 1
      if (missing(exp))
        exp <- TRUE
    }else{
      if (missing(zero))
        zero <- 0
      if (missing(exp))
        exp <- FALSE
    }
  }
  
  
  if (is.null(zero)){
    if (isFitCoxPH(regr.obj[[1]]) ||
      isFitLogit(regr.obj[[1]]))
      zero = 1
    else
      zero = 0
  }
  
  checkValidRegrObject <- function(fit){
    if (!is.object(fit)){
      stop("You must provide a valid model fit of either Cox PH, CRR or GLM type")
    }
    
    if (any(
      c(
        class(fit) %in% "glm",
        class(fit) %in% "crr", 
        class(fit) %in% "coxph")) == FALSE){
      stop("Object not of Cox PH type, Comp. Risk Regression or GLM")
    }
    
    return(TRUE)
  }
  
  # Validate that all fits are OK objects that
  # this function can handle
  for(i in 1:length(regr.obj))
    checkValidRegrObject(regr.obj[[i]])
    
  checkIfCompatibleFits <- function(sd1, sd2){
    checkEqualNames <- function(base_names, cmpr_names){
      # A function for checking that the labels/row names of 
      # the two vectors are the same so that we compare the
      # right elements 
      
      # Check that the vectors contain the same elements
      if (all(cmpr_names %in% base_names) == FALSE){
        missing <- cmpr_names[!cmpr_names %in% base_names]
        extra_elements <- paste(missing, collapse=", ")
        cat(paste("\nThere are more elements in comparing vector, additional elements are:", extra_elements))
        return(FALSE)
      }else if (all(base_names %in% cmpr_names) == FALSE){
        missing <- base_names[!base_names %in% cmpr_names]
        extra_elements <- paste(missing, collapse=", ")
        cat(paste("\nThere are more elements in base vector, additional elements are:", extra_elements))
        return(FALSE)
      }
      
      # Check that the order in which the element appear is
      # the same
      retval <- TRUE
      for(i in 1:length(cmpr_names)){
        if (base_names[i] != cmpr_names[i]){
          cat(paste("\nInvalid order between the two variable name vectors\n",
              "On place", i, "comes",
              base_names[i], "that is not equal to", 
              cmpr_names[i]))
          retval <- FALSE
        }
      }
      
      return(retval)
    }
    
    if (nrow(sd1) != nrow(sd2)){
      stop("The number of variables in the two objects don't match")
    }
    if (checkEqualNames(rownames(sd1), rownames(sd2)) == FALSE){
      stop("The rownames failed to pass the name check! This is probably due to that the variable names differ")
    }
    
    return(TRUE)
  }
  
  # Use the first fit to validate all the other
  # models against
  fit_data <- prGetFpDataFromFit(regr.obj[[1]],
                                 conf.int = 0.95,
                                 exp = exp)
  models_fit_fp_data <- list(fit_data)
  if (length(regr.obj) > 1){
    for(i in 2:length(regr.obj)){
      new_fit_data <- prGetFpDataFromFit(regr.obj[[i]],
        conf.int = 0.95,
        exp = xlog)
      # Check that rownames, number of rows and order matches the first fit
      # The script stops on error
      checkIfCompatibleFits(fit_data, new_fit_data)
      # Add the fit data to the list
      models_fit_fp_data <- append(models_fit_fp_data, list(new_fit_data))
    }
  }
  
  # Change the order according to a regexp vector
  greps <- list()
  if (!missing(order.regexps)){
    if(is.vector(order.regexps) && 
      is.character(order.regexps)){
        for (i in 1:length(models_fit_fp_data)){
          rn <- rownames(models_fit_fp_data[[i]])
          # Use the first one since they all should have the same rownames
          for (expression in order.regexps){
            matches <- grep(expression, rn)
            new_vars <- setdiff(matches, 
                                unlist(greps, use.names = FALSE))
            if (length(new_vars) > 0){
              greps <- append(greps, 
                              list(new_vars))
            }
          }
          
          models_fit_fp_data[[i]] <- 
            models_fit_fp_data[[i]][unlist(greps, use.names = FALSE), ,drop=FALSE]
          
          if (!missing(order.addrows) && length(greps) > 0){
            line_row <- 0
            for (t in 1:(length(greps)-1)) {
              line_row <- line_row + length(greps[[t]])  + 1
              for(i in 1:length(models_fit_fp_data)){
                  dels_fit_fp_data[[i]] <- 
                    insertRowAndKeepAttr(models_fit_fp_data[[i]], 
                                         r=line_row, ,
                                         rName= "* EMPTY ROW *")
              }
            }
          }
        }
        
        no_rows <- sapply(models_fit_fp_data, nrow)
        if (any(no_rows != no_rows[1]))
          stop("Failed to get the same number of rows after ordering the rows",
               " according to the order.regexps/addrows arguments. The number of rows",
               " returned were ", paste(no_rows, collapse=", "), " rows for",
               " the different objects")
    }else{
      stop("Your order.regexps should be of a vector of characters.")
    }
  } 
  
  
  # Add empty rows if specified
  if (!missing(add.empty_row)){
    
    count = 0
    # Walk through the additions
    for(row_nr in add.empty_row){
      # Add to each model
      for(i in 1:length(models_fit_fp_data)){
        if (is.list(row_nr)){
          models_fit_fp_data[[i]] <- 
            insertRowAndKeepAttr(m = models_fit_fp_data[[i]], 
              r = row_nr[[1]]+count,
              rName= row_nr[[2]])
        }else{
          models_fit_fp_data[[i]] <- 
            insertRowAndKeepAttr(m = models_fit_fp_data[[i]], 
              r = row_nr+count,
              rName= "* EMPTY ROW *")
        }
        
      }
      count = count + 1
    }
  }
  
  getVariables2Keep <- function (fit_data, skip.variables){
    rn <- rownames(fit_data)
    #TODO: change to map variable fn
    keep.variables = rep(TRUE, length=length(rn))
    # Find those variables that are set to be skipped
    if (!missing(skip.variables)){
      for(sk in skip.variables){
        keep.variables = !(sk == rn) & keep.variables
      }
    }
    return(keep.variables)
  }
  
  geRownames <- function(fit_data, 
                         rn_translate_fn,
                         skip.variables){
    if (is.matrix(fit_data) == FALSE){
      stop("You must provide a proper matrix from the getForsetplotData functions")
    }
    
    keep.variables <- getVariables2Keep(fit_data,
                                        skip.variables = skip.variables)
    
    # Prepare the row names
    rn <- rownames(fit_data)
    
    keep.variables <- getVariables2Keep(fit_data,
                                        skip.variables = skip.variables)
    
    
    # Modify the row names so that you can have
    # true expressions if you so wish
    rn <- as.list(rn[keep.variables])
    if (!missing(rn_translate_fn)){
      if (is.character(rn_translate_fn))
        rn_translate_fn <- get(rn_translate_fn)
      
      for (i in 1:length(rn)){
        rn[[i]] <- rn_translate_fn(rn[[i]]) 
      }
    }
    
    # Add at the top of the list the header
    return(append("Variable", rn))
  }
  
  col1 <- geRownames(models_fit_fp_data[[1]], 
                     rn_translate_fn = rowname.fn,
                     skip.variables = skip.variables)
  
  # The top is the header and should be bold
  is.summary <- c(TRUE, rep(FALSE, length(col1)-1))
  
  # Select the variables to display
  # and add empty top row
  for(i in 1:length(models_fit_fp_data)){
    keep.variables <- getVariables2Keep(models_fit_fp_data[[i]],
                                        skip.variables = skip.variables)
    
    models_fit_fp_data[[i]] <- models_fit_fp_data[[i]][keep.variables, ,drop=FALSE]
    models_fit_fp_data[[i]] <- 
      insertRowAndKeepAttr(models_fit_fp_data[[i]], 
        r=1)
  }
  
  # Initiate the plot data
  t.coef <- models_fit_fp_data[[1]][, "beta"]
  t.low <- models_fit_fp_data[[1]][, "low"]
  t.high <- models_fit_fp_data[[1]][, "high"]
  
  # Make the box smaller if there are many
  # models that share the same space
  if (missing(box.default.size)){
    box.default.size <- .75/length(models_fit_fp_data)
  }
  
  variable_count <- nrow(models_fit_fp_data[[1]])-1
  b_size <- get_box_size(p_values = models_fit_fp_data[[1]][, "p_val"],
    variable_count = variable_count,
    box.default.size = box.default.size)
  

  if (length(models_fit_fp_data) > 1){
    # Create a matrix for the box etc
    for(i in 2:length(models_fit_fp_data)){
      # Addon the main variables
      t.coef <- cbind(t.coef, models_fit_fp_data[[i]][, "beta"])
      t.low <- cbind(t.low, models_fit_fp_data[[i]][, "low"])
      t.high <- cbind(t.high, models_fit_fp_data[[i]][, "high"])
      
      # Addon the box size for this fit
      b_size <- cbind(b_size, 
        get_box_size(p_values = models_fit_fp_data[[i]][, "p_val"],
          variable_count = variable_count,
          box.default.size = box.default.size))
    }
    
    # There are only the variable names that should be displayed
    rn <- list(col1)
  }else{	
    # To complex for most
    #col2 <- list(expression(plain(e)^beta));
    col2 <- list(estimate.txt);
    
    # TODO: should probably use options("digits") but
    # it defaults to 7 why this might be more annoying
    # than helpful
    # The first element is NA and should not be included 
    # as it's HR instead
    for (coef in models_fit_fp_data[[1]][, "beta"][-1])
      col2 <- append(col2, ifelse(is.na(coef), "", sprintf(" %.2f ", coef)))
    
    # The first is the names and the second the hazard rate
    rn <- list(col1, col2)
  }
    
  forestplot(rn, 
              mean                 = t.coef, 
              lower                = t.low, 
              upper                = t.high,
              boxsize              = b_size,
              xlab                 = xlab,
              xlog                 = xlog,
              is.summary           = is.summary,
              zero                 = zero,
              ...)
  
}

#' @param p_values The p-values that will work as the foundation for the box size
#' @param variable_count The number of variables
#' @param box.default.size The default box size
#' @param significant Level of significance .05
#' @rdname forestplotRegrObj
#' @export
fpBoxSize <- function(p_values, 
         variable_count, 
         box.default.size, 
         significant = .05){
  b_size <- c(NA, rep(box.default.size, variable_count))
  b_size[p_values < significant] <- box.default.size*1.5
  return(b_size)
}
