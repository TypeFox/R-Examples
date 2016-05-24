# Functions to support variable binning

# Contents:
#   - function to bin a vector into equal height, equal width, or custom bins
#   - function to take the columns from a dataframe and replace them with their
#     binned values; as well can optionally bin an additional "test" set into
#     those same values
#   - function to take a dataframe of binned data and return a list containing
#     the unique cut points for each bin


#==============================================================================#

#' Get the cutpoints from a single factor vector.
#'
#' get_vector_cutpoints returns a numeric vector giving the unique
#' cutpoints of a variable that has been discretized using
#' vector_bin- more generally, using \code{\link[base]{cut}} and any
#' functions that depend on \code{\link[base]{cut}}
#'
#' @param v vector to get cutpoints from
#' @return a vector containing the unique cutpoints in v
#' @details
#' This function is provided for convienience, and is built to work with the
#' exact format for discretized variables that is used by the cut family. Hence
#' it will work for cut, cut_number/cut_interval, and any of the binning
#' functions from modellingTools, but it is not guaranteed to work for
#' arbitrary factors with numeric levels
#' @family discretization
#' @seealso \code{\link[base]{cut}}, \code{\link[ggplot2]{cut_number}},
#' \code{\link[ggplot2]{cut_interval}}, \code{\link{vector_bin}}
#' @examples
#' x <- cut(rnorm(100),c(-1,0,1))
#' get_vector_cutpoints(x) # -1, 0, 1
#' @export
#' @import magrittr

get_vector_cutpoints <- function(v) {
  if (is.factor(v)) {
    lv <- levels(v)
  } else {
    lv <- v
  }

  cut_points <- lv %>%
                stringr::str_split(",") %>%
                #unlist() %>%
                stringr::str_extract_all("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?") %>%
                unlist() %>%
                unique() %>%
                as.numeric() %>%
                sort()



  return(cut_points)
}


#' Get the unique cutpoints of each appropriate column in a tbl.
#'
#' This function takes a dataframe where any number of columns have been binned
#' into factors using cut/vector_bin and returns a named list containing the
#' cutpoints for each variable.
#' This is useful for, for example, binning a new dataset into the same bins
#' as an older dataset- such as when making predictions on test data
#'
#' @param d a tbl
#' @param parallel logical. If TRUE, parallel foreach used. Must register
#' parallel beforehand. Default false
#' @return a named list containing one vector for each factor type variable.
#' Each vector contains the unique cut points of that variable
#' @family discretization
#' @seealso \code{\link{get_vector_cutpoints}}, \code{\link{simple_bin}}
#' @details
#' This function essentially calls \code{\link{get_vector_cutpoints}} on all
#' factor type columns of \code{d}. It is meant to be used to provide an output
#' format that works with the \code{bins} argument of \code{\link{simple_bin}},
#' for the purpose of defining cutpoints based on one dataset and then applying
#' them to other datasets. The basic functionality of binning on training data
#' and applying those bins to test data is built right in to
#' \code{\link{simple_bin}}, but this function allows the user total
#' flexibility.
#' @examples
#' x <- dplyr::data_frame(v1 = cut(rnorm(100),c(-1,0,1)),
#'                        v2 = cut(rnorm(100),c(-.5,0,.5)))
#' binned_data_cutpoints(x)
#' binned_data_cutpoints(x,parallel = TRUE)
#' @export
#' @import foreach

binned_data_cutpoints <- function(d,parallel = FALSE) {

  # Which variables are factors?
  f_list <- c()
  for (i in 1:ncol(d)) {
    if (is.factor(column_vector(d,i))) {
      f_list <- c(f_list,colnames(d)[i])
    }
  }

  # Get the cutpoints
  if (parallel) {
    cutpoints <- foreach::foreach(cl = f_list,
                         .final = function(x) stats::setNames(x,f_list),
                         .export = c("get_vector_cutpoints",
                                     "column_vector",
                                     "cl"),
                         .packages = c("stringr")) %dopar% {
                            get_vector_cutpoints(column_vector(d,cl))
                         }

  } else {
    cutpoints <- foreach::foreach(cl = f_list,
                         .final = function(x) stats::setNames(x,f_list)) %do% {
                           get_vector_cutpoints(column_vector(d,cl))
                         }
  }

  return(cutpoints)
}

#' Bin a vector into equal height, equal width, or custom bins
#'
#' This function essentially calls \code{\link[base]{cut}}/
#' \code{\link[ggplot2]{cut_interval}}/\code{\link[ggplot2]{cut_number}},
#' depending on the value of bins and type. The one major difference is in the
#' treatment of missing values; those functions return NA, while vector_bin has
#' the default option of returning a bin for the missing values
#'
#' @param x vector of numeric data to bin
#' @param bins numeric vector. If length 1, then this is taken to be the number of
#'         desired bins, computed according to "type". If length > 1, this is taken
#'         to be the actual cutpoints desired
#' @param type character, equal to "height" or "width". Only used if length(bins) == 1.
#'         If "height", then bins are computed to have roughly equal numbers of
#'         observations; else, bins are computed to be of roughly equal width
#' @param na_include logical. If TRUE, then a bin labelled "missing" will be included
#'               in the output. Else NA values are dropped
#' @return the input vector x, with values replaced by the appropriate bins.
#'         Type also changed to factor
#' @family discretization
#' @seealso \code{\link[base]{cut}}, \code{\link[ggplot2]{cut_number}},
#' \code{\link[ggplot2]{cut_interval}},
#' @examples
#' x <- rnorm(100)
#' y <- x; y[sample(1:100,20)] <- NA
#'
#' cut(x,c(-1,0,1))
#' vector_bin(x,bins = c(-1,0,1))
#' cut(y,c(-1,0,1))
#' vector_bin(y,bins = c(-1,0,1))
#' vector_bin(y,bins = c(-1,0,1),na_include = FALSE)
#'
#' ggplot2::cut_number(x,5)
#' vector_bin(x,5)
#'
#' ggplot2::cut_interval(x,5)
#' vector_bin(x,5,type = "width")
#' @export
#' @import magrittr

vector_bin <- function(x,bins,type = "height",na_include = TRUE) {

  if (!is.numeric(x)) return(x)

  if (length(bins) > 1) {
    binned_x <- try(as.character(cut(x,
                                     bins,
                                     right = FALSE,
                                     include.lowest = TRUE)),silent = TRUE)
  }
  else if (type == "height") {
    if (bins >= length(unique(x))) return(x)
    binned_x <- try(as.character(ggplot2::cut_number(x,bins,right = FALSE)),silent = TRUE)
  }
  else {
    if (bins >= length(unique(x))) return(x)
    binned_x <- try(as.character(ggplot2::cut_interval(x,bins,right = FALSE)),silent = TRUE)
  }

  if (class(binned_x) == "try-error") {
    print("Could not bin x")
    return(x)
  }

  if (na_include) {
    binned_x[is.na(binned_x)] <- "Missing"
  }
  else {
    binned_x <- binned_x[!is.na(binned_x)]
  }

  # Replace all open brackets, replace with closed
  binned_x <- binned_x %>%
              stringr::str_replace_all("\\(","\\[") %>%
              stringr::str_replace_all("\\)","\\]")

  # Removed factor class to deal with NAs: re-add it now
  binned_x <- factor(binned_x)

  return(binned_x)

}

#' Discretize variables in your training and test datasets
#'
#' Function to apply simple equal-width or equal-height binning to columns of a
#' training dataset, and then optionally bin the columns of a test set into bins
#' with the appropriate cutpoints
#'
#' @param train training set
#' @param test test set
#' @param exclude_vars variables to exclude (e.g. the target, or the row ID)
#' @param include_vars if you only want certain variables binned, you may specify them
#'                 directly instead of excluding all other variables
#' @param bins single number specifying the number of bins to create on each variable,
#'         or a named list specifying cut-points for each variable
#' @param type if bins is given as a number, then this determines whether to create
#'         bins with equal number of observations ("height") or of equal width
#'         ("width")
#' @param na_include logical. Give missing values their own bin?
#' @return if test is not NULL, a list containing two tbl_df objects, with appropriate
#'     columns replaced by their binned values and all other columns unchanged
#'     if test is NULL, returns the training set portion of the list
#'
#' @family discretization
#' @seealso \code{\link{vector_bin}}, \code{\link{get_vector_cutpoints}}
#' @details
#' This function was built as a convenience, to automate the process of binning
#' continuous variables into disrete levels, and also to provide a simple,
#' interpretible, unambiguous method of dealing with missing values in data
#' science problems.
#' @export
#' @import foreach magrittr

simple_bin <- function(train,
                       test = NULL,
                       exclude_vars = NULL,
                       include_vars = NULL,
                       bins,
                       type = "height",
                       na_include = TRUE) {

  if (length(exclude_vars) > 0 && length(include_vars > 0)) {
    stop("Cannot specify both include and exclude lists")
  }

  index <- c()

  # If bins is a list, then bin only the variables included
  # else, bin everything not in exclude_list, OR that is in include_list

  if (is.list(bins)) {
    index <- names(bins)
  } else if (length(exclude_vars > 0)) {
    for (i in 1:ncol(train)) {
      if (!(colnames(train)[i] %in% exclude_vars)) {
        index <- c(index,colnames(train)[i])
      }
    }
  } else if (length(include_vars) > 0) {
    for (i in 1:ncol(train)) {
      if (colnames(train)[i] %in% include_vars) {
        index <- c(index,colnames(train)[i])
      }
    }
  } else {
    index <- colnames(train)
  }

  # Initialize the output data frame
  # UPDATE: depricated, replaced with foreach()
#   if (colnames(train)[1] %in% index) {
#     nm <- colnames(train)[1]
#     if (is.list(bins)) {
#       tmp_bins <- bins[[1]]
#     } else {
#       tmp_bins <- bins
#     }
#     binned_train <- data_frame(var1 = vector_bin(column_vector(train,nm),tmp_bins,type,na_include))
#   } else {
#     binned_train <- data_frame(var1 = column_vector(train,1))
#   }

  # Bin variables on training set
  k <- 1
  binned_train <- foreach::foreach(i = 1:ncol(train),
                         .combine = dplyr::bind_cols,
                         .multicombine = TRUE,
                         .final = function(x) {
                           x <- dplyr::tbl_df(as.data.frame(x))
                           x <- stats::setNames(x,colnames(train))
                           return(x)
                           },
                         .export = c("column_vector",
                                     "vector_bin"),
                         .packages = c("dplyr")) %do% {
    nm <- colnames(train)[i]
    if (nm %in% index) {
      if (is.list(bins)) {
        tmp_bins <- bins[[nm]]
      } else {
        tmp_bins <- bins
      }
      dplyr::data_frame(var1 = vector_bin(column_vector(train,nm),
                                          tmp_bins,
                                          type,
                                          na_include))
    }
    else {
      dplyr::data_frame(var1 = column_vector(train,nm))
    }
  }

  colnames(binned_train) <- colnames(train)

  # Now bin on the test set, if present
  if (!is.null(test)) {
    # Get the cutpoints from the binned training set
    train_cutpoints <- binned_data_cutpoints(binned_train,parallel = TRUE)

    # Recursively apply the parent function, with train = test and test ~ NULL
    binned_test <- simple_bin(train = test,test = NULL,
                              bins = train_cutpoints,
                              exclude_vars = exclude_vars)

    return(list(train = binned_train,
                test = binned_test))
  } else {
    # If no test set, returned the binned training set
    return(binned_train)
  }

}


#' Create a usable model matrix from a data frame containing a mix of
#' continuous and categorical variables
#'
#' This function takes your dataframe of input variables and returns a new
#' dataframe (or matrix) with the categorical variables replaced by dummy
#' variables, using \code{\link[stats]{model.matrix}}
#'
#' @param dat a tbl
#' @param id character, naming the variable in dat which serves as the unique
#'            row identifier. If blank, will be created
#' @param matrix_out logical. Should the result be a matrix (\code{TRUE}), suitable
#'            for input into many modelling functions, or should the result be
#'            a tbl (\code{FALSE}), suitible for inspection and further analysis?
#'            Default \code{TRUE}
#' @param parallel logical. If \code{TRUE}, parallel
#'            \code{\link[foreach]{foreach}} is used to compute on each
#'            variable. Must register a parallel backend first. Default
#'            \code{FALSE}.
#' @return a matrix or a tbl, consisting of dummy columns with 0/1 indicators
#'         of membership in each factor level for each factor variable, and all
#'         other input variables unchanged.
#'
#' @details
#' The function will only alter variables which are type \code{factor}. Contrary
#' to how it may sound, this actually offers the user greater flexibility, for
#' two reasons: it allows you to keep character type variables intact, and it
#' forces you to think about the levels of each factor variable rather than
#' picking them straight from the input data
#'
#' @examples
#' x <- simple_bin(iris,bins = 3)
#' create_model_matrix(x)
#' create_model_matrix(x,matrix_out = FALSE)
#'
#' @export
#' @import magrittr foreach
create_model_matrix <- function(dat,
                                id = c(),
                                matrix_out = TRUE,
                                parallel = FALSE) {

  # If no provided row id, create one
  if (length(id) == 0) {
    dat %<>% dplyr::mutate(id = 1:nrow(dat))
    id <- c("id")
  }

  # Create a vector of names of variables to include
  include <- foreach::foreach(i = 1:ncol(dat),
                                  .combine = c,
                                  .export = c("i")) %do% {
                                    if (is.factor(column_vector(dat,i))) {
                                          colnames(dat)[i]
                                        }
                                      }

  # Create a list of data frames, each with one widened variable
  # plus a column for the id

  if (parallel) {
    wide <- foreach::foreach(nm = colnames(dat),
                             .combine = dplyr::bind_cols,
                             .final = function(x) {
                               stats::setNames(x,include)
                               return(x)
                             },
                             .multicombine = TRUE,
                             .packages = c("dplyr","modellingTools"),
                             .export = c("nm")) %dopar% {
                               ff <- stats::formula("~ + -1 + " + nm)
                               stats::model.matrix(ff,dat) %>%
                                 as.data.frame() %>%
                                 dplyr::tbl_df()
                             }
  } else {
    wide <- foreach::foreach(nm = include,
                             .combine = dplyr::bind_cols,
                             .final = function(x) {
                               stats::setNames(x,include)
                               return(x)
                             },
                             .multicombine = TRUE,
                             .export = c("nm")) %do% {
                               ff <- stats::formula("~ -1 + " + nm)
                               stats::model.matrix(ff,dat) %>%
                                 as.data.frame() %>%
                                 dplyr::tbl_df()
                             }
  }


  if (matrix_out) {
    # If using as input to a model, remove the id and target columns
    wide %<>% as.matrix()
    return(wide)
  } else {
    wide %<>% dplyr::bind_cols(dat[ ,id])
    return(wide)
  }

}

















