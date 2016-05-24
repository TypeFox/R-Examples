#
#  sparsebnUtils-functions.R
#  sparsebnUtils
#
#  Created by Bryon Aragam (local) on 1/22/16.
#  Copyright (c) 2014-2016 Bryon Aragam. All rights reserved.
#

#
# PACKAGE SPARSEBNUTILS: Functions
#
#   CONTENTS:
#     check_if_matrix
#     check_if_data_matrix
#     check_if_complete_data
#     count_nas
#     list_classes
#     check_list_class
#     col_classes
#     cor_vector
#

#' @name sparsebn-functions
#' @rdname sparsebn-functions
#'
#' @param x a compatible object.
#' @param m a \code{matrix}.
#' @param df a \code{data.frame}.
#' @param li a \code{list}.
#' @param check.class \code{character} class name to compare against.
#' @param check.names \code{character} names to compare against.
#' @param X a matrix.
#' @param string a \code{character} string.
#'
#' @title Utility functions
#'
#' @description Various utility functions for packages in the \code{sparsebn} family
#'
NULL

# Check if an object is EITHER matrix or Matrix object
#' @rdname sparsebn-functions
#' @export
check_if_matrix <- function(m){
    is.matrix(m) || inherits(m, "Matrix")
} # END .CHECK_IF_MATRIX

# Check if an object is a valid dataset
#' @rdname sparsebn-functions
#' @export
check_if_data_matrix <- function(df){
    is.data.frame(df) || check_if_matrix(df)
} # END .CHECK_IF_DATA_MATRIX

# Check if a dataset contains missing data
#' @rdname sparsebn-functions
#' @export
check_if_complete_data <- function(df){
    (count_nas(df) == 0)
} # END .CHECK_IF_COMPLETE_DATA

# Count missing values in a matrix or data.frame
#' @rdname sparsebn-functions
#' @export
count_nas <- function(df){
    if( !check_if_data_matrix(df)){
        stop("Input must be a data.frame or a matrix!")
    }

    sum(is.na(df))
} # END .COUNT_NAS

# Return the types for each element in a list
#' @rdname sparsebn-functions
#' @export
list_classes <- function(li){
    unlist(lapply(li, class))
} # END .LIST_CLASSES

# Return the number of levels for each column in a data.frame
#' @rdname sparsebn-functions
#' @export
auto_count_levels <- function(df){
    if( !check_if_data_matrix(df)){
        stop("Input must be a data.frame or a matrix!")
    }

    if(!is.data.frame(df)) df <- data.frame(df)
    lapply(df, function(x) length(unique(x)))
} # END .COUNT_NAS

# Return TRUE if every element of a list inherits check.class, FALSE otherwise
#' @rdname sparsebn-functions
#' @export
check_list_class <- function(li, check.class){
    if(length(li) == 0){
        warning("List contains no elements!")

        TRUE # default to true if empty
    }

    all(unlist(lapply(li, function(x) inherits(x, check.class))))
} # END .CHECK_LIST_CLASS

# Return TRUE if names(list) matches check.names, FALSE otherwise
#' @rdname sparsebn-functions
#' @export
check_list_names <- function(li, check.names){
    if(length(li) == 0){
        warning("List contains no elements!")

        TRUE # default to true if empty
    }

    (length(li) == length(check.names)) && (names(li) == check.names)
} # END .CHECK_LIST_NAMES

# Output the class of each column in X, return as a character vector
#' @rdname sparsebn-functions
#' @export
col_classes <- function(X){
    if( !check_if_data_matrix(X)){
        stop("Input must be a data.frame or a matrix!")
    }

    apply(X, 2, class)
} # END .COL_CLASSES

# Compute the correlation matrix of a dataset, and return the unduplicated elements (i.e. upper-triangular portions) as a vector
#  Used as the primary "carrier of information" in ccdr since the algorithms only depends on pairwise correlations
#' @rdname sparsebn-functions
#' @export
cor_vector <- function(X){
    check.numeric <- (col_classes(X) != "numeric")
    if( any(check.numeric)){
        not.numeric <- which(check.numeric)
        stop(paste0("Input columns must be numeric! Columns ", paste(not.numeric, collapse = ", "), " are non-numeric."))
    }

    if( any(dim(X) < 2)){
        stop("Input must have at least 2 rows and columns!") # 2-8-15: Why do we check this here?
    }

    cors <- stats::cor(X)
    cors <- cors[upper.tri(cors, diag = TRUE)]

    cors
} # END COR_VECTOR

# Utility to capitalize the first letter in a string
#  Borrowed verbatim from the 'Hmisc' package
#' @rdname sparsebn-functions
#' @export
capitalize <- function(string) {
    capped <- grep("^[^A-Z]*$", string, perl = TRUE)
    substr(string[capped], 1, 1) <- toupper(substr(string[capped],
        1, 1))
    return(string)
} # END CAPITALIZE
