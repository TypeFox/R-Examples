#
#  sparsebnUtils-messages.R
#  sparsebnUtils
#
#  Created by Bryon Aragam (local) on 3/29/16.
#  Copyright (c) 2014-2016 Bryon Aragam. All rights reserved.
#

#
# PACKAGE SPARSEBNUTILS: Messages
#
#   CONTENTS:
#       input_not_sparsebnData
#       alg_input_data_frame
#
#

#' @name sparsebn-messages
#' @rdname sparsebn-messages
#'
#' @title Messages
#'
#' @description Warning and error messages for use in the \code{sparsebn} family
#'
#' @param data data object.
#' @param count number of missing values.
#' @param pkg package name.
#'
NULL

### User inputs invalid data object
#' @rdname sparsebn-messages
#' @export
input_not_sparsebnData <- function(data){
    sprintf("Input data must be a valid sparsebnData object! <Current type: %s>", class(data))
}

### User inputs data.frame instead of sparsebnData: Implicit coercion will happen
#' @rdname sparsebn-messages
#' @export
alg_input_data_frame <- function(){
    sprintf("Data input as a data.frame: In order to coerce your data to a valid sparsebnData object, we will assume the data is purely observational. In the future, it's best to do this yourself to prevent loss of information.")
}

### Data has missing values
#' @rdname sparsebn-messages
#' @export
has_missing_values <- function(count){
    sprintf("Data contains %d missing values. Presently, all of the methods in this package require complete data. Please impute these missing values before running any of the learning algorithms.", count)
}

### Invalid package requested by user
#' @rdname sparsebn-messages
#' @export
invalid_pkg_specification <- function(){
    sprintf("Incorrect package specified. Must be one of: 'graph', 'igraph', 'network'.")
}

### Required (suggested) package not installed
#' @rdname sparsebn-messages
#' @export
pkg_not_installed <- function(pkg){
    sprintf("The %s package is required in order to use this method. Please install it first.", pkg)
}

### Warning about coercion of objects in the global environment
#' @rdname sparsebn-messages
#' @export
global_coerce_warning <- function(pkg){
    if(!is.null(pkg)){
        sprintf("coerce set to TRUE: All fitted objects will be converted to use objects from the %s package internally.", pkg)
    } else{
        sprintf("coerce set to TRUE: All fitted objects will be reverted back to using default edgeList format.")
    }
}
