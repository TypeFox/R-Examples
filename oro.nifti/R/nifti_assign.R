#' @name nifti_assign-methods 
#' @title Methods for Function [<- in Package 'base'
#' 
#' @description Methods for function \code{[<-} in Package 'base'
#' 
#' @name [<--methods
#' @aliases [<--methods [<-,ANY,ANY,ANY,ANY-method [<-,nifti,ANY,ANY,ANY-method
#' [<-,nifti,numeric,numeric,ANY-method [<-,nifti,ANY,missing,ANY-method
#' [<-,nifti,numeric,missing,ANY-method [<-,nifti,missing,missing,array-method
#' @docType methods
#' @section Methods: 
#' \describe{ 
#' \item{x = "nifti", i = "ANY", j = "ANY", value = "ANY"}{Replaces the data at 
#' the provided co-ordinates with the value provided and updates the header.} 
#' \item{x = "nifti", i = "numeric", j = "numeric", value = "ANY"}{Replaces the 
#' data at the provided co-ordinates with the value provided and updates the 
#' header.} 
#' \item{x = "nifti", i = "ANY", j = "missing", value = "ANY"}{Replaces the data
#' row i of the provided nifti object with the value provided and updates the 
#' header.} 
#' \item{x = "nifti", i = "numeric", j = "missing", value = "ANY"}{Replaces the 
#' data row i of the provided nifti object with the value provided and updates 
#' the header.} 
#' \item{x = "nifti", i = "missing", j = "missing", value = "array"}{Replaces 
#' the data of the provided nifti object with the array provided and updates 
#' the header.} 
#' }
#' @keywords methods
NULL
