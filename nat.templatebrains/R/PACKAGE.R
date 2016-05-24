#' NeuroAnatomy Toolbox add-on package for handling template brains
#'
#' This package provides a class \code{templatebrain} that stores key
#' information about reference brains along with helper functions to simplify
#' transformation of data between template brains (a.k.a bridging) and mirroring
#' of data within a template brain.
#'
#' @section Helper functions: Easy-to-use functions for transforming data from
#'   one template brain to another, displaying slices alongside 3D data, etc.
#'   are provided. See especially \code{\link{xform_brain}},
#'   \code{\link{mirror_brain}} and \code{\link{plot3d.templatebrain}}.
#'
#' @section Package options:
#'
#'   \itemize{
#'
#'   \item options('nat.templatebrains.regdirs') specifies a character vector of
#'   directories containing registrations.
#'
#'   }
#'
#'   It is strongly recommended \itemize{
#'
#'   \item to use the \code{\link{download_reg_repo}} to download
#'   registration folders where possible. This will automatically add the new
#'   folders to this option both in the current session and on package startup.
#'
#'   \item to use \code{\link{add_reg_folders}} to add additional directories
#'   to the search list.
#'
#'   }
#'
#'   Note that registration directories will be searched in the order that they
#'   are listed. It is therefore strongly recommended that individual
#'   registrations have globally unique names.
#'
#' @name nat.templatebrains-package
#' @aliases nat.templatebrains
#' @seealso \code{\link[nat]{nat}}
#' @docType package
#' @keywords package registration template
#' @examples
#' \dontrun{
#' #
#'
#' # manually add a new directory containing registrations to the search list
#' options(nat.templatebrains.regdirs=union(
#'   getOption('nat.templatebrains.regdirs'), "/my/new/path"))
#'
#'
#' ## same, but override any built-in registration by putitng the new path
#' # at the first position in the search list
#' options(nat.templatebrains.regdirs=union("/my/new/path"),
#'   getOption('nat.templatebrains.regdirs'))
#' }
NULL

#' Sample template brain: FlyCircuit Whole Brain
#'
#' This is a sample template brain for testing purposes which is equivalent to
#' the FCWB template brain defined by the \code{nat.flybrains}, which should be
#' considered the canonical version.
#' @name FCWB.demo
#' @docType data
NULL
