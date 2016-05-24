#' Metadata Prefixes
#'
#' These three options can be used in \code{dc_oai_*()} functions as input to
#' the \code{prefix} parameter.
#'
#' \itemize{
#'  \item oai_dc - OAI Dublin Core - As a minimum requirement for OAI-PMH
#'  compliance, metadata must be made available in the OAI Dublin Core format.
#'  For more information please see the OAI-PMH web site.
#'  \item oai_datacite - OAI DataCite - This metadata format has been specifically
#'  established for the dissemination of DataCite records using OAI-PMH. In addition
#'  to the original DataCite metadata, this format contains several other elements
#'  describing the version of the metadata, whether it is of reference quality, and
#'  the registering datacentre. For more information about this format and its
#'  schema please see the Datacite OAI schema web site.
#'  \item datacite - DataCite Direct - This metadata format contains only the
#'  original DataCite metadata without additions or alterations. Because there are
#'  multiple versions of DataCite metadata in the MDS, there is no one schema that
#'  they will all adhere to. Therefore the schema for this format does not exist
#'  and metadata will not validate against it. Please note that this format is not
#'  OAI-PMH version 2.0 compliant for the previously stated reasons.
#' }
#'
#' @name Prefixes
NULL
