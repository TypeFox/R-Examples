#' Get the metadata of a dataset from the BEFdata portal
#'
#' This function fetches the metadata associated with a dataset on a BEFdata portal.
#' You need to provide the function with an id of the dataset or a path to an eml
#' file you downloaded from the portal manually. You can find the ID in the URL of
#' the dataset page in the portal (e.g ID 6: .../datasets/6).
#'
#' @param dataset This is either a dataset ID or a file path to a eml metadata file.
#'
#' @return A list of metadata. Metadata That doesn't exist is represented as \code{NA}
#' @import XML
#' @export bef.portal.get.metadata bef.get.metadata bef.get.metadata_for bef.portal.get.metadata_for
#' @aliases bef.get.metadata bef.get.metadata_for bef.portal.get.metadata_for

bef.portal.get.metadata <-  bef.get.metadata <- bef.get.metadata_for <- bef.portal.get.metadata_for <- function(dataset) {
  if(is.character(dataset)) {
    object = dataset
  } else {
    object = dataset_url(dataset, "eml", separate_category_columns=TRUE)
  }

  metadata = xmlTreeParse(object, useInternalNodes = T)

  template_dataset = list(title = "//dataset/title",
                  abstract = "//dataset/abstract",
                  publicationDate = "//dataset/pubDate",
                  language = "//dataset/language",
                  creators = list(givenName = "//dataset/creator//givenName",
                                 surName = "//dataset/creator//surName",
                                 electronicMailAddress = "//dataset/creator//electronicMailAddress"),
                  authors = list(givenName = "//creator/individualName/givenName",
                                 surName = "//creator/individualName/surName",
                                 electronicMailAddress = "//creator/electronicMailAddress"),
                  intellectualRights = "//dataset/intellectualRights",
                  distribution = "//dataset/distribution/online/url",
                  keywords = "//keyword",
                  generalTaxonomicCoverage = "//dataset//generalTaxonomicCoverage",
                  samplingDescription = "//dataset//samplingDescription/para",
                  spatial_coverage = list(description = "//geographicCoverage/geographicDescription",
                                             coordinates = list(west = "//geographicCoverage/boundingCoordinates/westBoundingCoordinate",
                                                             east = "//geographicCoverage/boundingCoordinates/eastBoundingCoordinate",
                                                             north = "//geographicCoverage/boundingCoordinates/northBoundingCoordinate",
                                                             south = "//geographicCoverage/boundingCoordinates/southBoundingCoordinate")),
                  temporal_coverage = c(begin = "//temporalCoverage//beginDate",
                                        end = "//temporalCoverage//endDate"),
                  related_material = list(objectName = "//otherEntity/entityName",
                                          dataFormat = "//otherEntity/physical/dataFormat",
                                          distribution = "//otherEntity/physical/distribution")
                  )


  out = rapply(template_dataset, function(x) xmlNodesValue(path=x, doc=metadata), how="replace")

  out$creators = as.data.frame(out$creators, stringsAsFactors=F)
  out$authors = as.data.frame(out$authors, stringsAsFactors=F)
  out$spatial_coverage$coordinates = as.data.frame(out$spatial_coverage$coordinates, stringsAsFactors=F)
  out$related_material = as.data.frame(out$related_material, stringsAsFactors=F)

  attributeList = getNodeSet(metadata, path="//attributeList/attribute")
  column_template = list(header = "./attributeLabel",
                         description = "./attributeDefinition",
                         unit = ".//unit",
                         numericDomain = ".//numericDomain",
                         info = ".//attributeName"
                         )
  columns = lapply(column_template, function(c) {
    sapply(attributeList, function(d) {
      xmlNodesValue(doc=d, path=c)
    })
  })

  columns = as.data.frame(columns, stringsAsFactors=F)

  if (nrow(columns)) {
    columns$unit = as.character(columns$unit)
    # columns$unit[is.na(columns$unit) | columns$unit == "dimensionless"] = "dimensionless"
    out$columns = columns
  }

  return(out)
}


