#' Get TSPlib edge weight types.
#'
#' @export
getValidEdgeWeightTypes = function() {
  c("EUC_2D", "EUC_3D", "MAX_2D", "MAX_3D", "MAN_2D", "MAN_3D", "CEIL_2D", "GEO", "ATT", "EXPLICIT")
}
