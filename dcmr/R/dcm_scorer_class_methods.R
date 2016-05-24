#' summary of dcm.scorer.class
#'  
#' @param object \code{\link{dcm.scorer.class}} input object
#' @param verbose a logical. If TRUE, additional diagnostics are printed
#' @method summary dcm.scorer.class
#' @export
setMethod(
  f = "summary",
  signature = signature(object = 'dcm.scorer.class'),
  definition = function(object, verbose = TRUE, ...){
    attribute.result <- slot(slot(object, 'results'), 'attribute.result')
    attribute.profile.result <- slot(slot(object, 'results'), 'attribute.profile.result')
    inputs <- slot(object, 'inputs')
    if (verbose){
      cat(sprintf('Number of Observations: %d\n', inputs$nobs))
      cat(sprintf('Number of Items: %d\n',inputs$nitems))
      cat(sprintf('Number of Attributes: %d\n', inputs$nattributes))
      cat(sprintf('Number of Attribute profiles: %d\n', inputs$nclass))
      cat(sprintf('Qmatrix:\n')) 
      print(inputs$qmatrix)
      cat(sprintf('\n'))
      summary(attribute.result)
      summary(attribute.profile.result)
    }
    invisible(list(attribute.result = summary(attribute.result, verbose = FALSE),
                   attribute.profile.result = summary(attribute.profile.result, verbose = FALSE)))
  }
)

#' plot of dcm.scorer.class 
#' @param x \code{\link{dcm.scorer.class}} input object
#' @param type a string containing one of the result type
#' @method plot.dcm.scorer.class
#' @export
setMethod(
  f = "plot",
  signature = signature(x = 'dcm.scorer.class', y = 'missing'),
  definition = function(x, y, type = 'attr.means', ...)
  {
    if (type == 'attr.means'){
      attr.results <- slot(slot(x, "results"), "attribute.result")
      plot(attr.results, type = 'mean')
    } else if (type == 'attr.profiles'){
      attr.results <- slot(slot(x, "results"), "attribute.result")
      plot(attr.results, type = 'profile')
    } else if (type == 'attr.profile.means'){
      attr.results <- slot(slot(x, "results"), "attribute.profile.result")
      plot(attr.results, type = 'mean')
    } else if (type == 'attr.profile.profiles'){
      attr.results <- slot(slot(x, "results"), "attribute.profile.result")
      plot(attr.results, type = 'profile')
    }
  }
)

