#' headattribute.class
#'  
#' @param x \code{\link{attribute.class}} input object
#' @export
setMethod(
  f = "head",
  signature = signature(x = 'attribute.class'),
  definition = function(x)
  {
    head(x@results)
  }
)

#' print attribute.class
#'  
#' @param x \code{\link{attribute.class}} input object
#' @export
setMethod(
  f = "print",
  signature = signature(x = 'attribute.class'),
  definition = function(x)
  {
    print(x@results)
  }
)

#' summary attribute.class
#'  
#' @param object \code{\link{attribute.class}} input object
#' @param verbose a logical. If TRUE, additional diagnostics are printed.
#' @export
setMethod(
  f = "summary",
  signature = signature(object = 'attribute.class'),
  definition = function(object, verbose = TRUE, ...){
    results <- slot(object, "results")
    narrow.results <- melt(results, id.vars = "id", variable.name = 'attributes', measure.vars = grep('mean', names(results)))
    mean.mastery.attributes <- ddply(narrow.results, .(attributes), summarize, means = round(mean(value), 3))
    mean.mastery.attributes$attributes <- gsub('_means', '', mean.mastery.attributes$attributes )
    if (verbose){
      cat(sprintf("\nNumber of attributes: %d", nrow(mean.mastery.attributes)))
      cat("\nMean Mastery of each attributes: \n")
      print(as.data.frame(mean.mastery.attributes, row.names = NULL))
    }
    invisible(mean.mastery.attributes)
  }
)

#' plot attribute.class
#'  
#' @param x \code{\link{attribute.class}} input object
#' @param type a string containing either \code{mean} or \code{profile}
#' @export
setMethod(
  f = "plot",
  signature = signature(x = 'attribute.class', y = "missing"),
  definition = function(x, y, type = 'mean', ...)
  {
    results <- slot(x, "results")
    if (type == 'mean'){
      melted.attr <- melt(results, id.vars = "id", measure.vars = grep('mean', names(results)),
                          , value.name = "mean.attr", variable.name = "attr.number")
      means.attr <- ddply(melted.attr, .(attr.number), summarize, mean.attr = mean(mean.attr))
      means.attr$attr.number <- gsub('_means', '', means.attr$attr.number)
      print(ggplot(means.attr, aes(x = attr.number, y = mean.attr, fill = attr.number)) + geom_bar(stat = "identity") +
              scale_fill_discrete(name = "Attribute") +
              ylim(0, 1) + ylab("Mean Mastery Proportion") + xlab("Attribute") + ggtitle("Mean Attribute Mastery"))
    }
    if (type == 'profile'){
      results.attr <- cbind(results$id, results[, grep('mean', names(results))])
      names(results.attr) <- c("id", gsub('_means', '', names(results.attr)[2:ncol(results.attr)]))
      ngroups <- ncol(results.attr)
      PlotSkillMasteryTableplot(results.attr, ngroups, is.max.class = FALSE)
    }
  }
)