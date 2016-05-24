#' print attribute.profile.class
#'  
#' @param x \code{\link{attribute.profile.class}} input object
#' @method print attribute.profile.class
#' @export
setMethod(
  f = "print",
  signature = signature(x = 'attribute.profile.class'),
  definition = function(x)
  {
    print(x@results)
  }
)

#' head attribute.profile.class
#'  
#' @param x \code{\link{attribute.profile.class}} input object
#' @method head attribute.profile.class
#' @export
setMethod(
  f = "head",
  signature = signature(x = 'attribute.profile.class'),
  definition = function(x)
  {
    head(x@results)
  }
)

#' summary attribute.profile.class
#'  
#' @param object \code{\link{attribute.profile.class}} input object
#' @param verbose a logical. If TRUE, additional diagnostics are printed.
#' @method summary attribute.profile.class
#' @export
setMethod(
  f = "summary",
  signature = signature(object = 'attribute.profile.class'),
  definition = function(object, verbose = TRUE, ...){
    results <- slot(object, "results")
    nresults <- nrow(results)
    pmatrix <- slot(object, "attribute.profile.matrix")
    attribute.profile <- ddply(results, .(max.class), summarize, counts = length(max.class),
                               proportion = length(max.class))
    attribute.profile$proportion <- with(attribute.profile, round(proportion/nresults, 3))
    
    attribute.profile <- rename(attribute.profile, c('max.class' = 'attribute.profile'))
    attribute.profile <- cbind(attribute.profile, pmatrix[attribute.profile$attribute.profile,])
    attribute.profile <- attribute.profile[, c('attribute.profile', colnames(pmatrix), 'counts', 'proportion')]
    if (verbose){
      cat(sprintf("\nNumber of attribute profiles: %d", nrow(attribute.profile)))
      cat("\nAttribute profile counts and proportion: \n")
      print(as.data.frame(attribute.profile, row.names = NULL))
    }
    invisible(attribute.profile)
  }
)

#' plot attribute.profile.class
#'  
#' @param x \code{\link{attribute.profile.class}} input object
#' @param type a string containing either \code{mean} or \code{profile}
#' @method plot attribute.profile.class
#' @export
setMethod(
  f = "plot",
  signature = signature(x = 'attribute.profile.class', y = "missing"),
  definition = function(x, y, type = 'mean', ...)
  {
    results <- slot(x, "results")
    if (type == 'mean'){
      pmatrix <- slot(x, "attribute.profile.matrix")
      melted.attr.profile <- melt(results, id.vars = "id", measure.vars = grep('[0-9]+', names(results), )
                                  , value.name = "mean.attr.profile", variable.name = "attr.profile.number")
      means.attr.profile <- ddply(melted.attr.profile, .(attr.profile.number), summarize, mean.attr.profile = mean(mean.attr.profile))
      means.attr.profile$profile.labels <- sapply(as.numeric(levels(means.attr.profile$attr.profile.number)), function(x) paste(pmatrix[x, ], collapse = ","))
      print(ggplot(means.attr.profile, aes(x = attr.profile.number, y = mean.attr.profile, fill = attr.profile.number)) + geom_bar(stat = "identity") +
              scale_fill_discrete(name = "Attribute Profile", labels = means.attr.profile$profile.labels) +
              ylim(0,1) + ylab("Mean Mastery Proportion") + xlab("Attribute Profile") + ggtitle("Mean Attribute Profile Mastery"))
    }
    if (type == 'profile'){
      ngroups <- ncol(results) - 2
      PlotSkillMasteryTableplot(results, ngroups, is.max.class = TRUE)
    }
  }
)