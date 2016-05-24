is.PhyloValue <-
function (x) 
{
  inherits(x, "PhyloValue")
}


plot.PhyloValue <- 
function (x, xlab = expression(italic("T")), ylab = NULL, main = NULL, ...) 
{
  Entity <- ""
  # Entity
  if (is.PhyloEntropy(x)) {
    Entity <- "Entropy"
  } else {
    if (is.PhyloDiversity(x)) {
      Entity <- "Diversity"
    }
  }
  
  # ylab
  if (is.null(ylab))
    ylab <- Entity
  
  # main
  if (is.null(main))
    main <- paste(Entity, "along the tree")
  
  graphics::plot(x=c(0, names(x$Cuts)), y=c(x$Cuts, 0), type="s", xlab=xlab, ylab=ylab, main=main, ...)
  graphics::abline(h=x$Total, lty=2)
}


summary.PhyloValue <-
function(object, ...) 
{
  cat(object$Function, "applied to", object$Distribution, "along the tree:", object$Tree, fill=TRUE)
  cat("\nResults are", ifelse(object$Normalized, "normalized", "not normalized"), fill=TRUE)
  cat("\nThe average value is:", object$Total)
  cat("\n\nValues along the tree are:\n")
  print(object$Cuts)
  
  return(invisible(NULL))
}