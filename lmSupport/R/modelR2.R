modelR2 <- function(Model,  Print=TRUE)
{
  Results = modelSummary(Model, Print = FALSE)
  
  if(Print){
    print(Results$call)
    cat(sprintf('\nR-squared:  %.4f, Adjusted R-squared: %.4f\n', Results$r.squared, Results$adj.r.squared), sep='')
    cat(sprintf('F(%.0f,%.0f) = %.4f, p = %.4f\n', Results$fstatistic[2], Results$fstatistic[3], Results$fstatistic[1], pf(Results$fstatistic[1],Results$fstatistic[2], Results$fstatistic[3], lower.tail=FALSE )))
  }
  
  invisible(Results)
}