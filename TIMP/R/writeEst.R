"writeEst" <- function (multimodel, multitheta, plotoptions) 
{
  if (length(plotoptions@makeps) != 0) 
            filename <- paste(plotoptions@makeps, "_paramEst.txt", sep = "")
    else
      filename <- "currParamEst"	
  fileparam <- file(filename, "w")
  
  cat(plotoptions@title, "\n\n", file = fileparam)

  onls <- multimodel@fit@nlsres$onls
  if(plotoptions@algorithm=="nls") {
    cat("Sum square error:", onls$m$deviance(), "\n\n", 
        file = fileparam)
  }
  if(plotoptions@algorithm=="nls.lm") {
    cat("Sum square error:", onls$dev, "\n\n", 
        file = fileparam)
  }
  if(plotoptions@algorithm=="optim") {
    cat("Error:", onls$value, "\n\n", 
        file = fileparam)
  }
  if(plotoptions@sumnls && (plotoptions@algorithm=="nls.lm" ||
                            plotoptions@algorithm=="nls")) {
    s <- summary(onls, multimodel)
    cat("Residual standard error:", s$sigma, "on", s$df[2], "degrees
    of freedom\n\n\n", file = fileparam)
  }
  parEst(list(currModel=multimodel, currTheta=multitheta), file=fileparam)
  
  close(fileparam)
}


