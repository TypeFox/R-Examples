### Unit tests of function logHist

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make


graphicstest.logHist <- function()
{
  ## Purpose: Test the log histogram function logHist
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date: 26 Jan 2010, 16:43

  ## open file for graphical output
  graphicsOutput <- paste(pathReport, "logHist.pdf", sep = "")
  cat("Graphics output in file ", graphicsOutput, "\n")
  pdf(file = graphicsOutput)
  print(pathReport)

  ## generate data
  x <- rgamma(200, 1)

  ## default
  logHist(x)

  ## log histogram only
  logHist(x, htype = "h")
  ## points only, some options
  logHist(x, htype = "p", pch = 20, cex = 2, col = "steelblue")

  dev.off()


  return()
}
