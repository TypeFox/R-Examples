### Unit tests of function dghyp
### By Joyce Li, 9/5/2010

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make


### Graphical Test for Generalized Hyperbolic Distribution
graphicstest.dhyperb <- function()
{
  ##Purpose: Test the dghyp function with small range parameters
  n <- 1000
  data(ghypParam)
  smallparam <- ghypSmallShape
  smallmu <- smallparam[, 1]
  smalldelta <- smallparam[, 2]
  smallalpha <- smallparam[, 3]
  smallbeta <- smallparam[, 4]
  smalllambda <- smallparam[, 5]
  sden <- 0
  ## open file for graphical output
  graphicsOutput <- paste(pathReport, "dhyperb.pdf", sep = "")
  cat("Graphics output in file ", graphicsOutput, "\n")
  pdf(file = graphicsOutput, height = 7,width = 10)
  par(mfrow = c(1, 2), oma = c(5, 5, 5, 5))
  for (i in 1 : nrow(smallparam))
  {
      x <- rghyp(n, param = smallparam[i, ])
      sde <- density(x, bw = 0.1)$y
      sden <- dghyp(x,param = smallparam[i, ])
      hist(x, freq = FALSE, breaks = 20, ylim = c(0, max(sden, sde)),
           main = "", xlab = "sample")
      mtext(expression(bold("Graph Test of dghyp")), line = 3.5, cex = 1.15)
      mtext(bquote(paste(lambda ==.(smalllambda[i]), ",",
                         alpha ==.(smallalpha[i]), ",",
                         beta ==.(smallbeta[i]), ",",
                         delta ==.(smalldelta[i]), ",",
                         mu ==.(smallmu[i]),sep = "")),
            line = 2.25, cex = 1.15)
      curve(dghyp(x, param = smallparam[i, ]), add = TRUE, col = "red")
      logHist(x, main = "", breaks = 20, htype = "h")
      mtext(expression(bold("Log Graph Test of dghyp")),
            line = 3.5, cex = 1.15)
      mtext(bquote(paste(lambda ==.(smalllambda[i]), ",",
                         alpha ==.(smallalpha[i]), ",",
                         beta ==.(smallbeta[i]), ",",
                         delta ==.(smalldelta[i]), ",",
                         mu ==.(smallmu[i]),sep = "")),
            line = 2.25, cex = 1.15)
      curve(log(dghyp(x, param = smallparam[i, ])), add = TRUE, col = "red")
      i <- i + 1
  }



  ##Purpose: Test the dghyp function with large range parameters
  largeparam <- ghypLargeShape
  largemu <- largeparam[, 1]
  largedelta <- largeparam[, 2]
  largealpha <- largeparam[, 3]
  largebeta <- largeparam[, 4]
  largelambda <- largeparam[, 5]
  lden <- 0
  par(mfrow = c(1, 2), oma = c(5, 5, 5, 5))
  for (i in 1 : nrow(largeparam))
  {
      x <- rghyp(n, param = largeparam[i, ])
      lde <- density(x, bw = 0.1)$y
      lden <- dghyp(x, param = largeparam[i, ])
      hist(x, freq = FALSE, breaks = 20, ylim = c(0, max(lden, lde)),
           main = "", xlab = "sample")
      mtext(expression(bold("Graph Test of dghyp")),
            line = 3.5, cex = 1.15)
      mtext(bquote(paste(lambda ==.(largelambda[i]), ",",
                     alpha ==.(largealpha[i]), ",",
                     beta ==.(largebeta[i]), ",",
                     delta ==.(largedelta[i]), ",",
                     mu ==.(largemu[i]),sep = "")),
                     line = 2.25, cex = 1.15)
      curve(dghyp(x, param = largeparam[i, ]), add = TRUE, col = "red")
      logHist(x, main = "", breaks = 20, htype = "h")
      mtext(expression(bold("Log Graph Test of dghyp")),
            line = 3.5, cex = 1.15)
      mtext(bquote(paste(lambda ==.(largelambda[i]),",",
                     alpha ==.(largealpha[i]),",",
                     beta ==.(largebeta[i]),",",
                     delta ==.(largedelta[i]),",",
                     mu ==.(largemu[i]),sep = "")),
                     line = 2.25, cex = 1.15)
       curve(log(dghyp(x,param = largeparam[i, ])), add = TRUE, col = "red")
       i <- i + 1
  }

  dev.off()
  
  return()
}



