### Unit tests of functions tailPlot, normTailPlotLine, tTailPlotLine
### and gammaTailPlotLine

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make


graphicstest.tailPlot <- function()
{
    ## Purpose: Test the tail plot functions
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: David Scott, Date: 10 Mar 2010, 21:00

    ## open file for graphical output
    graphicsOutput <- paste(pathReport, "tailPlot.pdf", sep = "")
    cat("Graphics output in file ", graphicsOutput, "\n")
    pdf(file = graphicsOutput)
    print(pathReport)

    ## Draw tail plot of some data
    x <- rnorm(100, 1, 2)
    tailPlot(x)
    ## Add normal distribution line
    normTailPlotLine(x, mean = 1, sd = 2)
    ## Add t distribution line
    tTailPlotLine(x, df = 5, lty = 2)
    ## Use fitted values
    normTailPlotLine(x, mean = mean(x), sd = sd(x), lty = 3)

    ## Gamma distribution
    x <- rgamma(100, shape = 1, scale = 1)
    tailPlot(x)
    ## Add gamma distribution line
    gammaTailPlotLine(x, shape = 1, scale = 1)
    ## Left tail example
    tailPlot(x, side = "l")
    ## Add gamma distribution line
    gammaTailPlotLine(x, shape = 1, scale = 1, side = "l")
    ## Log scale on both axes
    tailPlot(x, side = "l", log = "xy")
    ## Add gamma distribution line
    gammaTailPlotLine(x, shape = 1, scale = 1, side = "l")


    ## Add line from a standard distribution with default parameters
    x <- rlnorm(100)
    tailPlot(x)
    tailPlotLine(x, distrFn = "lnorm")

    ## Add line from a distribution with 'param' argument
    require(VarianceGamma)
    param <- c(0,0.5,0,0.5)
    x <- rvg(100, param = param)
    tailPlot(x)
    tailPlotLine(x, distrFn = "vg", param = param)
    dev.off()


    return()
}
