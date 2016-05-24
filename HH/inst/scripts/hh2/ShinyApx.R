### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/ShinyApx.tex'

###################################################
### code chunk number 1: ShinyApx.tex:8-9
###################################################
library(HH)


###################################################
### code chunk number 2: ShinyApx.tex:26-30
###################################################
## hhcode("shinyNTplot.R", '
NTplot(mean0=8, mean1=8.411, sd=2, n=64, cex.prob=1.3,
      shiny=TRUE)
## ')


###################################################
### code chunk number 3: ShinyApx.tex:55-59
###################################################
## hhcode("shinyBivNormal.R", '
shiny::runApp(system.file("shiny/bivariateNormal",
                          package="HH"))
## ')


###################################################
### code chunk number 4: ShinyApx.tex:81-85
###################################################
## hhcode("shinyBivNormalScatterplot.R", '
shiny::runApp(system.file("shiny/bivariateNormalScatterplot",
                          package="HH"))
## ')


###################################################
### code chunk number 5: ShinyApx.tex:111-115
###################################################
## hhcode("shinyUSage.R", '
shiny::runApp(system.file("shiny/PopulationPyramid",
              package="HH"))
## ')


