require(hwriterPlus)
examplesDir <- file.path(system.file(package = 'hwriterPlus'),
                         'examples')
catsSVG <- file.path(examplesDir, "cats.svg")
### Raw html produced
hwriteSVG(catsSVG, page = NULL, height = 600, width = 600,
          center = FALSE, br = TRUE)
