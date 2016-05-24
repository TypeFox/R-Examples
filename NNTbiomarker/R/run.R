#' run
#'
#' Run a shiny app for this package.
#' @param shinyDir Current options are "shinyElicit" and "shinyCombinePlots". If not provided, a menu of the options is provided.
#' @details The selected shiny app is run. See the vignette Using_the_NNTbiomarker_package for details, and the vignette The_Biomarker_Crisis for an overview.
#'

run = function(shinyDir) {
  if(missing(shinyDir)) {
     shinyDirs = dir(pattern = "^shiny",
               path = file.path(find.package("NNTbiomarker")))
     shinyDir = shinyDirs[menu(shinyDirs)]
  }
  runApp(
    file.path(find.package("NNTbiomarker"), shinyDir)
  )
}

#' runElicit
#'
#' Run a shiny app outlining the process of specifying a design for a biomarker validation study.
#'
#' @seealso run
runElicit = function() run("shinyElicit")

#' runCombinePlots
#'
#' Run a shiny app connecting a visual scale for NNT quantities and a "contra-Bayes" plot for mapping from predictive values to sensitivity/specificity (Bayes theorem in reverse).

#' @seealso run
runCombinePlots = function() run("shinyCombinePlots")


