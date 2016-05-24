#' Save Plot Object as PDF
#'
#' Internal method that saves plots as PDF.Can be used with all plotting methods in the package. The
#' name of the plot object and it's name to be saved under are passed in and the plot is saved as a PDF.
#'
#' @param title new name for the plot, and name of the PDF file created
#' @param plotName name of the ggplot object
#'
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#'

savePDF = function(title, plotName){
  PDF <-   downloadHandler(
    filename = title,
    content = function(file) {
      ggsave(file, plot = plotName)
    })

  return(PDF)
}


#' Save Plot Object as .RData file
#'
#' Internal method that saves ggplot plots as .RData files.Can be used with all plotting methods in the package. The
#' name of the plot object and it's name to be saved under are passed in and the plot is saved as an RData file.
#'
#' @param title new name for the plot, and name of the RData file created.
#' @param plotName name of the ggplot object.
#'
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#'
savePlot = function(title, plotName){
  PLOT <- downloadHandler(
    filename = title,
    content = function(file) {
      save(plotName, file = file)
    }
  )
}

