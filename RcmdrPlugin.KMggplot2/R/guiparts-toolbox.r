#' Tool Box Subclass
#'
#' \code{toolbox} class is a subclass for GUI tool box frame.
#'
#' This class is a subclass which make GUI tool box frame.
#'
#' @section Fields:
#' \describe{
#' \item{\code{frame}: }{\code{tkwin} class object; parent of widget window.} 
#' \item{\code{length}: }{Integer; number of grids.} 
#' \item{\code{back_list}: }{List of \code{tkwin} class object; list of grids.} 
#' \item{\code{size}: }{\code{textfield} class object; the value of the font size frame.} 
#' \item{\code{family}: }{\code{listbox} class object; the value of the font family frame.} 
#' \item{\code{colour}: }{\code{listbox} class object; the value of the colour set frame.} 
#' \item{\code{goption}: }{\code{checkboxes} class object; values of other options frame.} 
#' \item{\code{theme}: }{\code{radioboxes} class object; the value of the graph theme frame.} 
#' }
#' @section Contains:
#' \code{gparts_base}
#' @section Methods:
#' \describe{
#' \item{\code{back(perline = 3)}: }{\code{back} method for \code{gparts_base} class.}
#' \item{\code{front(top,
#'   showcolourbox = TRUE, 
#'   fontSize = unlist(options("kmg2FontSize")), 
#'   fontSize = unlist(options("kmg2FontSize")), 
#'   fontFamily = unlist(options("kmg2FontFamily")), 
#'   colourSet = unlist(options("kmg2ColourSet")), 
#'   saveGraph = unlist(options("kmg2SaveGraph")), 
#'   themeBase = unlist(options("kmg2Theme"))
#' )}: }{
#'   \code{front} method for \code{toolbox} subclass.
#' }
#' }
#' @family guiparts
#'
#' @name toolbox-class
#' @aliases toolbox
#' @rdname guiparts-toolbox
#' @docType class
#' @keywords hplot
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom grDevices windowsFonts postscriptFonts
#' @export toolbox
toolbox <- setRefClass(

  Class = "toolbox",

  fields = c("size", "family", "colour", "goption", "theme"),

  contains = c("gparts_base"),

  methods = list(

    front = function(
      top, showcolourbox = TRUE,
      fontSize   = unlist(options("kmg2FontSize")),
      fontFamily = unlist(options("kmg2FontFamily")),
      colourSet  = unlist(options("kmg2ColourSet")),
      saveGraph  = unlist(options("kmg2SaveGraph")),
      themeBase  = unlist(options("kmg2Theme"))
    ) {

      frame <<- tkframe(top)

      size <<- textfield$new()
      size$front(
        top          = frame,
        initialValue = fontSize,
        boxwidth     = "10", 
        title        = gettextKmg2("Font size")
      )
      
      if (.Platform$OS.type == "windows") {
        osfontsname <- names(windowsFonts())
      } else {
        osfontsname <- names(postscriptFonts())
      }

      family <<- variableListBox(
        parentWindow     = frame,
        variableList     = osfontsname,
        listHeight       = 5,
        selectmode       = "single",
        initialSelection = fontFamily,
        title            = gettextKmg2("Font family")
      )
      
      if (showcolourbox) {
        brewercolours <- rownames(brewer.pal.info)
        brewercolours <- c("Set1", brewercolours[brewercolours!="Set1"])
        colour <<- variableListBox(
          parentWindow     = frame,
          variableList     = brewercolours,
          listHeight       = 5,
          selectmode       = "single",
          initialSelection = colourSet,
          title            = gettextKmg2("Colour pattern")
        )
      }

      goption <<- checkboxes$new()
      goption$front(
        top        = frame,
        initValues = list(saveGraph),
        labels     = list(gettextKmg2("Save graph")),
        title      = gettextKmg2("Graph options")
      )
      
      theme <<- variableListBox(
        parentWindow     = frame,
        variableList     = c(
          "theme_bw", "theme_simple", "theme_classic",
          "theme_gray", "theme_minimal",
          "theme_linedraw", "theme_light", "theme_dark",
          "theme_base", "theme_calc", "theme_economist", 
          "theme_excel", "theme_few", "theme_fivethirtyeight", 
          "theme_gdocs", "theme_hc", "theme_par",
          "theme_pander", "theme_solarized", "theme_stata",
          "theme_tufte", "theme_wsj2", "theme_igray"
        ),
        listHeight       = 5,
        selectmode       = "single",
        initialSelection = themeBase,
        title            = gettextKmg2("Theme")
      )
      
      if (showcolourbox) {
        back_list <<- list(
          size$frame,
          family$frame,
          colour$frame,
          goption$frame,
          theme$frame
        )
      } else {
        back_list <<- list(
          size$frame,
          family$frame,
          goption$frame,
          theme$frame
        )
      }

      return()

    }

  )
)



#' The \code{front} Method for \code{toolbox} Subclass
#'
#' @usage \S4method{front}{toolbox}(top, 
#'   showcolourbox = TRUE, 
#'   fontSize = unlist(options("kmg2FontSize")), 
#'   fontSize = unlist(options("kmg2FontSize")), 
#'   fontFamily = unlist(options("kmg2FontFamily")), 
#'   colourSet = unlist(options("kmg2ColourSet")), 
#'   saveGraph = unlist(options("kmg2SaveGraph")), 
#'   themeBase = unlist(options("kmg2Theme"))
#' )
#' @param top \code{tkwin} class object; top of widget window.
#' @param showcolourbox Boolean; whether the colour set frame is shown or not.
#' @param fontSize Character; the initialization value of the font size.
#' @param fontFamily Numeric; the initialization value of the font family.
#' @param colourSet Numeric; the initialization value of the colour set.
#' @param saveGraph Numeric; the initialization value of the save graph option.
#' @param themeBase Numeric; the initialization value of the theme.
#' @family guiparts
#'
#' @name front,toolbox-method
#' @rdname guiparts-toolbox-front
#' @docType methods
#' @keywords hplot
NULL
