#guiSetStyle.tcltk <- function (style = c("classic", "SciViews", "Rcmdr", "large"),
#ask = FALSE)
#{
#### TODO: allow a more flexible definition of styles, with config files
#### TODO: if ask then show it in a dialog box
#    ## These four "hardcoded" styles are just here for demonstration
#    if (is.null(style)) style <- "classic"# Default value
#    style <- match.arg(style)
#    ## Look if .guiStyle.tcltk is set and in the correct style
#    if (existsTemp(".guiStyle.tcltk")) {  # Get it
#        st <- getTemp(".guiStyle.tcltk")
#        if (st$style == style) return(st) else {
#            ## Delete current fonts
#            tkfont.delete(st$font.label)
#            tkfont.delete(st$font.emph)
#            tkfont.delete(st$font.text)
#            tkfont.delete(st$font.fixed)
#            remove(".guiStyle.tcltk", pos = 1)
#        }
#    }
#    ## (Re)create fonts
#    st <- switch(style,
#        "classic"  = list(
#            style = "classic",
#            fg = c("black", "black", "black", "blue", "darkred", "gray12"),
#                 # Foreground normal, OK, cancel, selected, emphasized, disabled
#            pads = c(5, 5, 2), # padx, pady (wide), pady (small)
#            relief = "sunken",
#            sep = TRUE,
#            multiselect = "extended",
#            font.label = tkfont.create(family = "arial", size = 9),
#            font.emph  = tkfont.create(family = "arial", size = 9, weight = "bold"),
#            font.text  = tkfont.create(family = "times", size = 11),
#            font.fixed = tkfont.create(family = "courier", size = 9)),
#        "SciViews" = list(
#             style = "SciViews",
#             fg = c("black", "darkgreen", "darkred", "darkblue", "darkred", "gray12"),
#            pads = c(5, 5, 2),
#            relief = "flat",
#            sep = TRUE,
#            multiselect = "extended",
#            font.label = tkfont.create(family = "tahoma", size = 8),
#            font.emph  = tkfont.create(family = "tahoma", size = 8),
#            font.text  = tkfont.create(family = "verdana", size = 8),
#            font.fixed = tkfont.create(family = "lucida console", size = 9)),
#        "Rcmdr"  = list(
#            style = "Rcmdr",
#            fg = c("black", "darkgreen", "red", "blue", "red", "gray12"),
#            pads = c(0, 0, 0),
#            relief = "sunken",
#            sep = FALSE,
#            multiselect = "multiple",
#            font.label = tkfont.create(family = "arial", size = 9),
#            font.emph  = tkfont.create(family = "arial", size = 9),
#            font.text  = tkfont.create(family = "times", size = 11),
#            font.fixed = tkfont.create(family = "courier", size = 9)),
#        "large"    = list(
#            style = "large",
#            fg = c("black", "black", "black", "blue", "red", "gray12"),
#            pads = c(7, 7, 3),
#            relief = "sunken",
#            sep = TRUE,
#            multiselect = "extended",
#            font.label = tkfont.create(family = "arial", size = 14),
#            font.emph  = tkfont.create(family = "arial", size = 14, weight = "bold"),
#            font.text  = tkfont.create(family = "times", size = 16),
#            font.fixed = tkfont.create(family = "courier", size = 14)))
#    ## We need also the measure of these fonts to size widgets properly
#    st$font.measure <- c(
#        as.numeric(tkfont.measure(tkfont.actual(st$font.label), "0")),
#        as.numeric(tkfont.measure(tkfont.actual(st$font.emph), "0")),
#        as.numeric(tkfont.measure(tkfont.actual(st$font.text), "0")),
#        as.numeric(tkfont.measure(tkfont.actual(st$font.fixed), "0")))
#    names(st$font.measure) <- c("label", "emph", "text", "fixed")
#    assignTemp(".guiStyle.tcltk", st)
#    return(st)
#}
