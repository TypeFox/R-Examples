## These items still need to be implemented!
#dlgColor <- function (...)
#{
#    ## A color selection dialog box
#    ## TODO: a color range selector?
#    stop("Not yet implemented!")
#}

#dlgFont <- function (...)
#{
#    ## A font selector dialog box
#    ## TODO...
#    stop("Not yet implemented!")
#}

#dlgAssistant <- function (...)
#{
#    ## This is a non modal assistant dialog box... could also display tips
#    ## TODO...
#    stop("Not yet implemented!")
#}

#dlgDoubleList <- function (list1, list2, title = "Select", default1 = "",
#default2 = "", multi = c(TRUE, TRUE), new = c(FALSE, FALSE), sort = c(TRUE, TRUE),
#transfer = FALSE, parent = 0, GUI = getOption("guiWidgets"))
#{
#    ## A 'dual list' dialog box. This list serves two purposes:
#    ## 1) select elements in the first list and place them in the second list
#    ##    (transfer = TRUE)
#    ## 2) make a double selection in two separate lists
#    ## This dialog box is always modal
#    stop("Not yet implemented!")
#}

#dlgFormula <- function (...)
#{
#    ## This dialog box helps to create S language formulas
#    ## R Commander has something like that in glm dialog box. Look with John Fox
#    ## for permission to reuse it
#    ## TODO...
#    stop("Not yet implemented!")
#}

#dlgGraphOptions <- function (...)
#{
#    ## A graph options dialog box
#    ## Idem as guiDlgOptions, but specific to graph parameters? Or is it possible
#    ## to reuse guiDlgOptions?
#    ## TODO...
#    stop("Not yet implemented!")
#}

#dlgGrid <- function (table, title = deparse(substitute(table)), edit = TRUE,
#edit.vars = TRUE, add.vars = TRUE, add.rows = TRUE, parent = -1,
#GUI = getOption("guiWidgets"))
#{
#    ## A 'grid' display of a one or two dimensional table (vector, matrix,
#    ## data.frame)
#    ## It is similar to the data editor in Rgui, but can be non modal (parent = -1)
#    ## and can be used just to display the content of a table
#    ## TODO: possibly tabbed presentation?
#    stop("Not yet implemented!")
#}

#dlgItemSel <- function (list, classes = NULL, title = "Select items",
#default = "", default.items = "",  all.names = FALSE, multi = FALSE,
#sort = TRUE, sort.items = FALSE, subset = FALSE, save.restore = FALSE,
#parent = 0, GUI = getOption("guiWidgets"))
#{
#    ## Idem than guiDlgVarSel, but allows also to select items in lists,
#    ## data.frames, S4 objects, ...
#    ## This is the most advanced var selection dialog box
#    ## The argument subset allows to also subset rows
#    ## The argument save.restore allows to save (on a file, or in variables?
#    ## sel()...) and restore the given selected
#    ## TODO: display a summary of the selected object
#    stop("Not yet implemented!")
#}

#dlgOptions <- function (...)
#{
#    ## An option change dialog box
#    ## At left, options are presented in a tree (if possible)
#    ## corresponding to a list (or list of lists...)
#    ## and at right, the user can define labels, entry boxes, check boxes,
#    ## options, sliders, etc...
#    ## TODO...
#    stop("Not yet implemented!")
#}

#dlgProgress <- function (value, range = c(0, 100),
#message = "Please wait...", title = "Progression", percent = TRUE,
#cancel = TRUE, icon = "none", parent = -1, GUI = getOption("guiWidgets"))
#{
#    ## This is an (usually) non modal progression dialog box. It can also display
#    ## the text online (or use another function progress()?)
#    stop("Not yet implemented!")
#}

#dlgText <- function (text, file = NULL, title = deparse(substitute(text)),
#edit = TRUE, submit = TRUE, parent = -1, GUI = getOption("guiWidgets"))
#{
#    ## A simple 'text' display/editor. It can display the content of a character
#    ## variable, or of a file
#    ## It can be used as simple script editor, and if submit = TRUE, it offers
#    ## the possibility to submit its content to R
#    ## It is similar to the script editor in Rgui 2.0.0, and of its pager at the
#    ## same time
#    ## It can be non modal (parent = -1)
#    ## TODO: possibly tabbed presentation?
#    stop("Not yet implemented!")
#}

#dlgVarSel <- function (list, classes = NULL, title = "Select a variable",
#default = "", all.names = FALSE, multi = FALSE, new = FALSE, sort = TRUE,
#parent = 0, GUI = getOption("guiWidgets"))
#{
#    ## This is similar to guiDlgList(), but it is specific for variable selection.
#    ## If list is NULL, all variables of .GlobalEnv whose class belong to classes
#    ## are indicated
#    ## all.names indicates if hidden items should be also displayed
#    ## This dialog box is always modal
#    ## TO DO: display a summary of the selected object
#    stop("Not yet implemented!")
#}

#dlgView <- function (file, CSSfile, title = "View", report = TRUE,
#parent = -1, GUI = getOption("guiWidgets"))
#{
#    ## A 'view' dialog box (views are HTML presentation of R objects introduced
#    ## in SciViews-R version 0.6-0)
#    ## It can be non modal (parent = -1)
#    ## TODO: possibly tabbed presentation?
#    ## Possibly remove this, cf, SciViews specific!?
#    stop("Not yet implemented!")
#}
